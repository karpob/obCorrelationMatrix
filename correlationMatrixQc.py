#!/usr/bin/env python3
import argparse, os, sys, glob, h5py
from multiprocessing import Pool
from functools import partial

#    module load other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
import numpy as np
import matplotlib.pyplot as plt
# these tools are Will's.  They are available from:
#    https://github.com/will-mccarty/das_tools
# which also is in 
#   /discover/nobackup/wrmccart/das_tools/

# gmao_tools.py is just a random set of tools to handle some dates & 
#    whatnot.  It isn't a good name.
# ncdiag.py is an interface to handle the netcdf4 obs diag files 
import gmao_tools as gt
import ncdiag as ncd

def main(files, ichans, nThreads, outpath, instrument):
    # loop through files and process.

    print("Initialize Pool with {} Workers".format(nThreads))
    p = Pool(nThreads)
    # pass ichans in as non-iterable via partial.
    print("initialized.")
    print( "Processing {} files".format( len(files) ) )
    #obStatsList =  processFile(files[0],ichans)
    obStatsList = p.map( partial(processFile, ichans = ichans), files)
    print("Swapping out dictionaries for arrays.")
    # convert list of dictionaries to dictionary of numpy arrays
    obStats = listOfDictsToDictOfLists(obStatsList)        
    # convert list to arrays
    print('obstat keys', obStats.keys() )
    for k in list(obStats.keys()): 
        obStats[k] = np.asarray( obStats[k] )

    print("Computing overall covariance.") 
    covarianceCombined = total_covariance(obStats['covariance'], obStats['mean'], obStats['count'])
    print("Done overall covariance.")
    overallMean = obStats['mean'].mean()
    observationCount = obStats['count'].sum()
    print("Computing Correlation.")
    correlationCombined = covarianceToCorrelation( covarianceCombined )

    print("Writing File.")
    with h5py.File( os.path.join(outpath, instrument+'.h5'),"w" ) as f:
        dset = f.create_dataset("correlationCombined",data = correlationCombined)
        dset = f.create_dataset("covarianceCombined",data = covarianceCombined)
        dset = f.create_dataset("overallMean",data = overallMean)
        dset = f.create_dataset("observationCount",data = observationCount)
        dset = f.create_dataset("channels",data = ichans)
    print("Done.")

def covarianceToCorrelation(covariance):
    """
    take covariance matrix and return correlation matrix.
    """
    d = np.sqrt(np.diag(covariance))
    dinv = np.matrix(np.linalg.inv(np.diag(d)))
    correlation =  dinv*np.matrix(covariance)*dinv
    return correlation

def listOfDictsToDictOfLists(listOfDicts):
    dictOfLists = {}
    keys = []
    for i in range(len(listOfDicts)):
        for k in list(listOfDicts[i].keys()):
            if k not in keys:
                keys.append(k)

    for k in keys: dictOfLists[k]=[]
    for item in listOfDicts:
        for k in list(item.keys()):
            dictOfLists[k].append(item[k])
    return dictOfLists


def processFile(f, ichans):
    print(f) 
    d1 = ncd.obs(f)
    sensor_chan = d1.v('sensor_chan')
    
    # if ichans is an empty list, use alll sensor_chans
    if(len(ichans) == 0): ichans = sensor_chan
    print(ichans)
    # get locations where all channels pass QC.
    ombList = []
    obList = []
    oList = []
    qcList = []
    waterList = [] 
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d})'.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        obList.append(d1.v('obs'))
        ombList.append(d1.v('Obs_Minus_Forecast_unadjusted'))
        waterList.append(d1.v('Water_Fraction'))
        qcList.append(d1.v('qcmark'))
    
    # make list into array, and make first dimension ob location, second dimension channel.
    ombArray = np.asarray(ombList).T
    obArray = np.asarray(obList).T
    waterArray = np.asarray(waterList).T
    qcArray = np.asarray(qcList).T
    badLocs = []

    # filtering manually in python by location. Any channel at a location bad, drop all channels at that observation point.
    for i in range(ombArray.shape[0]):
        if( any(waterArray[i,:] == 1.0) or any(qcArray[i,:] != 0) or any(obArray[i,:] <= 0) ):
            badLocs.append(i) 
    ombArray = np.delete(ombArray, badLocs, axis=0)                
                
    # do stats on file f.
    obStats = {}
    if(ombArray.shape[0] < 2):
        return obStats
    else:
        obStats['count'] = ombArray.shape[0]
        obStats['covariance'] = np.cov(ombArray,rowvar=False) 
        obStats['mean'] = ombArray.mean(axis=0)
        return obStats

def total_covariance(covs, means, N):
    """
    Stolen from https://github.com/serega/statistika/blob/master/TotalCovariance.ipynb
    But changed so ram footprint doesn't break on discover.
    """
    c_shape = covs.shape
    grand_mean = np.dot(means.T, N) / (np.sum(N))
    ess = np.sum((covs.reshape((covs.shape[0], -1)) * (N - 1).reshape(N.shape[0], 1)).reshape(c_shape), axis=0)
    print("N, means shape, grand mean shape",N, means.shape,grand_mean.shape)
    # too much ram with this bit, use a running sum instead!
    #tgss = np.sum([ np.outer(x, x) * n for x, n  in zip(means - grand_mean, N)], axis=0)
    tgss = np.zeros([grand_mean.shape[0],grand_mean.shape[0]])
    for i,count in enumerate(N):
        # this weirdness is because we're getting strange numpy behavior on Discover if we multiply scalars by arrays.
        a = np.outer(means[i] - grand_mean, means[i] - grand_mean)
        a = a*count
        tgss += a
    tot = np.sum(N) - 1
    val = ess + tgss 
    val = val/tot
    return val

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'read ncdiag files and create correlation matrix')
    parser.add_argument('--path', help = 'path to ncdiag', required = True, dest = 'path')
    parser.add_argument('--outpath', help = 'path to ncdiag', required = True, dest = 'outpath')
    parser.add_argument('--instrument',help = 'instrument name to process', required = True, dest='instrument')
    parser.add_argument('--all', help="use all channels in instrument",action="store_true",required=False)
    parser.add_argument('--nthreads', help="number of threads", dest='nthreads', type=int, required=False, default=2 )
    a = parser.parse_args()

    #if a.instrument not in list(sensorOzoneChannelSets.keys()): sys.exit("'{}' instrument unknown".format(a.instrument))
    if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    if not os.path.exists(a.outpath): sys.exit("'{} outpath does not exist.".format(a.path))

    h5 = h5py.File('etc/'+a.instrument+'_wavenumbers.h5','r')
    idxBufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    ichans = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    idxNucapsOzoneInInstrument = np.asarray(h5['idxNucapsOzoneInInstrument']).astype('int')
    for i in idxBufrSubset:
        if i >= min(idxNucapsOzoneInInstrument) and i <= max(idxNucapsOzoneInInstrument):
            ichans.append(i)

    ichans.sort()

    # use given path and grab all nc diag files with instrument name in them.
    files = glob.glob( os.path.join(a.path,'*'+a.instrument+'*.nc4') )

    main(files, ichans, a.nthreads, a.outpath, a.instrument)    
