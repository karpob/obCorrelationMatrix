#!/usr/bin/env python3
import argparse, os, sys, glob, h5py
from multiprocessing import Pool
from functools import partial

#    module load other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
# these tools are Will's.  They are available from:
#    https://github.com/will-mccarty/das_tools
# which also is in 
#   /discover/nobackup/wrmccart/das_tools/

# gmao_tools.py is just a random set of tools to handle some dates & 
#    whatnot.  It isn't a good name.
# ncdiag.py is an interface to handle the netcdf4 obs diag files 
import gmao_tools as gt
import ncdiag as ncd
from lib.gsiCovarianceFile import gsiCovarianceFile

def main(filesAnl, filesGes, ichans, igeos, igsi, nThreads, outpath, instrument, select_obs_omf_oma):
    # loop through files and process.
    fgroups = []
    for i,f in enumerate(filesAnl):
        fgroups.append([filesGes[i], f])
        
    print("Initialize Pool with {} Workers".format(nThreads))
    p = Pool(nThreads)
    # pass ichans in as non-iterable via partial.
    print("Initialized Pool.")
    print( "Processing {} files.".format( len(fgroups) ) )
    obStatsList = p.map( partial(processFile, ichans = ichans, igeos= igeos, select_obs_omf_oma = select_obs_omf_oma), fgroups)
    print("Swapping out dictionaries for arrays.")
    # convert list of dictionaries to dictionary of numpy arrays
    obStats = listOfDictsToDictOfLists(obStatsList)        

    # convert list to arrays
    for k in list(obStats.keys()): 
        obStats[k] = np.asarray( obStats[k] )
    print(obStats.keys())
    print("Computing overall covariance.") 
    covarianceCombinedOmf, overallMeanOmf = total_covariance(obStats['covariance_omf'], obStats['mean_omf'], obStats['count'])
    covarianceCombinedOma, overallMeanOma = total_covariance(obStats['covariance_oma'], obStats['mean_oma'], obStats['count'])
    combinedR = combinedDesroziersR( obStats['desroziersR'], obStats['count'] ) 
    print("Done overall covariance.")
    observationCount = np.sum(obStats['count'])
    print("Computing Correlation.")
    correlationCombinedOmf = covarianceToCorrelation( covarianceCombinedOmf )
    correlationCombinedOma = covarianceToCorrelation( covarianceCombinedOma )
    correlationCombinedR = covarianceToCorrelation( combinedR )
   
    print('Writing File: {}'.format( os.path.join(outpath, instrument+'.h5') ) )
    with h5py.File( os.path.join(outpath, instrument+'.h5'),"w" ) as f:
        dset = f.create_dataset("correlationCombinedOmf",data = correlationCombinedOmf)
        dset = f.create_dataset("covarianceCombinedOmf",data = covarianceCombinedOmf)
        dset = f.create_dataset("correlationCombinedOma",data = correlationCombinedOma)
        dset = f.create_dataset("covarianceCombinedOma",data = covarianceCombinedOma)
        dset = f.create_dataset("correlationCombinedR",data = correlationCombinedR)
        dset = f.create_dataset("covarianceCombinedR",data = combinedR)
        dset = f.create_dataset("overallMeanOma",data = overallMeanOma)
        dset = f.create_dataset("overallMeanOmf",data = overallMeanOmf)
        dset = f.create_dataset("observationCount",data = observationCount)
        dset = f.create_dataset("channels",data = ichans)

    print('Writing binary file for use in the GSI.')
    # Last step (which I forgot), and is probably necessary in most cases)- Symmetrize Desroziers estimate of R.
    mR = np.asmatrix(combinedR)
    mR = 0.5*(mR+mR.T)
    gsi = gsiCovarianceFile( os.path.join(outpath, instrument+'.bin') )
    gsi.set( igsi, np.asarray(mR) )
    gsi.write()
    gsi.plot()
    print("Done!")

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


def processFile(fs, ichans, igeos, select_obs_omf_oma ):
    print( 'Processing Files {} {}'.format(fs[0],fs[1]) ) 
    # input is list of files 1st file ges, second file anl
    d1 = ncd.obs(fs[0])
    d2 = ncd.obs(fs[1])
    sensor_chan = d1.v('sensor_chan')
    # if ichans is an empty list, use alll sensor_chans
    if(len(ichans) == 0): ichans = sensor_chan
    # get locations where all channels pass QC.
    ombList = []
    omaList = []
    obList = []
    oList = []
    qcList = []
    waterList = []
    qcListA = []
    waterListA = []
    # create idx associated with geos assimilated channel
    idx_geos_assim = []
    for ii,c in enumerate(ichans):
        if (c in igeos):
            idx_geos_assim.append(ii)
    # go through each channel, filter it using ncdiag API, 
    # create a big list of observations (spatial dimension) for a given channel, 
    # add the list for a given channel to the bigger list so obList[channel][spatial]
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d})'.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        obList.append(d1.v('obs'))
        ombList.append(d1.v('omf'))
        waterList.append(d1.v('Water_Fraction'))
        qcList.append(d1.v('qcmark'))

    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d})'.format(subsetChan)
        d2.set_mask(mask)
        d2.use_mask(True)
        omaList.append(d2.v('omf'))
        waterListA.append(d2.v('Water_Fraction'))
        qcListA.append(d2.v('qcmark'))
 
        
    
    # make list into array, and make first dimension ob location, second dimension channel.
    ombArray = np.asarray(ombList).T
    omaArray = np.asarray(omaList).T
    obArray = np.asarray(obList).T
    waterArray = np.asarray(waterList).T
    qcArray = np.asarray(qcList).T

    waterArrayA = np.asarray(waterListA).T
    qcArrayA = np.asarray(qcListA).T

    badLocs = []

    # filtering manually in python by location. Any channel at a location bad, drop all channels at that observation point.
    for i in range(ombArray.shape[0]):
        if( any( waterArrayA[i,:] < 1.0 ) or any(qcArrayA[i,idx_geos_assim] != 0) or\
            any( waterArray[i,:] < 1.0 ) or any(qcArray[i,idx_geos_assim] != 0) or any(obArray[i,:] <= 0) ):
            badLocs.append(i) 
    ombArray = np.delete(ombArray, badLocs, axis=0)                
    omaArray = np.delete(omaArray, badLocs, axis=0)                
                
    # do stats on file f.
    obStats = {}
    if(ombArray.shape[0] < 1):
        return obStats
    else:
        obStats['count'] = ombArray.shape[0]
        obStats['covariance_omf'] = np.cov(ombArray,rowvar=False) 
        obStats['covariance_oma'] = np.cov(omaArray,rowvar=False) 
        obStats['mean_omf'] = ombArray.mean(axis=0)
        obStats['mean_oma'] = ombArray.mean(axis=0)
        obStats['desroziersR'] = np.zeros( [ombArray.shape[1], omaArray.shape[1]] )
        for i in range(ombArray.shape[1]):
            for j in range(omaArray.shape[1]):
                obStats['desroziersR'][i,j] = np.sum( ombArray[:,i]*omaArray[:,j] ) 
            
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
    return val, grand_mean
def combinedDesroziersR( R, count):
    """
    Combined R (where R is really just the sum, not covariance matrix (yet)) Sum the sums and divide by count to get the R.
    """
    
    combinedR = R.sum(axis=0)/count.sum() # maybe since this is a covariance count-1 to make statisticians happy?

    return combinedR

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'read ncdiag files and create correlation matrix')
    parser.add_argument('--path', help = 'path to ncdiag', required = True, dest = 'path')
    parser.add_argument('--outpath', help = 'path to ncdiag', required = True, dest = 'outpath')
    parser.add_argument('--instrument',help = 'instrument name to process', required = True, dest='instrument')
    parser.add_argument('--all', help="use all channels in instrument",action="store_true",required=False)
    parser.add_argument('--nthreads', help="number of threads", dest='nthreads', type=int, required=False, default=2 )
    parser.add_argument('--select', help="select obs omf oma", dest='select', required=False, default='omf' )
    a = parser.parse_args()

    #if a.instrument not in list(sensorOzoneChannelSets.keys()): sys.exit("'{}' instrument unknown".format(a.instrument))
    if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    if not os.path.exists(a.outpath): sys.exit("'{} outpath does not exist.".format(a.outpath))

    h5 = h5py.File('etc/'+a.instrument+'_wavenumbers.h5','r')
    idxBufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    ichans = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    igeos = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    idxNucapsOzoneInInstrument = np.asarray(h5['idxNucapsOzoneInInstrument']).astype('int') 
    ibufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    ozoneChans = {}
    ozoneChans['iasi'] = [ 1427, 1479, 1536, 1579, 1585, 1626, 1643, 1671 ]
    ozoneChans['airs'] = [ 1012, 1024, 1088, 1111, 1120]
    ozoneChans['cris-fsr'] = [ 596, 626, 646, 659 ] 
    ozoneChans['cris'] = [ 577, 607, 626, 650, 667 ]
    if(a.instrument in list(ozoneChans.keys())):
        for c in ozoneChans[a.instrument]:
            ichans.append(c)
            igeos.append(c)
           
    #for i in idxBufrSubset:
    #    if i >= min(idxNucapsOzoneInInstrument) and i <= max(idxNucapsOzoneInInstrument):
    #        ichans.append(i)
    ichans.sort()
    igeos.sort()
    # generate gsi style index for gsi R covariance file.
    igsi = []
    for ig in igeos: 
        igsi.append(np.where(ig == ibufrSubset)[0][0]+1)
 
    # use given path and grab all nc diag files with instrument name in them.
    filesAnl = glob.glob( os.path.join(a.path,'anl','*'+a.instrument+'*.nc4') )
    filesGes = glob.glob( os.path.join(a.path,'ges','*'+a.instrument+'*.nc4') )
    filesAnl.sort()
    filesGes.sort() 
    main(filesAnl,filesGes, ichans, igeos, igsi, a.nthreads, a.outpath, a.instrument, a.select)    
