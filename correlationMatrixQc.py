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
    for k in list(obStats.keys()): obStats[k] = np.asarray( obStats[k] )    

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
    #plt.matshow(correlationCombined)
    #plt.colorbar()
    #plt.show()
    #import pickle
    #pickle.dump( correlationCombined, open( "save2multi.p", "wb" ) )

def covarianceToCorrelation(covariance):
    """
    take covariance matrix and return correlation matrix.
    """
    d = np.sqrt(np.diag(covariance))
    dinv = np.matrix(np.linalg.inv(np.diag(d)))
    correlation =  dinv*np.matrix(covariance)*dinv
    return correlation

def getIdxAllChansPassQc(d, ichans, sensor_chan ):
    """
    Iterate through all channels, and get the subset of obs where all channels pass QC.
    """
    qcVals = {}
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d}) & (Water_Fraction == 1.0)  '.format(subsetChan)
        d.set_mask(mask)
        d.use_mask(True)
        qcVals[i] = d.v('QC_Flag')

    first = True
    for i in list(qcVals.keys()):
        if(first):
            combinedQc = qcVals[i]
            first = False
        else: combinedQc += qcVals[i]

    allGoodObsIdx, = np.where( combinedQc == 0)

    return allGoodObsIdx 

def listOfDictsToDictOfLists(listOfDicts):
    dictOfLists = {}
    for k in list(listOfDicts[0].keys()): dictOfLists[k]=[]
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
    idxUse = getIdxAllChansPassQc(d1, ichans, sensor_chan)
    obList = []
    #Read in obs, and append only locations where all channels pass QC
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d}) & (Water_Fraction == 1.0)'.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        #obList.append(d1.v('obs')[idxUse])
        obList.append(d1.v('Obs_Minus_Forecast_unadjusted')[idxUse])
    # make list into array, and make first dimension ob location, second dimension channel.
    obArray = np.asarray(obList).T
    # do stats on file f.
    obStats = {}
    obStats['count'] = obArray.shape[0]
    obStats['covariance'] = np.cov(obArray,rowvar=False) 
    obStats['mean'] = obArray.mean(axis=0)
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

    sensorOzoneChannelSets = {}
    sensorOzoneChannelSets['cris-fsr'] = [ 592, 594, 596, 598, 600, 602, 604, 611, 614, 616, 618, 620, 622, 626, 638, 646, 648, 652, 659]
    sensorOzoneChannelSets['cris'] =  [556, 557, 558, 560, 561, 562, 564, 565, 566, 569, 573, 574, 577, 580, 581, 584, 585, 587, 590, 591, 594,\
                                       597, 598, 601, 604, 607, 611, 614, 616, 617, 619, 622, 626, 628, 634, 637, 638, 640, 641, 642, 644, 646,\
                                       647, 650, 651, 652, 654, 655, 657, 659, 663, 667, 670] 
    sensorOzoneChannelSets['iasi']  = [1409, 1414, 1420, 1424, 1427, 1430, 1434, 1440, 1442, 1445, 1450, 1454, 1460, 1463, 1469, 1474, 1479, 1483, 1487, 1494,\
                                      1496, 1502, 1505, 1510, 1513, 1518, 1521, 1526, 1529, 1532, 1537, 1541, 1545, 1548, 1553, 1560, 1568, 1574, 1579, 1583,\
                                      1587, 1606, 1626, 1639, 1643, 1652, 1659, 1666, 1671, 1675, 1681, 1694, 1697]
    sensorOzoneChannelSets['airs'] =[ 1003, 1012, 1019, 1024, 1030, 1038, 1069, 1115, 1116, 1119, 1123 , 1130]

    parser = argparse.ArgumentParser( description = 'read ncdiag files and create correlation matrix')
    parser.add_argument('--path', help = 'path to ncdiag', required = True, dest = 'path')
    parser.add_argument('--outpath', help = 'path to ncdiag', required = True, dest = 'outpath')
    parser.add_argument('--instrument',help = 'instrument name to process', required = True, dest='instrument')
    parser.add_argument('--all', help="use all channels in instrument",action="store_true",required=False)
    parser.add_argument('--nthreads', help="number of threads", dest='nthreads', type=int, required=False, default=2 )
    a = parser.parse_args()

    if a.instrument not in list(sensorOzoneChannelSets.keys()): sys.exit("'{}' instrument unknown".format(a.instrument))
    if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    if not os.path.exists(a.outpath): sys.exit("'{} outpath does not exist.".format(a.path))

    # use given path and grab all nc diag files with instrument name in them.
    files = glob.glob( os.path.join(a.path,'*'+a.instrument+'*.nc4') )
    if a.all: ichans = []
    else: ichans = sensorOzoneChannelSets[a.instrument]   

    main(files, ichans, a.nthreads, a.outpath, a.instrument)    
