#!/usr/bin/env python3
import argparse
# these should be part of the anaconda module.  
# Right now (23 Oct 2018), I use this module:
#    other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
# which can be loaded by executing:
#    module load other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
import numpy as np
import matplotlib.pyplot as plt

# these tools are mine.  They are available from:
#    https://github.com/will-mccarty/das_tools
# which also is in 
#   /discover/nobackup/wrmccart/das_tools/

# gmao_tools.py is just a random set of tools to handle some dates & 
#    whatnot.  It isn't a good name.
# ncdiag.py is an interface to handle the netcdf4 obs diag files 
import gmao_tools as gt
import ncdiag as ncd
import pandas as pd
def covarianceToCorrelation(covariance):
    """
    take covariance matrix and return correlation matrix.
    """
    d = np.sqrt(np.diag(covariance))
    dinv = np.matrix(np.linalg.inv(np.diag(d)))
    correlation =  dinv*np.matrix(covariance)*dinv
    return correlation

def total_covariance(covs, means, N):
    """
    Stolen from https://github.com/serega/statistika/blob/master/TotalCovariance.ipynb
    """
    c_shape = covs.shape
    grand_mean = np.dot(means.T, N) / (np.sum(N))
    ess = np.sum((covs.reshape((covs.shape[0], -1)) * (N - 1).reshape(N.shape[0], 1)).reshape(c_shape), axis=0)
    tgss = np.sum([ np.outer(x, x) * n for x, n  in zip(means - grand_mean, N)], axis=0) 
    return (ess + tgss) / (np.sum(N) - 1)

 
def getIdxAllChansPassQc(d, ichans, sensor_chan ):
    """
    Iterate through all channels, and get the subset of obs where all channels pass QC.
    """
    qcVals = {}
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d})'.format(subsetChan)
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
if __name__ == "__main__":
    ichans = [ 592, 594, 596, 598, 600, 602, 604, 611, 614, 616, 618, 620, 622, 626, 638, 646, 648, 652, 659]

    #plan for tomorrow. Make this so it calculates running covariance matrix, then move on to handling multiple ncdiags for say a month....

    parser = argparse.ArgumentParser( description = 'read ncdiag file and create correlation matrix')
    parser.add_argument('--ncfile', help = 'ncdiag', required = True, dest = 'ncfile')
    a = parser.parse_args()
    files = ['x0035_Ana.diag_cris-fsr_n20_ges.20180626_00z_001.nc4', 'x0035_Ana.diag_cris-fsr_n20_ges.20180626_00z_002.nc4']
    
    obStats = {}
    obStats['count'] = []
    obStats['covariance'] = []
    obStats['mean'] = []
    # loop through files and process.
    for f in files:
        d1 = ncd.obs(f)
        sensor_chan = d1.v('sensor_chan')

        # get locations where all channels pass QC.
        idxUse = getIdxAllChansPassQc(d1, ichans, sensor_chan)
        obList = []
 
        for i in list(ichans):
            subsetChan, = np.where(sensor_chan == i)[0] + 1
            mask = '(ichan == {:d})'.format(subsetChan)
            d1.set_mask(mask)
            d1.use_mask(True)
            obList.append(d1.v('obs')[idxUse])
        obArray = np.asarray(obList).T
        obStats['count'].append(obArray.shape[0])
        obStats['covariance'].append( np.cov(obArray,rowvar=False) )
        obStats['mean'].append(obArray.mean(axis=0))
    # convert list to arrays
    for k in list(obStats.keys()):
        obStats[k] = np.asarray(obStats[k])    
    
    covarianceCombined = total_covariance(obStats['covariance'], obStats['mean'], obStats['count'])
    correlationCombined = covarianceToCorrelation( covarianceCombined )

    plt.matshow(correlationCombined)
    plt.colorbar()
    plt.show()
    import pickle
    pickle.dump( correlationCombined, open( "save2.p", "wb" ) )

