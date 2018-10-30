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
    parser = argparse.ArgumentParser( description = 'read ncdiag file and create correlation matrix')
    parser.add_argument('--ncfile', help = 'ncdiag', required = True, dest = 'ncfile')
    a = parser.parse_args()
    d = ncd.obs(a.ncfile)
    df = pd.DataFrame()
    firstPass = True
    df = pd.DataFrame()
    sensor_chan = d.v('sensor_chan')

    ichans = [ 592, 594, 596, 598, 600, 602, 604, 611, 614, 616, 618, 620, 622, 626, 638, 646, 648, 652, 659]

    idxUse = getIdxAllChansPassQc(d, ichans, sensor_chan)
    Aa = [] 
    for i in list(ichans):
        print(i)
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d})'.format(subsetChan)
        d.set_mask(mask)
        d.use_mask(True)
        if(firstPass):
            df = pd.DataFrame({'{:d}'.format(i):d.v('obs')[idxUse]})
            Aa.append(d.v('obs')[idxUse])
            firstPass = False
        else:
            obs = d.v('obs')
            df['{:d}'.format(i)] = pd.Series(d.v('obs')[idxUse], index = df.index)
            Aa.append(d.v('obs')[idxUse])

    X = np.asarray(Aa).T
    totalObs = np.asarray(Aa).shape[1]
    
    N1 = np.array([407,407,407])
    print(N1)
    groups = [ X[0:N1[0],:], X[N1[0]:N1[0] + N1[1],:], X[N1[0]+N1[1]:,:]]
    covs = np.array([np.cov(g,rowvar=False) for g in groups])
    means = np.array([g.mean(axis=0) for g in groups])
    covCombined = total_covariance(covs, means, N1)
    d = np.sqrt(np.diag(covCombined))
    dinv = np.matrix(np.linalg.inv(np.diag(d)))
    corCombined =  dinv*np.matrix(covCombined)*dinv
    # print(c1-c2)
    #plt.matshow(c1-c2)
    plt.matshow(corCombined-df.corr())
    plt.colorbar()
    import pickle
    mycorr = df.corr()
    pickle.dump(mycorr , open( "save1.p", "wb" ) )

    #plt.matshow(covA)
    #plt.colorbar()

    #plt.matshow(covB)
    #plt.colorbar()

    plt.show()
    print(A.shape)
