import numpy as np
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt 
from matplotlib import cm 
import h5py 

def doPCA(corr):
    """
    For a given covariance/correlation matrix compute eigenvalues & eigenvectors, 
    sort compute % of variance explained, and cumulative % of variance explained.

    In: 
        corr: covariance/correlation matrix
    Out: 
        pairsOfValsAndVects: list of eigenvalues/eigenvectors
        explainedVariance: list of explained variance for each increasing eigenvalue
        cumulativeExplained: list of cumulative variance explained (starting with largest/first eigenvalue). 
    """
    vals, vecs = np.linalg.eigh(corr)

    pairsOfValsAndVecs = [(np.abs(vals[i]), vecs[:,i]) for i in range(len(vals))]

    # sort by increasing eigenvalue
    pairsOfValsAndVecs.sort()

    # make largest eigenvalue first
    pairsOfValsAndVecs.reverse()
    totalVariance = sum(vals)
    explainedVariance = []
    for p in pairsOfValsAndVecs:
        explainedVariance.append( (p[0]/totalVariance) * 100 )
    cumulativeExplained = np.cumsum( np.asarray(explainedVariance) )

    return pairsOfValsAndVecs, explainedVariance, cumulativeExplained

def getCorrelationMatrixOzoneSubset( instrument ):

    h5 = h5py.File('etc/'+instrument+'_wavenumbers.h5','r')
    idxBufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    ichans = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    geos = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    wavenumbers = np.asarray(h5['wavenumbers'])
    idxNucapsOzoneInInstrument = np.asarray(h5['idxNucapsOzoneInInstrument']).astype('int')
    
    keepEm = []
    for i in idxBufrSubset:
        if i >= min(idxNucapsOzoneInInstrument) and i <= max(idxNucapsOzoneInInstrument): 
            ichans.append(i)

    ichans.sort()
    dropEm = []
    for c,ii in enumerate(ichans):
        if(ii in geos):
            dropEm.append(c) 
        else:
            keepEm.append(ii)    
    
    h5 = h5py.File(instrument+'.h5','r')
    corr = np.asarray(h5['covarianceCombined'])
    corr = np.delete(corr,dropEm, axis=0)
    corr = np.delete(corr,dropEm, axis=1)

    print("read instrument: {}".format(instrument))

    listToKeep = [] 
    geosAssimilated = np.asarray(keepEm)-1.0 
    wavenumbersRounded = []
    for w in wavenumbers[geosAssimilated.astype(int)]:
        wavenumbersRounded.append('{:10.3f}'.format(w))

    return corr, wavenumbersRounded

def computeInstrument(instrument):

    corr, wavenumbersRounded = getCorrelationMatrixOzoneSubset( instrument )
    pairsOfValsAndVecs, explainedVariance, cumulativeExplained = doPCA( corr )
    eigenVals = [] 
    eigenVecs = [] 
    for i,p in enumerate(pairsOfValsAndVecs):
        plt.figure()
        plt.title(instrument.upper().replace('CRIS','CrIS')\
                  +' Eigenfuction {:d} Explains {:.4f}% of Variance, Cumulative {:.4f}%'.format(i+1, explainedVariance[i], cumulativeExplained[i]),\
                  fontsize = 18)
        plt.plot(p[1],'kx')
        fig = plt.gcf()
        fig.set_size_inches(5*4,4*4)
        plt.xticks(np.arange(len(wavenumbersRounded)),wavenumbersRounded,rotation='vertical',fontsize=16)
        plt.ylabel('Eigenvector')
        plt.xlabel('Wavenumber [cm$^{-1}$]')
        plt.savefig(instrument+'PC{:04d}.png'.format(i+1))
        plt.close()
        eigenVals.append(p[0])
        if(cumulativeExplained[i] < 95.0 or len(eigenVecs) < 5):
            eigenVecs.append( p[1] )
    plt.figure()
    plt.title('Eigenvalues for {}'.format(instrument.upper().replace('CRIS','CrIS')))
    plt.plot(np.arange(1,len(eigenVals)+1), eigenVals)
    plt.ylabel('Eigenvalue')
    plt.xlabel('PCA')
    plt.savefig(instrument+'_eigenvals.png')
    plt.close()

    plt.figure()
    plt.matshow( np.asarray(eigenVecs),aspect=1, cmap = cm.coolwarm, vmax = np.asarray(eigenVecs).max(), vmin = -1.0*np.asarray(eigenVecs).max()  )
    plt.title(instrument.upper().replace('CRIS','CrIS')+' PCA Eigenvectors vs. Cumulative Explained Variance %')
    fig = plt.gcf()
    R =float(len(eigenVecs))/float(len(wavenumbersRounded))
    fig.set_size_inches(10,10*R+1.75)
    explainedVarianceRounded = []
    for c in cumulativeExplained:
        if(c < 95.0 or len(explainedVarianceRounded) < 5): explainedVarianceRounded.append('{:10.3f}'.format(c))
    plt.yticks(np.arange(len(explainedVarianceRounded)),explainedVarianceRounded,fontsize=8)
    plt.xticks(np.arange(len(wavenumbersRounded)),wavenumbersRounded,fontsize=8,rotation='vertical')
    plt.ylabel('Explained Variance [%]')
    plt.xlabel('Wavenumber [cm$^{-1}$]')
    ax = fig.gca()
    ax.xaxis.set_ticks_position('bottom')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('PCA loading/Eigenvector')
    plt.savefig(instrument+'_pca_loadings.png')

if __name__ == "__main__":
    computeInstrument('airs')
    computeInstrument('cris')
    computeInstrument('iasi')
    computeInstrument('cris-fsr')
