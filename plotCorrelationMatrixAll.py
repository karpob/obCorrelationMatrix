import numpy as np
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt 
import h5py 


def plotInstrument(instrument):
    h5 = h5py.File('etc/'+instrument+'_wavenumbers.h5','r')
    wavenumbersAll = np.asarray(h5['wavenumbers'])
    h5.close() 
    h5 = h5py.File('cris-fsr'+'.h5','r')
    corr = np.asarray(h5['correlationCombinedR'])
    ichan= np.asarray(h5['channels'])-1
    wavenumbers = wavenumbersAll[ichan]
    idxBufrSubset = ichan
    print("read instrument: {}".format(instrument))
    plt.matshow(corr)
    fig = plt.gcf()
    fig.set_size_inches(20,20)
    wavenumbersRounded = []
    for w in wavenumbers:
        wavenumbersRounded.append('{:10.3f}'.format(w))
    plt.xticks(np.arange(len(idxBufrSubset)),wavenumbersRounded,fontsize=18,rotation='vertical')
    plt.yticks(np.arange(len(idxBufrSubset)),wavenumbersRounded,fontsize=18)
    plt.clim(-1,1) 
    plt.colorbar()
    plt.savefig(instrument+'_all.png')
    plt.close()

    #do just SWIR that Erin assimilated
    idx, = np.where( (wavenumbers>=2380.) & (wavenumbers<=2510.0))
    cor1 = corr[idx,:]
    cor2 = cor1[:,idx]

    plt.matshow(cor2)
    fig = plt.gcf()
    fig.set_size_inches(20,20)
    plt.xticks(np.arange(len(idx)),np.asarray(wavenumbersRounded)[idx],fontsize=12,rotation='vertical')
    plt.yticks(np.arange(len(idx)),np.asarray(wavenumbersRounded)[idx],fontsize=12)
    plt.clim(-1,1) 
    plt.colorbar()
    plt.savefig(instrument+'_swir.png')
    plt.close()
    """
    uncomment if you want to look at error (and change from correlationCombinedR to covarianceCombinedR)
    plt.plot(np.arange(len(idx)),np.sqrt(np.diag(cor2)),'kx')
    fig = plt.gcf()
    fig.set_size_inches(20,20)
    plt.xticks(np.arange(len(idx)),np.asarray(wavenumbersRounded)[idx],fontsize=8,rotation='vertical')
    plt.savefig(instrument+'_swir_jones_error.png')
    for i,w in enumerate(np.asarray(wavenumbersRounded)[idx]):
        print(w,np.sqrt(np.diag(cor2)[i]))
    """
if __name__ == "__main__":
    plotInstrument('cris-fsr')
