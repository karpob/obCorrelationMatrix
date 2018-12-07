import numpy as np
import matplotlib
matplotlib.use('Agg')   
matplotlib.rcParams.update({'font.size': 32}) 

 
from matplotlib import pyplot as plt 
import h5py 


def plotInstrument(instrument):
    h5 = h5py.File('etc/'+instrument+'_wavenumbers.h5','r')
    idxBufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    ichans = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    geos = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    print('geos',geos)
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
    
    print('additional', keepEm)
    wavenumbers = np.asarray(h5['wavenumbers'])
    h5 = h5py.File(instrument+'.h5','r')
    corr = np.asarray(h5['correlationCombined'])
    corr = np.delete(corr,dropEm, axis=0)
    corr = np.delete(corr,dropEm, axis=1)

    print("read instrument: {}".format(instrument))
    #print("corr shape {}, geosAssimilated {}, idxBufrSubset {}".format(corr.shape,len(geosAssimilated),len(idxBufrSubset), len(geosAssimilated) ) )

    listToKeep = [] 
    geosAssimilated = np.asarray(keepEm)-1.0 
    plt.matshow(corr)
    fig = plt.gcf()
    fig.set_size_inches(20,20)
    wavenumbersRounded = []
    for w in wavenumbers[geosAssimilated.astype(int)]:
        wavenumbersRounded.append('{:10.3f}'.format(w))
    plt.xticks(np.arange(len(geosAssimilated)),wavenumbersRounded,fontsize=32,rotation='vertical')
    plt.yticks(np.arange(len(geosAssimilated)),wavenumbersRounded,fontsize=32)
    plt.colorbar()
    print(corr.shape,len(wavenumbersRounded),len(geosAssimilated))
    plt.savefig(instrument+'_ozone_only.png')
if __name__ == "__main__":
    plotInstrument('airs')
    plotInstrument('cris')
    plotInstrument('iasi')
    plotInstrument('cris-fsr')
