import numpy as np
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt 
import h5py 


def plotInstrument(instrument):
    h5 = h5py.File('etc/'+instrument+'_wavenumbers.h5','r')
    idxBufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    geosAssimilated = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    idxNucapsOzoneInInstrument = np.asarray(h5['idxNucapsOzoneInInstrument']).astype('int')
    wavenumbers = np.asarray(h5['wavenumbers'])
    h5 = h5py.File(instrument+'.h5')
    corr = np.asarray(h5['correlationCombined'])
    print("read instrument: {}".format(instrument))
    print("corr shape {}, geosAssimilated {}, idxBufrSubset {}".format(corr.shape,len(geosAssimilated),len(idxBufrSubset), len(geosAssimilated) ) )
    listToRemove = []
    print('assimilatedInGeos',geosAssimilated)
    print('BufrSubset',idxBufrSubset)
    for i in idxNucapsOzoneInInstrument:
        if i in idxBufrSubset:
            print('adding',i)
            geosAssimilated.append(i)

    geosAssimilated.sort()
    listToKeep = [] 
    for n,i in enumerate(idxBufrSubset):
        print(n,i)
        if i not in geosAssimilated:
            
            listToRemove.append(n)
        else:
            listToKeep.append(i)
    for ll in geosAssimilated:
        if ll not in listToKeep:
            print('what?',ll)
    #print("listToKeep", listToKeep)
    #print("Off",len(listToRemove) )
    #print("On",geosAssimilated)
    listToRemove.sort()
    corr = np.delete(corr,listToRemove, axis=0)
    corr = np.delete(corr,listToRemove, axis=1)
    geosAssimilated = np.asarray(geosAssimilated)-1
    plt.matshow(corr)
    fig = plt.gcf()
    fig.set_size_inches(20,20)
    wavenumbersRounded = []
    for w in wavenumbers[geosAssimilated.astype(int)]:
        wavenumbersRounded.append('{:10.3f}'.format(w))
    plt.xticks(np.arange(len(geosAssimilated)),wavenumbersRounded,fontsize=8,rotation='vertical')
    plt.yticks(np.arange(len(geosAssimilated)),wavenumbersRounded,fontsize=8)
    plt.colorbar() 
    print(corr.shape,len(wavenumbersRounded),len(geosAssimilated))
    plt.savefig(instrument+'.png')
if __name__ == "__main__":
    #plotInstrument('airs')
    plotInstrument('cris')
    #plotInstrument('iasi')
    plotInstrument('cris-fsr')
