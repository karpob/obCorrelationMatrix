#!/usr/bin/env python3
import argparse, os, sys, glob, h5py
import numpy as np
from lib.gsiCovarianceFile import gsiCovarianceFile
from matplotlib import pyplot as plt
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'read ncdiag files and create correlation matrix')
    parser.add_argument('--path', help = 'path to binary file to read in.', required = True, dest = 'path')
    parser.add_argument('--instrument',help = 'instrument name to process', required = True, dest='instrument')
    a = parser.parse_args()

    if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    if not os.path.exists(a.outpath): sys.exit("'{} outpath does not exist.".format(a.outpath))

    h5 = h5py.File('etc/'+a.instrument+'_wavenumbers.h5','r')
    idxBufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
    ichans = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    igeos = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    idxNucapsOzoneInInstrument = np.asarray(h5['idxNucapsOzoneInInstrument']).astype('int') 
    if (a.instrument == 'airs'):
        ibufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
        chans2add = [1012, 1024, 1088, 1111, 1120]
        chans2add_gsi =[]
        for c in chans2add:
            chans2add_gsi.append(np.where(c == ibufrSubset)[0][0]+1)
        new_variances = {}
        new_variances[chans2add_gsi[0]] = 1.7
        new_variances[chans2add_gsi[1]] = 1.0
        new_variances[chans2add_gsi[2]] = 1.0
        new_variances[chans2add_gsi[3]] = 1.4
        new_variances[chans2add_gsi[4]] = 1.4
    elif(a.instrument == 'iasi'):
        ibufrSubset = np.asarray(h5['idxBufrSubset']).astype('int')
        chans2add = [1427, 1479, 1536, 1579, 1585, 1626, 1643, 1671]
        chans2add_gsi =[]
        for c in chans2add:
            chans2add_gsi.append(np.where(c == ibufrSubset)[0][0]+1)
        new_variances = {}
        new_variances[chans2add_gsi[0]] = 1.6
        new_variances[chans2add_gsi[1]] = 1.4
        new_variances[chans2add_gsi[2]] = 1.6
        new_variances[chans2add_gsi[3]] = 1.5
        new_variances[chans2add_gsi[4]] = 1.4
        new_variances[chans2add_gsi[5]] = 1.4
        new_variances[chans2add_gsi[6]] = 1.7
        new_variances[chans2add_gsi[7]] = 1.6

    igeos_new = np.asarray(h5['geosAssimilated']).astype('int').tolist()
    for c in chans2add:
        if(c not in igeos):
            igeos_new.append(c)

    igeos_new.sort()

    igsi = []
    igsi_new = []
    for ig in igeos: 
        igsi.append(np.where(ig == ibufrSubset)[0][0]+1)

    for ig in igeos_new: 
        igsi_new.append(np.where(ig == ibufrSubset)[0][0]+1)
        if(ig in chans2add): print(ig,np.where(ig == ibufrSubset)[0][0]+1)
 
    #print('Reading binary file for use in the GSI.')
    gsi = gsiCovarianceFile( a.path )
    _, R = gsi.get()
    #R = 0.5*np.random.rand(len(igsi),len(igsi))
    igsi_positions = {}
    for i,c in enumerate(igsi):
        igsi_positions[c] = i
    newR = np.zeros([len(igsi_new),len(igsi_new)])
    for i,chn_i in enumerate(igsi_new):
        for j, chn_j in enumerate(igsi_new):
            if(chn_i in list(igsi_positions.keys()) and chn_j in list(igsi_positions.keys())):
                iold, jold = igsi_positions[chn_i],igsi_positions[chn_j] 
                newR[i,j] = R[iold,jold]
            elif i==j:
                newR[i,j] = new_variances[chn_i]  
    gsi.set(igsi_new, newR)
    gsi.setName( os.path.join(os.path.split(a.path)[0], os.path.split(a.path)[1]+'_ozone_added.bin'))
    gsi.write()
    plt.matshow(newR)
    plt.colorbar()
    plt.savefig(os.path.join(os.path.split(a.path)[0], a.instrument+'.png'), dpi=1200)
    
