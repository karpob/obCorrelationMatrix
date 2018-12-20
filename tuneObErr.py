#!/usr/bin/env python3
import argparse, os, sys, glob
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

def main(files, ichans, select_obs_omf_oma, target):
    for f in list(files):    
        obStatsList = processFile(f, ichans, select_obs_omf_oma, target) 

def processFile(f, ichans, select_obs_omf_oma, target):
    print( 'Stats for: {}'.format(f) ) 
    d1 = ncd.obs(f)
    sensor_chan = d1.v('sensor_chan')
    
    # go through each channel, filter it using ncdiag API, 
    print("Channel  mean(omf_uncorrected)    mean(omf_corrected)  cpen             RMS(omf)         STD(omf)        mean(sigo) ")
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d}) & (qcmark == 0)'.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        omf = d1.v('omf')
        omfnbc = d1.v('omfnbc')
        cpen = d1.stat('cpen','omf')
        sigoDesired = np.sqrt(np.mean((omf**2))/target)
        print( "{:d}     {:10.7f}               {:10.7f}           {:10.7f}       {:10.7f}       {:10.7f}     {:10.7f}  {:10.7f}".format(\
                  i,    omfnbc.mean(),          omf.mean(),        cpen,          np.sqrt( np.mean( (omf)**2 ) ) , np.std(omf), d1.v('sigo').mean(), sigoDesired ) )
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'read ncdiag files and create correlation matrix')
    parser.add_argument('--path', help = 'path to ncdiag', required = True, dest = 'path')
    parser.add_argument('--instrument',help = 'instrument name to process', required = True, dest='instrument')
    parser.add_argument('--select', help="select obs omf oma", dest='select', required=False, default='omf' )
    parser.add_argument('--target', help="select obs omf oma", dest='target', required=False, default='0.1' )
    a = parser.parse_args()

    #if a.instrument not in list(sensorOzoneChannelSets.keys()): sys.exit("'{}' instrument unknown".format(a.instrument))
    if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    
    instrumentChan = {}
    instrumentChan['iasi'] = [1427, 1479, 1536, 1579, 1585, 1626, 1643, 1671]
    instrumentChan['airs'] = [1012, 1024, 1088, 1111, 1120, 1669] 
    instrumentChan['cris'] = [577, 607, 626, 650, 667, 945, 991, 994]
    instrumentChan['cris-fsr'] = [596, 626, 646, 659]  

    ichans = instrumentChan[a.instrument]
    ichans.sort()
    
    # use given path and grab all nc diag files with instrument name in them.
    files = glob.glob( os.path.join(a.path,'*'+a.instrument+'*.nc4') )

    main(files, ichans, a.select,float(a.target))    
