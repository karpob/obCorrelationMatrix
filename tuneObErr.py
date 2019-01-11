#!/usr/bin/env python3
import argparse, os, sys, glob
import numpy as np
import matplotlib 
matplotlib.use('Agg')
from datetime import timedelta, date
import matplotlib.pyplot as plt
from lib.maps import plotMapHist 
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

def dateRange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)
 
def getFiles(start, end, instrument, opsPath, experimentName, anlOrGes, ncOrBin):
    files = []
    # pull out years, etc from args.
    startYear, startMonth, startDay, startHour = int(start[0:4]), int(start[4:6]), int(start[6:8]), int(start[8:10])
    endYear, endMonth, endDay, endHour = int(end[0:4]), int(end[4:6]), int(end[6:8]), int(end[8:10])

    # basepath
    pathInit =  os.path.join(opsPath, experimentName, 'obs')
    
    # get start/end date, put in datetime
    startDate = date(startYear, startMonth, startDay)
    endDate = date(endYear, endMonth, endDay)

    for today in dateRange(startDate, endDate):
        for hour in ['00','06','12','18']:
            path = os.path.join(pathInit, today.strftime("Y%Y/M%m/D%d"), 'H'+hour)
            if not os.path.exists(path): print(path +'does not exist.' )
            else:
                if( len( glob.glob(path+'/*'+instrument+'*'+anlOrGes+'*.'+ncOrBin) ) > 0): 
                    files.append(glob.glob(path+'/*'+instrument+'*'+anlOrGes+'*.'+ncOrBin)[0])
    return files            


def processFile(f, ichans, select_obs_omf_oma, target):
    print( 'Stats for: {}'.format(f) ) 
    d1 = ncd.obs(f)
    sensor_chan = d1.v('sensor_chan')
    
    # go through each channel, filter it using ncdiag API, 
    print("Channel  mean(omf_uncorrected)    mean(omf_corrected)  cpen             RMS(omf)         STD(omf)        mean(sigo) ")
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d}) & (qcmark == 0) & ( ob > 0.0 )'.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        omf = d1.v('omf')
        #oma = d1.v('oma')
        lat = d1.v('lat')
        lon = d1.v('lon')
        omfnbc = d1.v('omfnbc')
        cpen = d1.stat('cpen','omf')
        sigoDesired = np.sqrt(np.mean((omf**2))/target)
        plotMapHist(lat, lon, omf, 'Channel {:d} OMF for file {}'.format(i,os.path.basename(f)), os.path.basename(f)+'_Chan{:d}'.format(i), units='Kelvin') 
        print( "{:d}     {:10.7f}               {:10.7f}           {:10.7f}       {:10.7f}       {:10.7f}     {:10.7f}  {:10.7f}".format(\
                  i,    omfnbc.mean(),          omf.mean(),        cpen,          np.sqrt( np.mean( (omf)**2 ) ) , np.std(omf), d1.v('sigo').mean(), sigoDesired ) )
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'read ncdiag files and create correlation matrix')
    parser.add_argument('--instrument',help = 'instrument name to process', required = True, dest='instrument')
    parser.add_argument('--experiment', help = 'Experiment name in ops', required = True,dest = 'experiment')
    parser.add_argument('--start', help = 'start dtg YYYYMMDDhh', required = True, dest = 'start')
    parser.add_argument('--end', help = 'end dtg YYYYMMDDhh', required = True, dest = 'end')
    parser.add_argument('--ops', help = 'Optional arg to specify ops archive.', required = False, dest = 'ops',default="/archive/u/bkarpowi")
    parser.add_argument('--diagtype', help = 'specify to copy ges or anl.', required = False, dest = 'diagtype',default="anl")
    parser.add_argument('--select', help="select obs omf oma", dest='select', required=False, default='omf' )
    parser.add_argument('--target', help="select obs omf oma", dest='target', required=False, default='0.1' )
    a = parser.parse_args()

    #if a.instrument not in list(sensorOzoneChannelSets.keys()): sys.exit("'{}' instrument unknown".format(a.instrument))
    #if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    
    instrumentChan = {}
    instrumentChan['iasi'] = [1427, 1479, 1536, 1579, 1585, 1626, 1643, 1671]
    instrumentChan['airs'] = [1012, 1024, 1088, 1111, 1120, 1669] 
    instrumentChan['cris_npp'] = [577, 607, 626, 650, 667, 945, 991, 994]
    instrumentChan['cris-fsr'] = [596, 626, 646, 659]  

    ichans = instrumentChan[a.instrument]
    ichans.sort()
    
    # use given path and grab all nc diag files with instrument name in them.
    #files = glob.glob( os.path.join(a.path,'*'+a.instrument+'*.nc4') )

    files = getFiles(a.start, a.end, a.instrument, a.ops, a.experiment, a.diagtype, 'nc4')
    main(files, ichans, a.select, float(a.target))    
