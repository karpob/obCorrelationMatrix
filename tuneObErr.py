#!/usr/bin/env python3
import argparse, os, sys, glob
import numpy as np
import matplotlib 
matplotlib.use('Agg')
from datetime import timedelta, date, datetime
import matplotlib.pyplot as plt
from lib.graphics.maps import plotMapHist 
# these tools are Will's.  They are available from:
#    https://github.com/will-mccarty/das_tools
# which also is in 
#   /discover/nobackup/wrmccart/das_tools/

# gmao_tools.py is just a random set of tools to handle some dates & 
#    whatnot.  It isn't a good name.
# ncdiag.py is an interface to handle the netcdf4 obs diag files 
import gmao_tools as gt
import ncdiag as ncd
import pandas as pd

def main(files, ichans, select_obs_omf_oma, target, mapsOn):
    chanStatList = []
    experiment = os.path.basename(files[0]).split('.')[0]
    instrument = os.path.basename(files[0]).split('.')[1].split('_')[1]
    start = os.path.basename(files[0]).split('.')[2]
    end = os.path.basename(files[-1]).split('.')[2]
    # loop through files generate maps
    for f in list(files): 
        chanStatList.extend( processFile(f, ichans, select_obs_omf_oma, target, mapsOn) )
    chanStatArray = np.asarray(chanStatList)
    timeIndex = chanStatArray[:,0]
    df = pd.DataFrame(chanStatArray, index = timeIndex, columns= ['date', 'channel', 'Bias','OmF', 'Penalty', 'rmse', 'std', 'sigo','sigo_desired'])
    grouped = df.groupby('channel')
    f, ax = plt.subplots(nrows=len(ichans), ncols=1, figsize=(14,3*len(ichans)) )
    i = 0
    for n, g in grouped:
        g.plot(ax = ax[i], y=['Bias','OmF','Penalty'],legend=False)
        dd = g['date'].values
        omf = np.array(g['OmF'].values,dtype=np.float64)
        std = np.array(g['std'].values,dtype = np.float64)
        ax[i].fill_between(dd, omf+std, omf-std, alpha=0.3, facecolor='orange')
        ax[i].set_title('Channel {:d} Bias = {:.4f} Mean OmF = {:.4f} Std= {:.4f} Penalty= {:.4f} RMSE = {:.4f} $\sigma_o$= {:.4f} $\sigma_o(desired)$= {:4f}'.format(ichans[i], g['Bias'].mean(), g['OmF'].mean(), g['std'].mean(), g['Penalty'].mean(), g['rmse'].mean(), g['sigo'].mean(), g['sigo_desired'].mean() ) )
        ax[i].set_ylabel('Kelvin')
        ax[i].legend(bbox_to_anchor= (1.1,1),loc='upper right')
        i+=1
    plt.tight_layout()
    plt.savefig(experiment+'_'+instrument+'_'+start+'_'+end+'_stats_timeseries.png')
  
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
    print('pathInit',pathInit) 
    # get start/end date, put in datetime
    startDate = date(startYear, startMonth, startDay)
    endDate = date(endYear, endMonth, endDay)
    print(startDate, endDate)
    for today in dateRange(startDate, endDate):
        for hour in ['00','06','12','18']:
            path = os.path.join(pathInit, today.strftime("Y%Y/M%m/D%d"), 'H'+hour)
            print(path)
            if not os.path.exists(path): print(path +'does not exist.' )
            else:
                if( len( glob.glob(path+'/*'+instrument+'*'+anlOrGes+'*.'+ncOrBin) ) > 0): 
                    files.append(glob.glob(path+'/*'+instrument+'*'+anlOrGes+'*.'+ncOrBin)[0])
    return files            


def processFile(f, ichans, select_obs_omf_oma, target, mapsOn):
    print( 'Stats for: {}'.format(f) )
    dtg = os.path.basename(f.split('.')[-2])
    dateTag = datetime(int(dtg[0:4]), int(dtg[4:6]), int(dtg[6:8]), int(dtg[9:11]))
    d1 = ncd.obs(f)
    sensor_chan = d1.v('sensor_chan')
    chanStats = []
    # go through each channel, filter it using ncdiag API, 
    print("Channel  mean(omf_uncorrected)    mean(omf_corrected)  cpen             RMS(omf)         STD(omf)        mean(sigo) ")
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        mask = '(ichan == {:d}) & (qcmark == 0) & ( ob > 0.0 )'.format(subsetChan)
        #mask = '(ichan == {:d})'.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        omf = d1.v('omf')
        obs = d1.v('obs')
        #oma = d1.v('oma')
        lat = d1.v('lat')
        lon = d1.v('lon')
        omfnbc = d1.v('omfnbc')
        cpen = d1.stat('cpen','omf')
        sigoDesired = np.sqrt(np.mean((omf**2))/target)
        
        # Generate and save some stats.
        s ={}
        s['Channel'] = i 
        s['bias'] = omfnbc.mean()
        s['omf'] = omf.mean()
        s['penalty'] = cpen
        s['rmse'] = np.sqrt( np.mean( (omf)**2 ) )
        s['std'] = np.std(omf)
        s['sigo'] = d1.v('sigo').mean()
        s['sigo_desired'] = sigoDesired
        if(mapsOn):
            plotMapHist(lat, lon, omfnbc, 'Channel {:d} OMF for file {}'.format(i,os.path.basename(f)), os.path.basename(f)+'_Chan{:d}'.format(i), units='Kelvin') 
        print( "{:d}     {:10.7f}               {:10.7f}           {:10.7f}       {:10.7f}       {:10.7f}     {:10.7f}  {:10.7f}".format(\
                  i,    s['bias'],          s['omf'],        s['penalty'],         s['rmse'] , s['std'], s['sigo'], s['sigo_desired'] ) )
    
        chanStats.append([dateTag, i, s['bias'], s['omf'], s['penalty'], s['rmse'] , s['std'], s['sigo'], s['sigo_desired'] ])
    return chanStats
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
    parser.add_argument('--no-maps', help="turn of generating OmF maps/histograms.", dest='maps', action='store_false' )
    a = parser.parse_args()

    #if a.instrument not in list(sensorOzoneChannelSets.keys()): sys.exit("'{}' instrument unknown".format(a.instrument))
    #if not os.path.exists(a.path): sys.exit("'{} path does not exist.".format(a.path))
    
    instrumentChan = {}
    #instrumentChan['iasi'] = [1427, 1479, 1536, 1579, 1585, 1626, 1643, 1671]
    #instrumentChan['airs'] = [1012, 1024, 1088, 1111, 1120, 1669] 
    #instrumentChan['cris_npp'] = [266, 577, 607, 626, 650, 667, 945, 991, 994]
    """ Based on NOAA selection.
    instrumentChan['cris-fsr_n20'] = [  95,  99, 103, 107, 109, 127, 131, 125, 121, 123, 153,
                                      1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959] 
    """
    #Based on GMAO x0041 selection
    #instrumentChan['cris-fsr_n20'] = [  95,  99, 104, 107, 104, 126, 130, 126, 120, 123, 153,
    #                                  1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959] 
     
    instrumentChan['cris-fsr_n20'] = [  95, 104, 104, 130,
                                      1949,1951,1953,1955,1945,1947,2143] 
    #instrumentChan['cris-fsr_npp'] = [596, 626, 646, 659]  

    ichans = instrumentChan[a.instrument]
    ichans.sort()
    
    # use given path and grab all nc diag files with instrument name in them.
    #files = glob.glob( os.path.join(a.path,'*'+a.instrument+'*.nc4') )

    files = getFiles(a.start, a.end, a.instrument, a.ops, a.experiment, a.diagtype, 'nc4')
    main(files, ichans, a.select, float(a.target), a.maps)    
