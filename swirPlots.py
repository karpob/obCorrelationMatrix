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
        mask = '(ichan == {:d}) & ( ob > 0.0 ) '.format(subsetChan)
        d1.set_mask(mask)
        d1.use_mask(True)
        omf = d1.v('omf')
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
            plotMapHist(lat, lon, omf, 'Channel {:d} OMF for file {}'.format(i,os.path.basename(f)), os.path.basename(f)+'_Chan{:d}'.format(i), units='Kelvin') 
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
    instrumentChan['iasi'] = [ 5446 ,5455 ,5472 ,5480 ,5483 ,5485 ,5492 ,5497 ,5502 ,5507 ,5509 ,5517 ,5528 ,5558 ,5697 ,5714 ,5749 ,\
                               5766 ,5785 ,5798 ,5799 ,5801 ,5817 ,5833 ,5834 ,5836 ,5849 ,5851 ,5852 ,5865 ,5869 ,5881 ,5884 ,5897 ,5900 ,5916 ,\
                               5932 ,5948 ,5963 ,5968 ,5978 ,5988 ,5992 ,5994 ,5997 ,6003 ,6008 ,6023 ,6026 ,6039 ,6053 ,6056 ,6067 ,6071 ,6082 ,6085 ,\
                               6098 ,6112 ,6126 ,6135 ,6140 ,6149 ,6154 ,6158 ,6161 ,6168 ,6174 ,6182 ,6187 ,6205 ,6209 ,6213 ,6317 ,6339 ,6342 ,6366 ,6381 ,\
                               6391 ,6489 ,6962 ,6966 ,6970 ,6975 ,6977 ,6982 ,6985 ,6987 ,6989 ,6991 ,6993 ,6995 ,6997 ,6999 ,7000 ,7004 ,7008 ,7013 ,7016 ,7021 ,\
                               7024 ,7027 ,7029 ,7032 ,7038 ,7043 ,7046 ,7049 ,7069 ,7072 ,7076 ,7081 ,7084 ,7089 ,7099 ,7209 ,7222 ,7231 ,7235 ,7247 ,7267 ,7269 ,7284 ,\
                               7389 ,7419 ,7423 ,7424 ,7426 ,7428 ,7431 ,7436 ,7444 ,7475 ,7549 ,7584 ,7665 ,7666 ,7831 ,7836 ,7853 ,7865 ,7885 ,7888 ,7912 ,7950 ,7972 ,7980 ,\
                               7995 ,8007 ,8015 ,8055 ,8078 ]

    instrumentChan['airs'] = [1865 ,1866 ,1868 ,1869 ,1872 ,1873 ,1876 ,1881 ,1882 ,1883 ,1911 ,1917 ,1918 ,1924 ,1928 ,1937 ,1941 ,2099 ,2100 ,2101 ,2103 ,2104 ,2106 ,2107 ,\
                              2108 ,2109 ,2110 ,2111 ,2112 ,2113 ,2114 ,2115 ,2116 ,2117 ,2118 ,2119 ,2120 ,2121 ,2122 ,2123 ,2128 ,2134 ,2141 ,2145 ,2149 ,2153 ,2164 ,2189 ,2197 ,\
                              2209 ,2226 ,2234 ,2280 ,2318 ,2321 ,2325 ,2328 ,2333 ,2339 ,2348 ,2353 ,2355 ,2357 ,2363 ,2370 ,2371 ,2377]

    instrumentChan['cris_npp'] = [1147 ,1148 ,1149 ,1150 ,1151 ,1152 ,1153 ,1154 ,1155 ,1156 ,1157 ,1158 ,1159 ,1160 ,1161 ,1162 ,1163 ,\
                                  1164 ,1165 ,1166 ,1167 ,1168 ,1169 ,1170 ,1171 ,1172 ,1173 ,1174 ,1175 ,1177 ,1178 ,1179 ,1180 ,1181 ,1187 ,\
                                  1189 ,1190 ,1192 ,1193 ,1194 ,1196 ,1197 ,1198 ,1199 ,1200 ,1202 ,1203 ,1204 ,1206 ,1207 ,1208 ,1210 ,1212 ,1214 ,\
                                  1215 ,1217 ,1218 ,1220 ,1222 ,1224 ,1226 ,1228 ,1229 ,1231 ,1232 ,1234 ,1235 ,1236 ,1237 ,1238 ,1239 ,1241 ,1242 ,1243 ,\
                                  1244 ,1245 ,1247 ,1250 ,1270 ,1271 ,1282 ,1285 ,1288 ,1290 ,1293 ,1298 ,1301]

    instrumentChan['cris-fsr'] = [1596,1602 ,1619 ,1624 ,1635 ,1939 ,1940 ,1941 ,1942 ,1943 ,1944 ,1946 ,1947 ,1948 ,1949 ,1950 ,1951 ,\
                                  1952 ,1953 ,1954 ,1955 ,1956 ,1957 ,1958 ,1959 ,1960 ,1961 ,1962 ,1963 ,1964 ,1965 ,1966 ,1967 ,1968 ,\
                                  1969 ,1970 ,1971 ,1972 ,1973 ,1974 ,1975 ,1976 ,1977 ,1978 ,1979 ,1980 ,1981 ,1982 ,1983 ,1984 ,1985 ,\
                                  1986 ,1987 ,2119 ,2140 ,2143 ,2147 ,2153 ,2158 ,2161 ,2168 ,2171 ,2175 ,2182]

    ichans = instrumentChan[a.instrument]
    ichans.sort()
    
    # use given path and grab all nc diag files with instrument name in them.
    #files = glob.glob( os.path.join(a.path,'*'+a.instrument+'*.nc4') )

    files = getFiles(a.start, a.end, a.instrument, a.ops, a.experiment, a.diagtype, 'nc4')
    main(files, ichans, a.select, float(a.target), a.maps)    
