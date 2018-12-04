#!/usr/bin/env python3

import argparse, os, glob, shutil 
from datetime import timedelta, date

def dateRange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)
 
def getFiles(start, end, instrument, opsPath, experimentName):
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
            print("path I'm looking at",path)
            if not os.path.exists(path): print(path +'does not exist.' )
            else: files.append(glob.glob(path+'/*'+instrument+'*ges*')[0])
    return files            

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'Script to pull, and convert NETCDF files from OPS. python traverseObs.py\
                                                     --experiment x0035_Ana --start 2018060100 --end 2018060300 --instrument cris_npp')
    parser.add_argument('--experiment', help = 'Experiment name in ops', required = True,dest = 'experiment')
    parser.add_argument('--start', help = 'start dtg YYYYMMDDhh', required = True, dest = 'start')
    parser.add_argument('--end', help = 'end dtg YYYYMMDDhh', required = True, dest = 'end')
    parser.add_argument('--instrument', help = 'instrument name', required = True, dest = 'instrument')
    parser.add_argument('--ops', help = 'Optional arg to specify ops archive.', required = False, dest = 'ops',default="/archive/u/dao_it/")
    parser.add_argument('--exe', help = 'Binary converter.', required = False, dest = 'exe',\
                        default="/discover/nobackup/wrmccart/progress_cvs/./EnADAS-5_17_0p9/Linux/bin/gsidiag_rad_bin2nc4.x")
    parser.add_argument('--outpath', help ='path to store output', required=True, dest ='outpath')
    a = parser.parse_args()

    destination = a.outpath
    fullPath = os.path.abspath(destination)
    files = getFiles(a.start, a.end, a.instrument, a.ops, a.experiment)
    print("User input: start, end, instrument, ops path, Experiment Name, Destination")
    print(a.start, a.end, a.instrument, a.ops, a.experiment, fullPath)
    #if destination doesn't exist create it.
    if ( not os.path.exists(destination) ):
        os.makedirs(destination)
        print( "Making directory: {}".format(fullPath) )

    here = os.path.dirname(os.path.abspath(__file__))  
    convertList = open( os.path.join( here, 'convert.list'), 'w' )
    for f in files:
        print('copying',f+' to '+ os.path.join( fullPath, os.path.split(f)[1] ) )

        # binary file copied to output path
        # copy the file
        shutil.copyfile(f, os.path.join(destination, os.path.split(f)[1] ) )
        convertList.write('{}/convert.csh {}\n'.format( here , os.path.abspath( os.path.join(destination, os.path.split(f)[1]) ) ) )
    convertList.close()
