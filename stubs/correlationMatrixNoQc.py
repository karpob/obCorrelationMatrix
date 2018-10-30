#!/usr/bin/env python3
import argparse
# these should be part of the anaconda module.  
# Right now (23 Oct 2018), I use this module:
#    other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
# which can be loaded by executing:
#    module load other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
import numpy as np
import matplotlib.pyplot as plt

# these tools are mine.  They are available from:
#    https://github.com/will-mccarty/das_tools
# which also is in 
#   /discover/nobackup/wrmccart/das_tools/

# gmao_tools.py is just a random set of tools to handle some dates & 
#    whatnot.  It isn't a good name.
# ncdiag.py is an interface to handle the netcdf4 obs diag files 
import gmao_tools as gt
import ncdiag as ncd
import pandas as pd 

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'read ncdiag file and create correlation matrix')
    parser.add_argument('--ncfile', help = 'ncdiag', required = True, dest = 'ncfile')
    a = parser.parse_args()
    d = ncd.obs(a.ncfile)
    df = pd.DataFrame()
    ichans = [ 592, 594, 596, 598, 600, 602, 604, 611, 614, 616, 618, 620, 622, 626, 638, 646, 648, 652, 659]
    firstPass = True
    df = pd.DataFrame()
    sensor_chan = d.v('sensor_chan')
    for i in list(ichans):
        subsetChan, = np.where(sensor_chan == i)[0] + 1
        print(subsetChan)
        mask = '(ichan == {:d})'.format(subsetChan)
        d.set_mask(mask)
        d.use_mask(True)
        if(firstPass):
            df = pd.DataFrame({'{:d}'.format(i):d.v('obs')})
            firstPass = False
        else:
            obs = d.v('obs')
            df['{:d}'.format(i)] = pd.Series(d.v('obs'), index = df.index)
    #print(df.corr())
    plt.matshow(df.corr())
    plt.show()
