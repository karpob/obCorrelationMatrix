# code to come up with a correlation matrix based on observations.
you need to make sure to checkout the "feature/swir" branch, along with submodules submodules, do this by:
```shell
git clone https://github.com/karpob/obCorrelationMatrix.git
git checkout feature/swir
git submodule update --init --recursive
```

Main script of interest is probably desroziersR.py which computes the so-called R matrix (observation error covariance). This can be run on gsi netcdf "diag" files that are placed into the following structure, for example:
```
/discover/nobackup/projects/gmao/obsdev/bkarpowi/tmp_swir_diags/anl
```
and
```
/discover/nobackup/projects/gmao/obsdev/bkarpowi/tmp_swir_diags/ges
```
Where files like ```x41_swir.diag_cris-fsr_n20_anl.20191221_18z.nc4``` would go into the ```anl``` directory and  files like ```x41_swir.diag_cris-fsr_n20_ges.20191221_18z.nc4``` go into the ges direcotry. 

To run the script:
```source sourceMe.csh``` 
Then the usage for the script itself:
usage: desroziersR.py [-h] --path PATH --outpath OUTPATH --instrument
                      INSTRUMENT [--all] [--nthreads NTHREADS]
                      [--select SELECT]
desroziersR.py: error: the following arguments are required: --path, --outpath, --instrument

So for a given set of diag files as given in the directory structure above:
```
./desroziersR.py --path /discover/nobackup/projects/gmao/obsdev/bkarpowi/tmp_swir_diags --outpath $PWD --instrument cris-fsr --nthreads 2 
```
If you connect to a compute node and have more than 2 processors you can set ```--nthreads``` according to the number of processors you have.

Once the script is completed it will produce the following:

```
cris-fsr.h5 
cris-fsr.bin
```
The first file can be used to generate plots, the second cand be read in by the GSI.
Then to plot the correlation matrix:
```
./plotCorrelationMatrixAll.py
```
This will produce two pngs
```
cris-fsr_all.png
cris-fsr_swir.png
```
The first will be the entire correlation matrix of assimilated channels, the second will just be the correlation matrix of the SWIR assimilated channels.

