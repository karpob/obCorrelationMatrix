#module load other/SSSO_Ana-PyD/SApd_4.2.0_py3.5
module use -a /home/mathomp4/modulefiles
module load python/anaconda/5.2.0/3.6
# bit to add Will's library
if ( ! ($?PYTHONPATH) ) then
    setenv PYTHONPATH $PWD/das_tools/python_tools
else
    setenv PYTHONPATH $PYTHONPATH\:$PWD/das_tools/python_tools
endif 
