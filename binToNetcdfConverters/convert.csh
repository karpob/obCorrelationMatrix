#!/bin/csh -f

set file = $argv[1]

unset argv
setenv argv

source /discover/nobackup/wrmccart/progress_cvs/./EnADAS-5_17_0p9/src/g5_modules

/discover/nobackup/wrmccart/progress_cvs/./EnADAS-5_17_0p9/Linux/bin/gsidiag_rad_bin2nc4.x -npred 12  $file
