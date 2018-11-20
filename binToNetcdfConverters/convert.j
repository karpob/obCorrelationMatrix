#!/bin/csh -fx
# ------------------------------
#SBATCH -A g0613
#SBATCH --export=NONE
#
#PBS -N cnvrt
#PBS -o cnvrt.log.o%j
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=4
#SBATCH --constraint=hasw
##SBATCH --qos=debug
##PBS -l walltime=1:00:00
##SBATCH --partition=preops
##SBATCH --qos=dastest
##SBATCH --qos=obsdev
#SBATCH --qos=debug
#PBS -l walltime=1:00:00
##PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -j eo
#BSUB -J m2m1c
#BSUB -n 384
#BSUB -W 5:00
#BSUB -o m2m1c.log.o%J
#BSUB -e m2m1c.log.o%J
# ------------------------------


/usr/local/other/PoDS/PoDS/pods.py -x convert.list -n 4

