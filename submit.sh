#!/bin/bash 
#SBATCH -p long
#SBATCH -J yelmo 
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#SBATCH --mem=18000
#SBATCH -t 15-02:30:15

##SBATCH -p long_nolimits
##SBATCH -q qnolimits
##SBATCH -J ismip6_AIS 
##SBATCH -o yelmo.out
##SBATCH -e yelmo.err
##SBATCH -t 14-23:30:15

# Run the job
#./yelmox.x yelmo_Antarctica.nml
#./yelmox.x yelmo_Antarctica_sia.nml
./yelmox.x yelmo_Antarctica_restart.nml
#./yelmox.x yelmo_Antarctica_ctrl.nml

# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ytopo="kt_ref=0.0025,0.005,0.01,0.05,0.1"
