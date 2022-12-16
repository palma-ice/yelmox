#!/bin/bash 
# Submission script for Brigit

## No limit
##SBATCH -p long_nolimits
##SBATCH -q qnolimits
##SBATCH -J pyrenees
##SBATCH -n 12
##SBATCH -o yelmo.out
##SBATCH -e yelmo.err
##SBATCH --mem=2Gb
##SBATCH -t 14-23:30:15

## Max 30
#SBATCH -p long
#SBATCH -J deglaciation
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#SBATCH --mem=2Gb
#SBATCH -t 6-23:59:59

# Run the job
#./yelmox.x yelmo_pyrenees_1cycle.nml 
#./yelmox.x yelmo_Antarctica_lgm.nml
./yelmox.x yelmo_Antarctica_deglaciation.nml
#./yelmox_ismip6.x yelmo_ismip6_Antarctica_spinup.nml
#./yelmox_ismip6.x yelmo_ismip6_Antarctica.nml

# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ydyn="cb_z0=-400 cb_z1=500 beta_u0=100"
