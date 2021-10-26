#!/bin/bash 
#SBATCH -p long
#SBATCH -J yelmo 
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#SBATCH --mem=1800
#SBATCH -t 7-02:30:15

# Run the job
#./yelmox.x yelmo_Antarctica_pico.nml
#./yelmox.x yelmo_Antarctica_restart.nml 
#./yelmox.x yelmo_Antarctica_lgm.nml
#./yelmox.x yelmo_Antarctica_melt_ensemble.nml
./yelmox.x yelmo_Antarctica_deglaciation.nml

# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ydyn="cb_z0=-400 cb_z1=500 beta_u0=100"
