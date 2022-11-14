#!/bin/bash 
# Submission script for Dragon2
#SBATCH --job-name=yelmox_v1.8
##SBATCH --partition=long  # batch or long
##SBATCH --time=7-23:00:00 # hh:mm:ss
#SBATCH --partition=batch  # batch or long
#SBATCH --time=0-4:00:00 # hh:mm:ss
#
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000 # megabytes
#
##SBATCH --mail-user=javier.blasco.navarro@ulb.be
##SBATCH --mail-type=ALL

# Run the job
#./yelmox.x yelmo_Antarctica_pico.nml
#./yelmox.x yelmo_Antarctica_restart.nml 
#./yelmox.x yelmo_Antarctica_lgm.nml
#./yelmox.x yelmo_Antarctica_melt_ensemble.nml
./yelmox.x yelmo_Antarctica_deglaciation.nml

# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ydyn="cb_z0=-400 cb_z1=500 beta_u0=100"
