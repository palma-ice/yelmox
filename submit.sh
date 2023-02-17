#!/bin/bash 
# Submission script for Dragon2
#SBATCH --job-name=pyrenees
#SBATCH --partition=long  # batch or long
#SBATCH --time=15-23:00:00 # hh:mm:ss
##SBATCH --partition=batch  # batch or long
##SBATCH --time=4-23:00:00 # hh:mm:ss
#
#SBATCH -o yelmo.out
#SBATCH -e yelmo.err
#
##SBATCH --mail-user=javier.blasco.navarro@ulb.be
##SBATCH --mail-type=ALL
#SBATCH --mem=18000

# Run the job
./yelmox.x yelmo_pyrenees_lgm.nml 
#./yelmox.x yelmo_Antarctica_lgm.nml
#./yelmox.x yelmo_Antarctica_deglaciation.nml
#./yelmox_ismip6.x yelmo_ismip6_Antarctica_spinup.nml
#./yelmox_ismip6.x yelmo_ismip6_Antarctica.nml

# Run ensemble
# python2.7 yelmo_ensemble.py -l -f -a output/deglaciation_32km/ \&ydyn="cb_z0=-400 cb_z1=500 beta_u0=100"
