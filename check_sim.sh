#!/bin/bash

fldr=$1
file_out=${fldr}/ens_err.txt 

echo "Processing results in: ${fldr}"

#for D in ${fldr}/* ; do ./check_sim.x $D ; done > ${file_out}

# Reset file to be empty 
cat <<EOF > ${file_out} 

# Loop over subfolders in fldr 
for D in ${fldr}/*
do
    echo $D 
    ./check_sim.x $D >> ${file_out}
done 

echo "Results written to: ${file_out}"