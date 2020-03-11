#!/bin/bash

fldr=$1
file_out=${fldr}/ens_err.txt 

echo "Processing results in: ${fldr}"

for D in ${fldr}/* ; do ./check_sim.x $D ; done > ${file_out}

echo "Results written to: ${file_out}"