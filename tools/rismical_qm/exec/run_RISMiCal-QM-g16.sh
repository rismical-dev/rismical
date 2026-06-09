#!/bin/bash
### 
#PBS -q cpu
#PBS -l select=1:ncpus=20
#PBS -j oe

export GAUSS_PDEF=$NCPUS

cd $PBS_O_WORKDIR

#
# usage: qsub -v arg1=inputfile.inp run_RISMiCal-QM.sh 
#

output=`basename $arg1 .inp`.out

python3 ~/software/RISMiCal-local/RISMiCal-QM/exec/RISMiCal-QM-g16.py $arg1 $arg2 &> $output
