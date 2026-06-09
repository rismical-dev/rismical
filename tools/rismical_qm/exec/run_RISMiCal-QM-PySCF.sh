#!/bin/bash
### 
#PBS -q cpu
#PBS -l select=1:ncpus=20
#PBS -j oe

export OMP_NUM_THREADS=$NCPUS

cd $PBS_O_WORKDIR

#
# usage: qsub -v arg1=inputfile.inp run_RISMiCal-QM-PySCF.sh 
# usage: qsub -v arg1=inputfile.inp,arg2=-FC run_RISMiCal-QM-PySCF.sh 
#

output=`basename $arg1 .inp`.out

source /home/noriwo/pyscf-env/bin/activate
python3 ~/software/RISMiCal-local/rismical/tools/rismical_qm/exec/RISMiCal-QM-PySCF.py $arg1 $arg2 &> $output
