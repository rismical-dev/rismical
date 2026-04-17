#!/bin/bash
### sample script for pbs
#PBS -q gpu
#PBS -l select=1:ncpus=20:ngpus=1

export RISMICALHOME=~/software/RISMiCal-dev/rismical/
export PATH=$RISMICALHOME/bin:$PATH

cd $PBS_O_WORKDIR

~/software/RISMiCal-dev/rismical/bin/rismical.x 3d 4ow1.inp &> 4ow1.out
