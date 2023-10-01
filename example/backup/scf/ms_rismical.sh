#/bin/bash

#
# Set environment variables
#
export RISMICALHOME="/Users/noriwo/software/RISMiCal/"
export PATH=$RISMICALHOME/bin/:$PATH

#
# set VQE
#
vqerun="python3 vqe_hogehoge"

#
# Set input
#
# RISM solvent
solvent="h2o_vv"
#
# RISM solute 
solute_3d="h2o_3d"
#
# VQE solute
solute_qm="test.json"
#
# solute mol
solute="h2o"
#
# maximum 3D-RISM-SCF interation
maxit=300     
#
# 3D-RISM-SCF convergence criteria
conv=0.000001


#
# RISM for solvent
#
# Input : ${solvent}.inp
# Output: ${solvent}.xvk
rismical.x vv ${solvent}.inp

#
# VQE in gas phase
#
# Input : ${solute_qm}
# Output: ${solute_3d}.esp
$vqerun ${solute_qm}

#
# 3D-RISM-SCF iteration
#
for itr in `seq 1 "$maxit"` ; do 

#
# 3D-RISM
#
# Input : ${solute_3d}.inp ${solute_3d}.esp ${solute}.xyz ${solute}.lj ${solute_3d}.xvk
# Output: ${solute_3d}.qv ${solute_3d}.guv ${solute_3d}.log
rismical.x 3d ${solute_3d}.inp &> ${solute_3d}.log

#
# VQE
#
# Input : ${solute_qm} ${solute}.qv 
# Output: ${solute_3d}.esp 
$vqerun ${solute_qm}

#
# Check convergence
#
#if [ hoge_hoge ]
#then
#  echo "3D-RISM-SCF converged"
#  break
#fi

done
