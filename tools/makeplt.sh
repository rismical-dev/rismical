#!/bin/bash
#
# $1 number of 3dgrid 
# $2 size of 3d grid width
# $3 number of symmetry uniq solvent atom
# $4 type of UVDATA format
# $5 Name of UVDATA (Usually "UVDATA")
# 

gOpenMolbin=/home/noriwo/software/gOpenMol-3.00/bin/pltfile
makepltbin=/home/noriwo/software/gamess-2016/tools/rism/makeplt.x



for i in `seq 1 $3`
do

$makepltbin $1 $2 $3 $i $4  < $5 > temp_mkplt
$gOpenMolbin -fu -itemp_mkplt -o$5.$i.plt

done

rm temp_mkplt
