# RISMiCal
The Reference Interaction-Site Model integrated Calculator

## INSTALL

Go to install directory, and
```
git clone https://github.com/rismical-dev/rismical.git
```
Set enviromental valuable RISMICALHOME to your installation directory
```
export RISMICALHOME="/foo/bar/RISMiCal"
```
Go to source directory
```
cd $RISMICALHOME/source
```
Make your binary
```
make
```
`Makefile` should be modified to fit your system before `make`.

## QUICK START
Go to exmaple directory
```
cd $RISMICALHOME/example
```
Usually, any RISM/3D-RISM calculation begins with a calculation of the solvent (vv) system.
Let's do the solvent calculation for water.
```
cd h2o_vv
$RISMICALHOME/source/rismical.x vv h2o_vv.inp &> h2o_vv.log
```
1st argument of `rismical.x` specifies the system to be computed.
You can obtain the solvent susptibility function `h2o_vv.xsv`.
For the system having a solute molecule immersed in solvent at infinite dilution, solute-solvent RISM or 3D-RISM is performed. 
To run the solute-solvent system, `.xsv` file is needed. 
For example, 
```
cd $RISMICALHOME/example/h2o_3d
cp $RISMICALHOME/example/h2o_vv/h2o_vv.xsv ./
$RISMICALHOME/source/rismical.x 3d h2o_3d.inp &> h2o_3d.log
```
You can find some examples in `$RISMICALHOME/example/`.
The details of the input format is explained in the manual in `$RISMICALHOME/doc/`.

# Copyright
Norio YOSHIDA, Nagoya University

# Contact Information
Norio YOSHIDA, Prof. of Graduate School of Informatics, Nagoya University
noriwo@nagoya-u.jp
https://sites.google.com/view/yoshida-group/homenoriwo@nagoya-u.jp
