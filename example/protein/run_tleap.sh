#!/bin/bash

prestep=step1
step=step2

cp ../${prestep}_modify-pdb/*_mod.pdb ./
tleap -f leap.in
