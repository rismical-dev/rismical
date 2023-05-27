#! /usr/bin/python3

import sys

args=sys.argv

# Check arguments

if len(args) <= 1:
    print ("Insufficient arguments.")
    print ("amb2rismcal.py [parm7] [pdb|crd] ([format])")
    exit()

# 
parm7file=args[1]
crdfile=args[2]
#
if len(args)>= 3:
    xformat=args[3]
else:
    xformat=0

