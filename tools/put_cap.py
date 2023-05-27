#!/usr/local/bin/python3

import numpy as np

# INPUT
# Coordinate 3 atoms
A=np.array([3.326,1.548,-0.000])  ## N  (C )
B=np.array([3.909,0.724,-0.000])  ## H  (O )
C=np.array([3.970,2.846,-0.000])  ## CA (CA)
r=1.3

# calc length of three sides
CB=C-B
CA=C-A
BA=B-A
a=np.linalg.norm(CB)
b=np.linalg.norm(CA)
c=np.linalg.norm(BA)

# Calc in-center 

E=1.0/(a+b+c)*(a*A+b*B+c*C)

# Vector AE and its norm
AE=A-E
ae=np.linalg.norm(AE)

# Coordinate of D

D=A+r/ae*AE

print(f"ATOM     XX  C   ACE     X    {D[0]:8.3f}{D[1]:8.3f}{D[2]:8.3f}  1.00  0.00")








