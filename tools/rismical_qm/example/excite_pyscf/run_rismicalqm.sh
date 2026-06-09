#/bin/bash

/foo/baa/RISMiCal-QM-PySCF.py pna.inp &> pna.out
cp pna.qv pna_fc.qv
/foo/baa/RISMiCal-QM-PySCF.py pna_fc.inp &> pna_fc.out
