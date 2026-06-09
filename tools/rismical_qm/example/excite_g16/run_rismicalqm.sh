#/bin/bash

/foo/baa/RISMiCal-QM-g16.py pna.inp &> pna.out

cp pna.qv pna_fc.qv
cp pna.chk pna_fc.chk
/foo/baa/RISMiCal-QM-g16.py pna_fc.inp &> pna_fc.out
