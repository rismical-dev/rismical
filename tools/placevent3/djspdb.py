'''
    Dan Sindhikara sindhikara@gmail.com
    Does PDB related stuff
    Todo: 
    read pdb
    atom class?
'''

def makepdbstring(atomserial,atomname,resname,chainid,resseqnum,x,y,z,occ,tfact,element,charge):
    '''
    Doesn't use altlocation,insertion code
    '''
    pdbstring="HETATM%5d %4s %3s%2c %3d    %8.3f%8.3f%8.3f%6.2f %05.2f           %1s%2d\n" % (atomserial,atomname,resname,chainid,resseqnum,x,y,z,occ,tfact,element,charge)
    return(pdbstring)

