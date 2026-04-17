'''Daniel J. Sindhikara
    sindhikara@gmail.com
    Does stuff related to dx files
    '''

def readdx(filenames):
    '''
    Reads one or more dx files into memory
    returns grid data and single 4D array containing dx data [index][xindex][yindex][zindex]
    '''
    import numpy as np
    dxfile = open(filenames[0],"r")
    dxlines = []
    for i in range(10): # only need the first few lines to get grid data
        dxlines.append(dxfile.readline())
    dxfile.close()
    gridcount = []
    origin = []
    delta = [0,0,0]
    startline=0
    for i,line in enumerate(dxlines) :
        splitline = line.split()
        if len(splitline) > 2 :
            if(splitline[1]=="1") :
                gridcount.append(int(splitline[5]))
                gridcount.append(int(splitline[6]))
                gridcount.append(int(splitline[7]))
                print(("# gridcounts ",gridcount))
            if(splitline[0]=="origin") :
                origin.append(float(splitline[1]))
                origin.append(float(splitline[2]))
                origin.append(float(splitline[3]))
            if(splitline[0]=="delta") :
                if(float(splitline[1])>0) :
                    delta[0]=float(splitline[1])
                if(float(splitline[2])>0) :
                    delta[1]=float(splitline[2])
                if(float(splitline[3])>0) :
                    delta[2]=float(splitline[3])
            if(splitline[1]=="3") :
                numpoints=int(splitline[9])
                print(("# Total number of gridpoints is ",numpoints))
                startline=i+1
        if(startline>1) :
            break
    distributions = [[[[0.0 for x in range(gridcount[2])] for y in range(gridcount[1])] for z in range(gridcount[0])]for w in range(len(filenames))]
    gridvolume=delta[0]*delta[1]*delta[2]
    print(("# I have to read",len(filenames)*gridcount[2]*gridcount[1]*gridcount[0],"values."))
    for i, dxfilename in enumerate(filenames): 
        dxfile = open(dxfilename,"r")
        dxtext = dxfile.read()
        dxfile.close()
        splittext=dxtext.split()
        del splittext[0:splittext.index("follows")+1] # get rid of header text, last word is "follows"
        del splittext[-4:-1]
        del splittext[-1] # dont know why i cant just use a 0 instead of -1 above
        floats = [float(element) for element in splittext]
        index=0
        for x in range(gridcount[0]):
            for y in range(gridcount[1]):
                for z in range(gridcount[2]):
                    distributions[i][x][y][z]=floats[index]
                    index+=1 # there's probably a more pythonic way to do this
    print(("# corr[0][0][0][0] = ",distributions[0][0][0][0],"[0][-1][-1][-1]=",distributions[0][-1][-1][-1]))
    distributions=np.array(distributions)
    return(distributions,origin,delta,gridcount)    


def printdx(values,indices,origin,delta,gridcounts,filename):
    ''' Print a dx file'''
    '''first construct full grid as 1d list (z fast, y med, z slow)
    note values do not span entire grid, so just replce them if they exist
    values = 1D array containing all values
    indices[value index][dimension] = index in the dimension specified (as opposed to 1D index)
    origin = length 3 float list of origin coordinate
    delta = length 3 float list of grid spacing
    gridcounts = length 3 integer list of number of grid points in each dimension
    filename = string of output file name
    '''
    gridvalues=[0.0]*gridcounts[0]*gridcounts[1]*gridcounts[2]
    for i,value in enumerate(values):
        gridindex=indices[i][2]+indices[i][1]*gridcounts[2]+indices[i][0]*gridcounts[2]*gridcounts[1]
        gridvalues[gridindex]=value

    f = open(filename,"w")
    f.write("#DX file from Dan's program\n")
    f.write("object 1 class gridpositions counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("origin {0} {1} {2}\n".format(origin[0],origin[1],origin[2]))
    f.write("delta {0} 0 0\n".format(delta[0]))
    f.write("delta 0 {0} 0\n".format(delta[1]))
    f.write("delta 0 0 {0}\n".format(delta[2]))
    f.write("object 2 class gridconnections counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("object 3 class array type double rank 0 items {0} data follows\n".format(gridcounts[0]*gridcounts[1]*gridcounts[2]))
    for gridvalue in gridvalues:
        f.write("{0}\n".format(gridvalue))
    f.write("object {0} class field\n".format(filename))
    f.close()

