##
###################################################################
##                                                               ##
##  Makefile for RIMIiCal Package                                ##
##                                                               ##
###################################################################
##
##  RISMICALDIR  RISMiCal Distribution Directory
##  BINDIR       Hard Copies of RISMiCal Executables
##  LINKDIR      Linked Copies of RISMiCal Executables
##

RISMICALDIR = ~/software/RISMiCal
BINDIR = $(RISMICALDIR)/bin

###################################################################
#
#   Machine dependent compile option
#
###################################################################


# PGI Compiler
#FC = pgf95
#FCFLAGS = -fastsse -O4 -Mprefetch=distance:8,nta -Minline=size:50 -mp -Mfixed -mcmodel=medium 
#LDFLAGS = -lacml -lacml_mv 

# GNU Compiler
FC =  gfortran
FCFLAGS = -O3 -fomit-frame-pointer -ffixed-form 
LDFLAGS = -lblas -llapack

# Intel Compiler
# FC =  ifort
# #FCFLAGS = -O3 -shared-intel -mcmodel=medium -fixed
# FCFLAGS = -O3 -shared-intel -mcmodel=medium -openmp -fixed
# LDFLAGS =  -mkl 


FCLD = $(FC)

###################################################################

SRCS = calxvk.f  output3d.f  readguess3d.f  cl_oz1d.f  outputvv.f  readhvk.f \
       closure3d.f  oz3d.f  readinput.f  drism.f  potential3duv.f  readuvdata.f\
       fbond3duv.f  potentialuv.f  readvvdata.f  fbonduv.f  potentialvv.f\
       rismical.f  fbondvv.f  printlogo.f  rismical1d.f  makewx.f  prop1duv.f\
       rismical3d.f  mathlib.f  prop1dvv.f  rismicalvv.f  mdiis1d.f  prop3duv.f \
       rismiofile.f  mdiis3d.f  setup1dvx.f  misc.f  read3duvdata.f\
       setuparraysize.f  output1d.f  readguess1d.f  uvbind3d.f

OBJS = $(SRCS:.f=.o)

PROGS =rismical.x 


all: $(PROGS)

$(PROGS) : $(OBJS)
	$(FCLD) $(FCFLAGS) $(LDFLAGS) $(OBJS) -o $(PROGS) 

clean:
	-rm -f $(OBJS)  $(PROGS)

.SUFFIXES: .f .o

.f.o:
	$(FC) $(FCFLAGS) -c $< -o $@ 

