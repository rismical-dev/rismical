##
###################################################################
##                                                               ##
##  Makefile for RIMIiCal Tools                                  ##
##                                                               ##
###################################################################
##
##  RISMICALDIR  RISMiCal Distribution Directory
##  BINDIR       Hard Copies of RISMiCal Executables
##  LINKDIR      Linked Copies of RISMiCal Executables
##

RISMICALDIR = ~/software/RISMiCal
BINDIR = $(RISMICALDIR)/bin
OBJDIR = $(RISMICALDIR)/object

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
FCFLAGS = -O3 -fomit-frame-pointer -ffixed-form -Wno-argument-mismatch -std=legacy 
# FCFLAGS = -O3 -fomit-frame-pointer -ffixed-form -fallow-argument-mismatch -std=legacy 
LDFLAGS = -lblas -llapack

# Intel Compiler
#FC =  ifort
#FCFLAGS = -O3 -shared-intel -mcmodel=medium -fixed
#FCFLAGS = -O3 -shared-intel -qopenmp -fixed
# # LDFLAGS =  -mkl 


FCLD = $(FC)

###################################################################

SRCS = sdf2rdf.f

OBJS = $(SRCS:.f=.o)
#OBJS = $(OBJDIR)/$(SRCS:.f=.o)

PROGS =sdf2rdf.x 


all: $(PROGS)

$(PROGS) : $(OBJS)
	$(FCLD) $(FCFLAGS) $(LDFLAGS) $(OBJS) -o $(PROGS) 

clean:
	-rm -f $(OBJS)  $(PROGS)

.SUFFIXES: .f .o

.f.o:
	$(FC) $(FCFLAGS) -c $< -o $@ 

copy:
	cp $(PROGS)  $(BINDIR)/$(PROGS)

