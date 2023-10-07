CC      = g++

CFLAGS  = -O3 -Xpreprocessor  -fopenmp -lomp -fomit-frame-pointer -std=c++11
#CFLAGS  = -O3 -Xpreprocessor -fopenmp -lomp -fomit-frame-pointer -std=c++11

OBJS	= main.o set_fname.o

PROGRAM = guv2dx
.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

all:    $(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(INCLUDES) $(LIBS) -o $(PROGRAM)

clean:; @rm -f $(OBJS) core

$(OBJS) :
