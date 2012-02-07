# Makefile

#############################################################
## modify to match your compiler defaults 

F90 = gfortran 
CC = gcc

# modify to point to your own SAC installation path
SAC_INC = /usr/local/sac/sac/include
SAC_LIB = /usr/local/sac/sac/lib

#############################################################

# function library
F90FLAGS = -O2
CCFLAGS = -O2

CSRC = rw_ascfile_c
F90SRC = rw_ascfile_f

F90OBJ_LIB = $(patsubst %,%.o,$(F90SRC))
COBJ_LIB = $(patsubst %,%.o,$(CSRC))
OBJ_LIB = $(F90OBJ_LIB) $(COBJ_LIB)

#
# program asc2sac
#
PROG = asc2sac

SRCS = $(PROG).c
OBJS = ${SRCS:%.c=%.o}

CVS_INCDIR = $(SAC_INC) #formats
CVS_LIBDIR = $(SAC_LIB) #formats

INCS = -I$(CVS_INCDIR) -DSAC_NULL
CFLAGS = $(INCS)
LIBS = -L$(CVS_LIBDIR) -lsac -lsacio -lm #-lformats 


all :  $(PROG)

lib : $(OBJ_LIB)

#
# convolve_stf
#
convolve_stf: convolve_stf.c
	$(CC) $(CFLAGS) -o convolve_stf convolve_stf.c $(LIBS)

#
# program asc2sac
#
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS) $(LIBS)

#
# library routines
#
$(F90OBJ_LIB) : %.o : %.f90
	$(F90)  $(F90FLAGS) -c $*.f90

$(COBJ_LIB): %.o : %.c
	$(CC) $(CCFLAGS) -c $*.c 


clean: 
	\rm -f $(OBJS) *~ *.o $(PROG) $(OBJS) *~ core convolve_stf

cleanall:
	\rm -f $(OBJS)  *.o *~ $(PROG) $(OBJS) *~ core convolve_stf
