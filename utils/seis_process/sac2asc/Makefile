#########################################################
# Makefile
#
# needs SAC libraries installed 
########################################################

##
## modify to include your SAC library path
##
INCLUDE_SACPATH=-I/opt/seismo-util/lib/Sacio/single
LIBRAY_SACPATH=-L/opt/seismo-util/lib

CFLAGS = -g -Wall $(INCLUDE_SACPATH)

PROG=sac2asc

OBJS=$(PROG).o

LIBS=$(LIBRARY_SACPATH) -lSacio -lSacTools


$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) core *~ 
