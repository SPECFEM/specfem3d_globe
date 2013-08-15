#########################################################
# Makefile
########################################################

INCLUDE_SACPATH=-I.

CFLAGS = -g -Wall $(INCLUDE_SACPATH)

PROG=sac2asc

OBJS=$(PROG).o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) core *~ 
