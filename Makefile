CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function
DFLAGS=		
OBJS=		khmm.o kmin.o cli.o core.o em.o aux.o
PROG=		beta-psmc
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

beta-psmc:$(OBJS) main.o
		$(CC) $(CCFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)


khmm.o:khmm.h
kmin.o:kmin.h
cli.o core.o aux.o em.o:psmc.h khmm.h

