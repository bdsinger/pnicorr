# Makefile for pnicorr
# edit Make.user to specify your BLAS_{C,LD}FLAGS
include Make.user

CC ?= cc
MPICC ?= mpicc
DBGFLAGS ?= -DDEBUG -O0 -g
OPTFLAGS ?= -O3
CFLAGS ?= -Wall -Wextra $(OPTFLAGS) -std=c11 $(BLAS_CFLAGS)
LDFLAGS ?= $(BLAS_LDFLAGS) -lz
PROG ?= pnicorr
MPIPROG ?= pnicorr_mpi

MPISRC= pnicorr.c
OTHSRC = pnicorr_io.c pnicorr_debug.c opts2struct/opts2struct.c
SOURCES = $(MPISRC) $(OTHSRC)
OBJECTS = $(SOURCES:.c=.o)
MPIOBJ = $(MPIPROG).o
OTHOBJ = $(OTHSRC:.c=.o)

all: $(SOURCES) $(PROG)
clean:
	rm -f $(OBJECTS) $(MPIOBJ) $(PROG) $(MPIPROG) *~
$(PROG): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
$(MPIOBJ): $(MPISRC)
	$(MPICC) -c $(CFLAGS) -DHAVE_MPI $< -o $@
$(MPIPROG): $(MPIOBJ) $(OTHOBJ)
	$(MPICC) $(MPIOBJ) $(OTHOBJ) $(LDFLAGS) -o $@
.PHONY: clean

