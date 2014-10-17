# Makefile for pnicorr
# edit Make.user to specify your BLAS_{C,LD}FLAGS
include Make.user

CC ?= cc
CFLAGS ?= -Wall -Wextra -O3 -std=c99 $(BLAS_CFLAGS)
LDFLAGS ?= $(BLAS_LDFLAGS) -lz
PROG ?= pni_corr

SOURCES = pnicorr.c pnicorr_io.c pnicorr_debug.c opts2struct/opts2struct.c
OBJECTS = $(SOURCES:.c=.o)

all: $(SOURCES) $(PROG)
clean:
	rm -f *.o $(PROG) *~
$(PROG): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@
.%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
.PHONY: clean

