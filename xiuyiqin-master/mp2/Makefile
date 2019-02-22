SHELL = /bin/sh

CC = mpiicc
CFLAGS = -O3 -ipo -std=gnu99 -Wall -fopenmp
LDFLAGS = -O3 -ipo -fopenmp

TARGETS = mm
OBJS = mm.o mm1.o mm2.o mm3.o support.o stats.o

.PHONY: all
all: $(TARGETS)

mm: $(OBJS)

.PHONY: clean
clean:
	rm $(OBJS) $(TARGETS)
