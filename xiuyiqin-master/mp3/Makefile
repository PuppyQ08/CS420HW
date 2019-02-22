SHELL = /bin/sh

CC = cc
CFLAGS = -O3 -h std=c99 -h acc_model=auto_async_none
LDFLAGS = -O3

TARGETS = mm
OBJS = mm.o mm1.o mm2.o mm3.o support.o stats.o
OTHER = mm2_1.cub mm2_1.ptx mm3_1.cub mm3_1.ptx

.PHONY: all
all: $(TARGETS)

mm: $(OBJS)

.PHONY: clean
clean:
	rm $(OBJS) $(TARGETS) $(OTHER)
