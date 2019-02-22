SHELL = /bin/sh

CC = icc
CFLAGS = -O3 -ipo -std=gnu99 -fopenmp -Wall
LDFLAGS = -O3 -ipo -fopenmp

TARGETS = mm mm-papi
OBJS = mm.o mm1.o mm2.o mm3.o mm4.o support.o
mm_OBJS = stats.o
mm_papi_OBJS = stats_papi.o

.PHONY: all
all: $(TARGETS)

mm: LDLIBS += -lrt
mm: $(OBJS) $(mm_OBJS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

mm-papi: LDFLAGS += -L/usr/local/papi/5.6.0/lib
mm-papi: CFLAGS += -I/usr/local/papi/5.6.0/include
mm-papi: LDLIBS += -lpapi
mm-papi: $(OBJS) $(mm_papi_OBJS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

.PHONY: clean
clean:
	rm $(OBJS) $(mm_OBJS) $(mm_papi_OBJS) $(TARGETS)
