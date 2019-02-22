SHELL = /bin/sh

CC = icc
# Change the following line to change the optimization flags
CFLAGS = -std=gnu99 -O3 -qopt-report=5 -S
LDLIBS = -lrt

TARGETS = mp0_b

all: $(TARGETS)

.PHONY: clean
clean:
	rm $(TARGETS)
