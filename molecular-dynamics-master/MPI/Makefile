SHELL = /bin/sh 
CC = mpiicc
CFLAGS = -O3 -ipo -std=c++11 -Wall 

dmpi : mpi.cpp mdmpi.h
	$(CC) $(CFLAGS) mpi.cpp -o dmpi


clean : 
	-rm -f *.o dmpi
