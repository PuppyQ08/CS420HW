SHELL = /bin/sh 
CC = mpiicc
CFLAGS = -O3 -ipo -std=c++11 -Wall 

mdmpi : main_mpi.cpp mdmpi.h
	$(CC) $(CFLAGS) main_mpi.cpp -o mdmpi

md : main_seq.cpp sequential.h
	$(CC) $(CFLAGS) main_seq.cpp -o md

clean : 
	-rm -f *.o mdmpi md
