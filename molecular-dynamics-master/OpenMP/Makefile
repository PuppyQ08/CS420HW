CC = icc
CFLAGS = -O3 -ipo -std=c++11 -fopenmp -Wall
LDFLAGS = -O3 -ipo -fopenmp

mdop : main_OpenMP.cpp functions.h 
	$(CC) $(CFLAGS) main_OpenMP.cpp -o mdop

md : main_seq.cpp functions.h
	$(CC) $(CFLAGS) main_seq.cpp -o md
        

        
clean:
	-rm -f *.o mdop md
