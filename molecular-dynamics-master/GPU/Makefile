CC = CC
CFLAGS = -O3 -h std=c++11 -h acc_model=auto_async_none

gpu: gpumd.cpp functions.h 
	$(CC) $(CFLAGS) gpumd.cpp -o gpu

        
clean:
	-rm -f *.o gpu 
