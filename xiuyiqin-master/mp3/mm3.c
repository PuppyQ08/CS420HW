#include <omp.h>
#include <stddef.h>
#include <stdio.h>

#include "mm.h"
#include "stats.h"
#include "support.h"

/* lower triangular matrix multiplication */
void mm3(float *restrict A, const float *restrict B,
           const float *restrict C, const size_t N)
{
	const size_t tN = triangular(N);
	start_stats();
	double send_data = omp_get_wtime();
	#pragma omp target enter data map(to:A[0:tN],B[0:tN],C[0:tN])
//	#pragma omp target teams distribute simd map(tofrom:A[0:tN]) map(to:B[0:tN],C[0:tN]) 
	send_data = omp_get_wtime() - send_data;
	printf("Send Data time: %fs \n", send_data);
	//Next is on GPU
	#pragma omp target teams distribute simd 
	for (size_t i = 0; i < N; i++) {
		size_t ii = triangular(i);
		for (size_t j = 0; j <= i; j++) {
			float sum = 0.0;
			for (size_t k = j; k <= i; k++) {
				size_t kk = triangular(k);
				sum += B[ii + k] * C[kk + j];
			}
			A[ii + j] = sum;
		}
	}
	//next is back to CPU
	//#pragma omp target map(delete:B[0:tN],C[0:tN])
	double get_data = omp_get_wtime();
	#pragma omp target exit data map(from: A[0:tN]) map(delete:B[0:tN],C[0:tN]) 	
	get_data = omp_get_wtime() - get_data;
	printf("Get Data time: %fs \n", get_data);
	//#pragma omp target map(delete:B[0:tN],C[0:tN])
	stop_stats();
}
