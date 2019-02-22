#include <stddef.h>
#include <omp.h>
#include <stdio.h>
#include "mm.h"
#include "stats.h"

/* lower triangular matrix multiplication */
void mm4(double *restrict A, const double *restrict B,
           const double *restrict C, size_t N)
{
   // omp_set_num_threads(2);
    #pragma omp parallel
	{
               // printf("The number of threads: %d \n", omp_get_num_threads());
                printf("The number of thread: %d \n", omp_get_thread_num());
		start_stats();
                #pragma omp for schedule(dynamic)
		for (size_t j = 0; j < N; j++) {
			for (size_t i = j; i <N; i++) {
			        size_t ii = i * (i+1) / 2;
				double sum = 0.0;
				for (size_t k = j; k <= i; k++) {
					size_t kk = k * (k+1) / 2;
					sum += B[ii + k] * C[kk + j];
				}
				A[ii + j] = sum;
			}
		}

		stop_stats();
	}
}
