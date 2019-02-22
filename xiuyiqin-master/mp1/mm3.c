#include <stddef.h>

#include "mm.h"
#include "stats.h"

/* lower triangular matrix multiplication */
void mm3(double *restrict A, const double *restrict B,
           const double *restrict C, size_t N)
{
	{
		start_stats();

		for (size_t j = 0; j < N; j++) {
			for (size_t i = j; i < N; i++) {
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
