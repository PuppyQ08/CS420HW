#include <stddef.h>

#include "mm.h"
#include "stats.h"

/* lower triangular matrix multiplication */
void mm1(double *restrict A, const double *restrict B,
           const double *restrict C, size_t N)
{
	{
		start_stats();

		for (size_t i = 0; i < N; i++) {
			size_t ii = i * (i+1) / 2;
			for (size_t j = 0; j <= i; j++) {
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
