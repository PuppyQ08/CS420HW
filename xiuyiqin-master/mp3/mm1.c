#include <stddef.h>

#include "mm.h"
#include "stats.h"
#include "support.h"

/* lower triangular matrix multiplication */
void mm1(float *restrict A, const float *restrict B,
           const float *restrict C, const size_t N)
{
	#pragma omp parallel
	{
		start_stats();

		#pragma omp for schedule(dynamic)
		for (size_t j = 0; j < N; j++) {
			for (size_t i = j; i < N; i++) {
				size_t ii = triangular(i);
				float sum = 0.0;
				for (size_t k = j; k <= i; k++) {
					size_t kk = triangular(k);
					sum += B[ii + k] * C[kk + j];
				}
				A[ii + j] = sum;
			}
		}

		stop_stats();
	}
}
