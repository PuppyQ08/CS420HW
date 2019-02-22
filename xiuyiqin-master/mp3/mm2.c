#include <stddef.h>
#include <omp.h>
#include "mm.h"
#include "stats.h"
#include "support.h"
#include <stdio.h>
/* lower triangular matrix multiplication */
void mm2(float *restrict A, const float *restrict B,
           const float *restrict C, const size_t N)
{
	const size_t tN = triangular(N);
	start_stats();
	#pragma omp target teams distribute simd map(tofrom:A[0:tN]) map(to:B[0:tN],C[0:tN])
	//#pragma omp target teams distribute parallel for simd schedule(dynamic)
	for (size_t i = 0; i < N; i++) {
		size_t ii = triangular(i);
		for (size_t j = 0; j <=i; j++) {
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
