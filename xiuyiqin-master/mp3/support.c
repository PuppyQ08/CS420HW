#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "support.h"

#define EPSILON (1.0e-6)

void init(float *restrict B, float *restrict C, size_t N)
{
	srand(time(NULL));
//	srand(0xBAD5EED);
	/* initialize B and C to random values */
	for (size_t i = 0; i < triangular(N); i++) {
		B[i] = ((float)rand()) / ((float)rand());
		C[i] = ((float)rand()) / ((float)rand());
	}
}

void check(const float *restrict A, const float *restrict B, size_t N)
{
	for (size_t i = 0; i < N; i++) {
		size_t ii = triangular(i);
		for (size_t j = 0; j <= i; j++) {
			float a_val = A[ii + j];
			float b_val = B[ii + j];
			if (fabs(1.0 - b_val/a_val) > EPSILON) {
				printf("TEST FAILED at (%zu, %zu): %f != %f\n",
				       i, j, a_val, b_val);
				return;
			}
		}
	}
	printf("TEST PASSED\n");
}

void clear(float *restrict A, size_t N)
{
	for (size_t i = 0; i < triangular(N); i++) {
		A[i] = 0.0;
	}
}
