#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "support.h"

#define EPSILON (1.0e-6)

void init(double *restrict B, double *restrict C, size_t N)
{
	srand(time(NULL));
//	srand(0xBAD5EED);
	/* initialize B and C to random values */
	for (size_t i = 0; i < N*(N+1)/2; i++) {
		B[i] = ((double)rand()) / ((double)rand());
		C[i] = ((double)rand()) / ((double)rand());
	}
}

void check(const double *restrict A, const double *restrict B, size_t N)
{
	for (size_t i = 0; i < N; i++) {
		size_t ii = i * (i+1) / 2;
		for (size_t j = 0; j <= i; j++) {
			double a_val = A[ii + j];
			double b_val = B[ii + j];
			if (fabs(1.0 - b_val/a_val) > EPSILON) {
				printf("TEST FAILED at (%zu, %zu): %f != %f\n",
				       i, j, a_val, b_val);
				return;
			}
		}
	}
	printf("TEST PASSED\n");
}

void clear(double *restrict A, size_t N)
{
	for (size_t i = 0; i < N*(N+1)/2; i++) {
		A[i] = 0.0;
	}
}
