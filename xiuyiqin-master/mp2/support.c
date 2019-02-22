#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "support.h"

#define EPSILON (1.0e-6)

int root_printf(const char *restrict format, ...)
{
	int rank;
	int retval = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		va_list args;
		va_start(args, format);
		retval = vprintf(format, args);
		va_end(args);
	}

	return retval;
}

void init(double *restrict B, double *restrict C, size_t N)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		srand(time(NULL));
//		srand(0xBAD5EED);
		/* initialize B and C to random values */
		for (size_t i = 0; i < triangular(N); i++) {
			B[i] = (double)rand() / (double)rand();
			C[i] = (double)rand() / (double)rand();
		}
	}
}

void check(const double *restrict A, const double *restrict B, size_t N)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		for (size_t i = 0; i < N; i++) {
			size_t ii = triangular(i);
			for (size_t j = 0; j <= i; j++) {
				double a_val = A[ii + j];
				double b_val = B[ii + j];
				if (fabs(1.0 - b_val/a_val) > EPSILON) {
					printf("TEST FAILED "
					       "at (%zu, %zu): %f != %f\n",
					       i, j, a_val, b_val);
					return;
				}
			}
		}
		printf("TEST PASSED\n");
	}
}

void clear(double *restrict A, size_t N)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		for (size_t i = 0; i < triangular(N); i++) {
			A[i] = 0.0;
		}
	}
}
