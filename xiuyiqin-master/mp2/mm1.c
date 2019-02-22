#include <mpi.h>
#include <stddef.h>

#include "mm.h"
#include "stats.h"
#include "support.h"

/* lower triangular matrix multiplication */
void mm1(double *restrict A, const double *restrict B,
           const double *restrict C, size_t N)
{
	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	start_stats();

	if (rank == 0) {
		for (size_t j = 0; j < N; j++) {
			for (size_t i = j; i < N; i++) {
				size_t ii = triangular(i);
				double sum = 0.0;
				for (size_t k = j; k <= i; k++) {
					size_t kk = triangular(k);
					sum += B[ii + k] * C[kk + j];
				}
				A[ii + j] = sum;
			}
		}
	}

	stop_stats();
}
