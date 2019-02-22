#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mm.h"
#include "stats.h"
#include "support.h"

int main(int argc, char *argv[])
{
	int rank;
	const size_t N =4000;
	double *A = NULL;
	double *B = NULL;
	double *C = NULL;
	double *M = NULL;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		A = malloc(sizeof(double) * triangular(N));
		B = malloc(sizeof(double) * triangular(N));
		C = malloc(sizeof(double) * triangular(N));
		M = malloc(sizeof(double) * triangular(N));
	}

	init(B, C, N);
	clear(A, N);

	root_printf("MM1:\n");
	mm1(A, B, C, N);
	print_stats();
	root_printf("\n");

	/* Copy A to M */
	if (rank == 0) {
		memcpy(M, A, sizeof(double) * triangular(N));
	}
	
	clear(A, N);
	root_printf("MM2:\n");
	mm2(A, B, C, N);
	print_stats();
	check(M, A, N);
	root_printf("\n");
	
	clear(A, N);
	root_printf("MM3:\n");
	mm3(A, B, C, N);
	print_stats();
	check(M, A, N);
	root_printf("\n");

	free(A);
	free(B);
	free(C);
	free(M);

	MPI_Finalize();

	return 0;
}
