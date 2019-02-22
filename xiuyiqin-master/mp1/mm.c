#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mm.h"
#include "stats.h"
#include "support.h"

int main(int argc, char *argv[])
{
	const size_t N = 3000;
	double *A = malloc(sizeof(double) * N*(N+1)/2);
	double *B = malloc(sizeof(double) * N*(N+1)/2);
	double *C = malloc(sizeof(double) * N*(N+1)/2);
	double *M = malloc(sizeof(double) * N*(N+1)/2);

	clear(A, N);
	init(B, C, N);
	init_stats();

	printf("MM1:\n");
	mm1(A, B, C, N);
	print_stats();
	printf("\n");

	/* Copy A to M */
	memcpy(M, A, sizeof(double) * N*(N+1)/2);
	
	clear(A, N);
	printf("MM2:\n");
	mm2(A, B, C, N);
	print_stats();
	check(M, A, N);
	printf("\n");
	
	clear(A, N);
	printf("MM3:\n");
	mm3(A, B, C, N);
	print_stats();
	check(M, A, N);
	printf("\n");

	clear(A, N);
	printf("MM4:\n");
	mm4(A, B, C, N);
	print_stats();
	check(M, A, N);
	printf("\n");

	free(A);
	free(B);
	free(C);
	free(M);

	return 0;
}
