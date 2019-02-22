#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mm.h"
#include "stats.h"
#include "support.h"

int main(int argc, char *argv[])
{
	//const size_t N = 1000;
for(int N = 1000; N <=16000; N*=2){
	
	float *A = malloc(sizeof(float) * N*(N+1)/2);
	float *B = malloc(sizeof(float) * N*(N+1)/2);
	float *C = malloc(sizeof(float) * N*(N+1)/2);
	float *M = malloc(sizeof(float) * N*(N+1)/2);
	
	clear(A, N);
	init(B, C, N);
	init_stats();

	printf("total number is %d\n",N);
	printf("MM1:\n" );
	mm1(A, B, C, N);
	print_stats();
	printf("\n");
	
	/* Copy A to M */
	memcpy(M, A, sizeof(float) * N*(N+1)/2);
	/*	
	clear(A, N);
	printf("MM2:\n");
	mm2(A, B, C, N);
	print_stats();
	check(M, A, N);
	printf("\n");
	*/
	/*
	for(int i = 0; i < N*(N+1)/2; i++){
	printf("%f \n",M[i]-A[i]);
	}
	*/
	

	clear(A, N);
	printf("MM3:\n");
	mm3(A, B, C, N);
	print_stats();
	check(M, A, N);
	printf("\n");


	free(A);
	free(B);
	free(C);
	free(M);
}	
	return 0;
}
