#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define ITER 10

double get_clock(void);

void test1(void);
void test2(void);
void test3(void);
void test4(void);
void test5(void);

int main(int argc, char *argv[])
{
	test1();
	test2();
	test3();
	test4();
	test5();
	return 0;
}
// to run a*a 10E6 times repeat r times
void test1(void)
{
	double total_time = 0.0;

	for (size_t r = 0; r < ITER; r++) {
		int a = 2;
		int a2 = 1;

		double t = get_clock();
		for (size_t i = 0; i < 1000000; i++) {
			a2 = a*a;
		}
		t = get_clock() - t;

		total_time += t;
	}

	printf("Time taken for test1: %7.16fs\n", total_time/ITER);
}
//to run array adding 10E6 repeat r times
void test2(void)
{
	const size_t N = 1000000;
	int *A = malloc(sizeof(int) * N);
	int *B = malloc(sizeof(int) * N);
	int *C = malloc(sizeof(int) * N);
	double total_time = 0.0;

	srand(time(NULL));
	for(size_t i = 0; i < N; i++) {
		A[i] = rand() % N;
		B[i] = rand() % N;
	}

	for (size_t r = 0; r < ITER; r++) {

		double t = get_clock();
		for (size_t i = 0; i < N; i++) {
			C[i] = A[i] + B[i];
		}
		t = get_clock() - t;

		total_time += t;
	}

	free(A);
	free(B);
	free(C);

	printf("Time taken for test2: %7.16fs\n", total_time/ITER);
}
//to init array ~10E6 but in two loops. loop unrolling?
void test3(void)
{
	const size_t N = 1024;
	int *A = malloc(sizeof(int) * N * N);
	double total_time = 0.0;

	for (size_t r = 0; r < ITER; r++) {

		double t = get_clock();
		for (size_t j = 0; j < N; j++) {
			for (size_t i = 0; i < N; i++)
				A[i*N + j] = 1000;
		}
		t = get_clock() - t;

		total_time += t;
	}

	free(A);

	printf("time taken for test3: %7.16fs\n", total_time/ITER);
}
//to copy ~10E6 array also in two loops.
void test4(void)
{
	const size_t N = 1024;
	int *A = malloc(sizeof(int) * N * N);
	int *B = malloc(sizeof(int) * N * N);
	double total_time = 0.0;

	srand(time(NULL));
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			A[i*N + j] = rand() % 1000;

	for (size_t r = 0; r < ITER; r++) {

		double t = get_clock();
		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < N; j++)
				B[j*N + i] = A[i*N + j];
		}
		t = get_clock() - t;

		total_time += t;
	}

	free(A);
	free(B);

	printf("Time taken for test4: %7.16fs\n", total_time/ITER);
}
// two loops of 10000 array sequentilally
void test5(void)
{
	const size_t N = 10000;
	int *A = malloc(sizeof(int) * N);
	double total_time = 0.0;

	int t, i;
	for (size_t r = 0; r < ITER; r++) {

		double t = get_clock();
		for (size_t i = 0; i < N; i++) {
			A[i] = 25;
		}
		for (size_t i = 0; i < N; i++) {
			A[i] *= 2;
		}
		t = get_clock() - t;

		total_time += t;
	}

	free(A);

	printf("Time taken for test5: %7.16fs\n", total_time/ITER);
}

double get_clock(void)
{
	struct timespec tp;
	int ok = clock_gettime(CLOCK_MONOTONIC, &tp);
	if (ok != 0) {
		err(ok, "clock_gettime");
	}

	return (tp.tv_sec*1.0 + tp.tv_nsec*1.0E-9);
}
