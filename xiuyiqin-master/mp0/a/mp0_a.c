#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXSIZE (1UL << 26)
#define REPS (1UL << 27)
#define STRIDE 1000

void timerTest(void);
double get_clock(void);

void initValues(int *A, size_t *indirection, size_t size)
{
	srand(time(NULL));
	indirection[0] = STRIDE % size;
	for (size_t i = 1; i < size; i++) {
		A[i] = rand();
		indirection[i] = (indirection[i-1]+STRIDE) % size;
	}
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


int main(int argc, char *argv[])
{
	double t;
	int *A;
	size_t *indirection;

	A = malloc(sizeof(int) * MAXSIZE);
	indirection = malloc(sizeof(size_t) * MAXSIZE);
	timerTest();

	for (size_t size = 1<<6; size <= MAXSIZE; size <<= 1) {
//		size_t size = 1<<23;
                initValues(A, indirection, size);

		t = get_clock();
		for (size_t i = 0; i < REPS/size; i++) {
			for (size_t j = 0; j < size; j++) {
				A[j] += A[indirection[j]];
			}
		}
		t = get_clock() - t;

		printf("Time for size %lu: %lfs. Per access: %lfns\n",
		       size, t,  (t*1.0E9)/REPS);
	}
}

void timerTest(void)
{
	double times[10];

	for (size_t i = 0; i < 10; i++) {
		times[i] = get_clock();
	}

	for (size_t i = 0; i < 10; i++) {
		printf("Time at %lu: %.10fs\n", i, times[i] - times[0]);
	}
}
