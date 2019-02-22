#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "stats.h"

static inline double get_clock(void);

double local_sec;
#pragma omp threadprivate(local_sec)
double sec;

void init_stats(void) {
	sec = 0.0;
}

void start_stats(void)
{
	local_sec = get_clock();
}

void stop_stats(void)
{
	local_sec = get_clock() - local_sec;

	#pragma omp critical
	if (local_sec > sec)
		sec = local_sec;
}

void print_stats(void)
{
	printf("Time: %fs\n", sec);
	sec = 0.0;
}

static inline double get_clock(void)
{
	struct timespec tp;
	int ok = clock_gettime(CLOCK_MONOTONIC, &tp);
	if (ok != 0) {
		err(EXIT_FAILURE, "clock_gettime");
	}

	return (tp.tv_sec*1.0 + tp.tv_nsec*1.0E-9);
}

