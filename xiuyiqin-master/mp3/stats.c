#include <omp.h>
#include <stdio.h>

#include "stats.h"

double local_sec;
#pragma omp threadprivate(local_sec)
double sec;

void init_stats(void) {
	sec = 0.0;
}

void start_stats(void)
{
	local_sec = omp_get_wtime();
}

void stop_stats(void)
{
	local_sec = omp_get_wtime() - local_sec;

	#pragma omp critical
	if (local_sec > sec)
		sec = local_sec;
}

void print_stats(void)
{
	printf("Time: %fs\n", sec);
	sec = 0.0;
}

