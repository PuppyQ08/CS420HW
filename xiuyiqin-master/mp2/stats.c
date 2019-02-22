#include <mpi.h>
#include <stdio.h>

#include "stats.h"
#include "support.h"

double sec;

void start_stats(void)
{
	sec = MPI_Wtime();
}

void stop_stats(void)
{
	sec = MPI_Wtime() - sec;
}

void print_stats(void)
{
	root_printf("Time: %fs\n", sec);
	sec = 0.0;
}
