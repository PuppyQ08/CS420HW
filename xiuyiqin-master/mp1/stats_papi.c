#include <err.h>
#include <omp.h>
#include <papi.h>
#include <stdio.h>
#include <stdlib.h>

#include "stats.h"

#define NUM_EVENTS 2
int events[NUM_EVENTS] = {PAPI_L2_DCM, PAPI_L2_DCA};

int eventset = PAPI_NULL;
long long local_nsec;
#pragma omp threadprivate(eventset, local_nsec)
long long nsec;
long long counters[NUM_EVENTS];

void init_stats(void)
{
	int retval = PAPI_library_init(PAPI_VER_CURRENT);// PAPI checks the argument against the internal value of PAPI_VER_CURRENT when the library was compiled. This guards against portability problems when updating the PAPI shared libraries on your system.
	if (retval != PAPI_VER_CURRENT && retval > 0) {
// PAPI checks the argument against the internal value of PAPI_VER_CURRENT when the library was compiled. This guards against portability problems when updating the PAPI shared libraries on your system.
		errx(EXIT_FAILURE,
		     "PAPI_library_init: version mismatch, %d != %d",
		     PAPI_VER_CURRENT, retval);
	} else if (retval < 0) {
		PAPI_perror("PAPI_library_init");//perror is to convert the error code to error message.
		exit(EXIT_FAILURE);
	}

	retval = PAPI_thread_init((unsigned long (*)(void))omp_get_thread_num);
	if (retval != PAPI_OK) {
		PAPI_perror("PAPI_thread_init");
		exit(EXIT_FAILURE);
	}

	for (size_t i = 0; i < NUM_EVENTS; i++) {
		counters[i] = 0;
	}

	nsec = 0;
}

void start_stats(void)
{
	if (PAPI_create_eventset(&eventset) != PAPI_OK) {
		PAPI_perror("PAPI_create_eventset");
		exit(EXIT_FAILURE);
	}

	if (PAPI_add_events(eventset, events, NUM_EVENTS) != PAPI_OK) {
		PAPI_perror("PAPI_add_events");
		exit(EXIT_FAILURE);
	}

	if (PAPI_start(eventset) != PAPI_OK) {
		PAPI_perror("PAPI_start");
		exit(EXIT_FAILURE);
	}

	local_nsec = PAPI_get_virt_nsec();
}

void stop_stats(void)
{
	long long local_counters[NUM_EVENTS];

	local_nsec = PAPI_get_virt_nsec() - local_nsec;

	if (PAPI_stop(eventset, local_counters) != PAPI_OK) {
		PAPI_perror("PAPI_stop");
		exit(EXIT_FAILURE);
	}

	#pragma omp critical
	if (local_nsec > nsec) {
		nsec = local_nsec;
	}

	for (size_t i = 0; i < NUM_EVENTS; i++) {
		#pragma omp atomic update
		counters[i] += local_counters[i];
	}

	if (PAPI_cleanup_eventset(eventset) != PAPI_OK) {
		PAPI_perror("PAPI_cleanup_eventset");
		exit(EXIT_FAILURE);
	}

	if (PAPI_destroy_eventset(&eventset) != PAPI_OK) {
		PAPI_perror("PAPI_destroy_eventset");
		exit(EXIT_FAILURE);
	}

	if (omp_get_thread_num() != 0) {
		if (PAPI_unregister_thread() != PAPI_OK) {
			PAPI_perror("PAPI_unregister_thread");
			exit(EXIT_FAILURE);
		}
	}
}

void print_stats(void)
{
	printf("PAPI event 0: %lld\n", counters[0]);
	printf("PAPI event 1: %lld\n", counters[1]);
	printf("Time: %fs\n", ((double)nsec) * 1e-9);

	for (size_t i = 0; i < NUM_EVENTS; i++) {
		counters[i] = 0;
	}

	nsec = 0;
}
