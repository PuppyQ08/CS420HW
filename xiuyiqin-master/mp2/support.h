#include <math.h>
#include <stddef.h>

int root_printf(const char *restrict format, ...);
void init(double *restrict B, double *restrict C, size_t N);
void check(const double *restrict A, const double *restrict B, size_t N);
void clear(double *restrict A, size_t N);

static inline size_t triangular(size_t N)
{
	return N * (N+1) / 2;
}

static inline size_t f1(int i, int P, size_t N)
{
	return (N * i) / P;
}

static inline size_t f2(int i, int P, size_t N)
{
	return N * cbrt((double)i / (double)P);
}

static inline size_t min(size_t lhs, size_t rhs)
{
	return lhs < rhs ? lhs : rhs;
}

static inline size_t max(size_t lhs, size_t rhs)
{
	return lhs > rhs ? lhs : rhs;
}
