#include <stddef.h>

void init(double *restrict B, double *restrict C, size_t N);
void check(const double *restrict A, const double *restrict B, size_t N);
void clear(double *restrict A, size_t N);

static inline size_t min(size_t lhs, size_t rhs)
{
	return lhs < rhs ? lhs : rhs;
}

