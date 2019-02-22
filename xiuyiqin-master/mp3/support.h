#include <stddef.h>

void init(float *restrict B, float *restrict C, size_t N);
void check(const float *restrict A, const float *restrict B, size_t N);
void clear(float *restrict A, size_t N);

static inline size_t triangular(size_t i)
{
	return i * (i+1) / 2;
}
static inline size_t min(size_t lhs, size_t rhs)
{
	return lhs < rhs ? lhs : rhs;
}

static inline size_t max(size_t lhs, size_t rhs)
{
	return lhs > rhs ? lhs : rhs;
}
