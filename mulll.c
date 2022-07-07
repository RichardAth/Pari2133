#include <intrin.h>
#include <stdint.h>

#pragma intrinsic(_mul128)

/* return a*b  The top 64 bits of 128-bit product are in c */
int64_t mulll(int64_t a, int64_t b, int64_t* c) {
	int64_t r = _mul128(a, b, c);
	return r;
}