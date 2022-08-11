#include "pari.h"
#line 2 "..bfffo.c"
//typedef unsigned long long ulong;
//#  define BITS_IN_LONG 64

/* could use _BitScanReverse64 intrinsic Note that bit numbering 
is reverse of bfffo but should still be quicker */
/* return bit number of most significant bit. The most significant bit is 0,
 the least significant is 63, but if x = 0, 64 is returned. */
int bfffo(ulong x)
{
	static int tabshi[16] = { 4,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0 };
	int value = BITS_IN_LONG - 4;
	ulong arg1 = x;
	if (arg1 & ~0xffffffffULL) { 
		value -= 32; 
		arg1 >>= 32; }
	if (arg1 & ~0xffffULL) { 
		value -= 16; 
		arg1 >>= 16; }
	if (arg1 & ~0x00ffULL) { 
		value -= 8; 
		arg1 >>= 8; }
	if (arg1 & ~0x000fULL) { 
		value -= 4; 
		arg1 >>= 4; }
	return value + tabshi[arg1];
}
#include <intrin.h>
#include <assert.h>
int bfffo2(ulong x) {
	unsigned long index;
	unsigned char dst;
	if (x == 0)
		return 64;
	dst = _BitScanReverse64(&index, x);
	assert(dst == 1);
	return (63 - index);
}