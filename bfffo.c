//#include "pari.h"
#line 2 "..bfffo.c"
typedef unsigned long long ulong;
#  define BITS_IN_LONG 64

/* return bit number of most significant bit */
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