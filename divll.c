#include "pari.h"
#include "int.h"
#include "parisys.h"
#include <intrin.h>

/* external */
int bfffo(ulong x);
extern ulong hiremainder;


#define __GLUE(hi, lo) (((hi) << BITS_IN_HALFULONG) | (lo))
#define __SPLIT(a, b, c) b = HIGHWORD(a); c = LOWWORD(a)
#define __LDIV(a, b, q, r) q = a / b; r = a - q*b

/* divide (hiremainder * 2^BITS_IN_LONG + n0) by d; assume hiremainder < d.
 * Return quotient, set hiremainder to remainder */

#ifdef _WIN32
/* The _udiv128 intrinsic divides a 128-bit unsigned integer by a 64-bit 
unsigned integer. The return value holds the quotient, and the intrinsic 
returns the remainder through a pointer parameter. 
*/
int64_t divll(ulong n0, ulong d, ulong* hiremainder) {
    /* top 64 bits of dividend are in hiremainder, bottom 64 bits in n0,
    divisor is d. quotient won't overflow unless initial value of hiremainder >= d.
    remainder is returned in hiremainder. */
    ulong quotient = _udiv128(*hiremainder, n0, d, hiremainder);
    return quotient;
}
#else
int64_t divll(ulong n0, ulong d) {
    
    ulong __d1, __d0, __q1, __q0, __r1, __r0, __m, __n1, __n0;
    ulong __k, __d;

    __n1 = hiremainder; 
    __n0 = n0; 
    __d = d;

    if (__n1 == 0) { /* Only one division needed */
        __LDIV(__n0, __d, __q1, hiremainder);
    }
    else if (__d < LOWMASK) { /* Two half-word divisions  */
        __n1 = __GLUE(__n1, HIGHWORD(__n0));
        __LDIV(__n1, __d, __q1, __r1);
        __n1 = __GLUE(__r1, LOWWORD(__n0));
        __LDIV(__n1, __d, __q0, hiremainder);
        __q1 = __GLUE(__q1, __q0);
    }

    else { /* General case */
        if (__d & HIGHBIT) {
            __k = 0; 
            __SPLIT(__d, __d1, __d0);
        }
        else {
            __k = bfffo(__d);
            __n1 = (__n1 << __k) | (__n0 >> (BITS_IN_LONG - __k));
            __n0 = __n0 << __k;
            __d = __d << __k; 
            __SPLIT(__d, __d1, __d0);
        }

        __LDIV(__n1, __d1, __q1, __r1);
        __m = __q1 * __d0;
        __r1 = __GLUE(__r1, HIGHWORD(__n0));

        if (__r1 < __m) {
            __q1--, __r1 += __d;
            if (__r1 >= __d) /* we didn't get carry when adding to __r1 */
                if (__r1 < __m) __q1--, __r1 += __d;
        }

        __r1 -= __m;
        __LDIV(__r1, __d1, __q0, __r0);
        __m = __q0 * __d0;
        __r0 = __GLUE(__r0, LOWWORD(__n0));

        if (__r0 < __m) {
            __q0--, __r0 += __d;
            if (__r0 >= __d)
                if (__r0 < __m) __q0--, __r0 += __d;
        }
        hiremainder = (__r0 - __m) >> __k;
        __q1 = __GLUE(__q1, __q0);
    }
    return __q1;
}
#endif
