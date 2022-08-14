
#line 2 "level1.c"
/* Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/* This file defines "level 1" kernel functions. 
based mainly on level1.h */

#include <intrin.h>
#include <stdint.h>
#include <assert.h>
#include "pari.h"
#include "paripriv.h"
#include "int.h"
#include "parisys.h"

/* external */
int bfffo(ulong x);
int64_t divll(ulong x, ulong y, ulong *hiremainder);

/* based on src\kernel\none\addll.h */
/* note: __extension__ does not exist in visual studio C */

//ulong overflow = 0;
//ulong hiremainder = 0;
#define LOCAL_HIREMAINDER ulong hiremainder=0
#define LOCAL_OVERFLOW ulong overflow = 0;


/*
#define addll(a, b)                                             \
__extension__ ({                                                \
   ulong __arg1 = (a), __arg2 = (b), __value = __arg1 + __arg2; \
   overflow = (__value < __arg1);                               \
   __value;                                                     \
})
*/
int64_t addll(ulong a, ulong b, ulong *overflow) {
    ulong  __value = a + b;
    *overflow = (__value < a);
    return __value;
}

//#define addllx(a, b)                                          \
//__extension__({                                              \
//   ulong __arg1 = (a), __arg2 = (b), __value, __tmp = __arg1 + overflow;\
//   overflow = (__tmp < __arg1);                               \
//   __value = __tmp + __arg2;                                  \
//   overflow |= (__value < __tmp);                             \
//   __value;                                                   \
//})
int64_t addllx(ulong a, ulong b, ulong *overflow) {
    ulong __tmp = a + *overflow;
    *overflow = (__tmp < a);
    ulong  __value = __tmp + b;
    *overflow |= (__value < __tmp);
    return __value;
}

//#define subll(a, b)                                           \
//__extension__ ({                                              \
//   ulong __arg1 = (a), __arg2 = (b);                          \
//   overflow = (__arg2 > __arg1);                              \
//   __arg1 - __arg2;                                           \
//})
int64_t subll(ulong __arg1, ulong __arg2, ulong * overflow) {
    *overflow = (__arg2 > __arg1);
    return __arg1 - __arg2;
}


//#define subllx(a, b)                                  \
//__extension__({                                      \
//   ulong __arg1 = (a), __arg2 = (b), __value, __tmp = __arg1 - overflow;\
//   overflow = (__arg1 < overflow);                    \
//   __value = __tmp - __arg2;                          \
//   overflow |= (__arg2 > __tmp);                      \
//   __value;                                           \
//})
int64_t subllx(ulong a, ulong b, ulong* overflow) {
    ulong __tmp = a - *overflow;
    *overflow = (a < *overflow);
    ulong __value = __tmp - b;
    *overflow |= (b > __tmp);
    return __value;
}

/* end*/



/* return low 64 bits of x*y.
if product exceeds 64 bits top bits are in hiremainder,
otherwise hiremainder is set to zero */

#ifdef _WIN32
#pragma intrinsic(_umul128)
/* same function as original mulll but uses intrinsic to increase speed*/
int64_t mulll(ulong x, ulong y, ulong *hiremainder) {
    ulong result = _umul128(x, y, hiremainder);
    return result;
}
#else
/* taken from src64\kernel\none\mulll.h */
int64_t
mulll(ulong x, ulong y)
{
    const ulong xlo = LOWWORD(x), xhi = HIGHWORD(x);
    const ulong ylo = LOWWORD(y), yhi = HIGHWORD(y);
    ulong xylo, xymid, xyhi, xymidhi, xymidlo;
    ulong xhl, yhl;

    xylo = xlo * ylo;
    xyhi = xhi * yhi;
    xhl = xhi + xlo;
    yhl = yhi + ylo;
    xymid = xhl * yhl - (xyhi + xylo);

    xymidhi = HIGHWORD(xymid);
    xymidlo = xymid << BITS_IN_HALFULONG;

    xylo += xymidlo;
    hiremainder = xyhi + xymidhi + (xylo < xymidlo)
        + ((((xhl + yhl) >> 1) - xymidhi) & HIGHMASK);

    return xylo;
}
#endif

/* return x*y + old hiremainder 
if result exceeds 64 bits top bits are returned in hiremainder,
otherwise hiremainder is set to zero*/
#ifdef _WIN32
#pragma intrinsic(_umul128, _addcarry_u64)
/* _addcarry_u32(), _addcarry_u64()
Computes sum of two 32/64 bit wide unsigned integer values and a carry-in and 
returns the value of carry-out produced by the sum
Syntax:
extern unsigned char _addcarry_u64(unsigned char c_in, unsigned __int64 src1, unsigned __int64 src2, unsigned __int64 *sum_out);
Parameters
c_in     Value used for determining carry-in value
src1     32/64 bit source integer
src2     32/64 bit source integer
*sum_out Pointer to memory location where result is stored
*/
int64_t addmul(ulong x, ulong y, ulong *hiremainder) {
    ulong result, top;
    unsigned char c2;
    result = _umul128(x, y, &top);
    /* add hiremainder*/
    c2 = _addcarry_u64(0, result, *hiremainder, &result);
    *hiremainder = top + c2;
    assert(*hiremainder >= top); /* check for overflow*/
    return result;
}
#else
int64_t
addmul(ulong x, ulong y)
{
    const ulong xlo = LOWWORD(x), xhi = HIGHWORD(x);
    const ulong ylo = LOWWORD(y), yhi = HIGHWORD(y);
    ulong xylo, xymid, xyhi, xymidhi, xymidlo;
    ulong xhl, yhl;

    xylo = xlo * ylo; 
    xyhi = xhi * yhi;
    xhl = xhi + xlo; 
    yhl = yhi + ylo;
    xymid = xhl * yhl - (xyhi + xylo);

    /* add current value of hiremainder to xyhi, xylo*/
    xylo += hiremainder; 
    xyhi += (xylo < hiremainder);

    xymidhi = HIGHWORD(xymid);
    xymidlo = xymid << BITS_IN_HALFULONG;

    xylo += xymidlo;
    hiremainder = xyhi + xymidhi + (xylo < xymidlo)
        + ((((xhl + yhl) >> 1) - xymidhi) & HIGHMASK);

    return xylo;
}
#endif
/* end */


int64_t evallg(int64_t x) {
    if (x & ~LGBITS) 
        pari_err_OVERFLOW("lg()");
    return _evallg(x);
}

int64_t evalvalp(int64_t x) {
    int64_t v = _evalvalp(x);
    if (v & ~VALPBITS) 
        pari_err_OVERFLOW("valp()");
    return v;
}

int64_t evalexpo(int64_t x) {
    int64_t v = _evalexpo(x);
    if (v & ~EXPOBITS) 
        pari_err_OVERFLOW("expo()");
    return v;
}

int64_t evalprecp(int64_t x) {
    int64_t v = _evalprecp(x);
    if (x & ~((1ULL << (BITS_IN_LONG - VALPnumBITS)) - 1)) 
        pari_err_OVERFLOW("precp()");
    return v;
}

int varncmp(int64_t x, int64_t y) {
    if (varpriority[x] < varpriority[y]) 
        return  1;
    if (varpriority[x] > varpriority[y]) 
        return -1;
    return 0;
}

int64_t varnmin(int64_t x, int64_t y) {
    return (varpriority[x] <= varpriority[y]) ? x : y;
}

int64_t varnmax(int64_t x, int64_t y) {
    return (varpriority[x] >= varpriority[y]) ? x : y;
}

/* Inhibit some area gerepile-wise: declare it to be a non recursive
 * type, of length l. Thus gerepile won't inspect the zone, just copy it.
 * For the following situation:
 *   z = cgetg(t,a); av = avma; garbage(); ltop = avma;
 *   for (i=1; i<HUGE; i++) gel(z,i) = blah();
 *   stackdummy(av,ltop);
 * loses (av-ltop) words but save a costly gerepile. */
void stackdummy(pari_sp av, pari_sp ltop) {
    int64_t l = ((GEN)av) - ((GEN)ltop);
    if (l > 0) {
        GEN z = (GEN)ltop;
        z[0] = evaltyp(t_VECSMALL) | evallg(l);
#ifdef DEBUG
        { int64_t i; for (i = 1; i < l; i++) z[i] = 0; }
#endif
    }
}

void fixlg(GEN x, int64_t ly) {
    int64_t lx = lg(x), l = lx - ly;
    if (l > 0)
    { /* stackdummy(x+lx, x+ly) */
        GEN z = x + ly;
        z[0] = evaltyp(t_VECSMALL) | evallg(l);
        setlg(x, ly);
#ifdef DEBUG
        { int64_t i; for (i = 1; i < l; i++) z[i] = 0; }
#endif
    }
}

/* update lg(z) before affrr(y, z)  [ to cater for precision loss ]*/
void affrr_fixlg(GEN y, GEN z) { 
    fixlg(z, lg(y)); affrr(y, z); 
}

/*******************************************************************/
/*                                                                 */
/*                       ALLOCATE ON STACK                         */
/*                                                                 */
/*******************************************************************/
void set_avma(ulong av) { 
    avma = av; }

double gc_double(pari_sp av, double d) { 
    set_avma(av); return d; 
}

int64_t gc_long(pari_sp av, int64_t s) { 
    set_avma(av); return s; 
}

ulong gc_ulong(pari_sp av, ulong s) { 
    set_avma(av); 
    return s; 
}

int gc_bool(pari_sp av, int s) { 
    set_avma(av); 
    return s; 
}

int gc_int(pari_sp av, int s) { 
    set_avma(av); 
    return s; 
}

GEN gc_NULL(pari_sp av) { 
    set_avma(av); 
    return NULL; 
}

GEN gc_const(pari_sp av, GEN x) { 
    set_avma(av); 
    return x; 
}

GEN new_chunk(size_t x) /* x is a number of longs */
{
    GEN z = ((GEN)avma) - x;
    CHECK_CTRLC
        if (x > (avma - pari_mainstack->bot) / sizeof(int64_t))
            new_chunk_resize(x);
    set_avma((pari_sp)z);
#ifdef MEMSTEP
    if (DEBUGMEM > 1 && pari_mainstack->memused != DISABLE_MEMUSED) {
        int64_t d = (int64_t)pari_mainstack->memused - (int64_t)z;
        if (labs(d) > 4 * MEMSTEP)
        {
            pari_mainstack->memused = (pari_sp)z;
            err_printf("...%4.0lf Mbytes used\n",
                (pari_mainstack->top - pari_mainstack->memused) / 1048576.);
        }
    }
#endif
    return z;
}

char* stack_malloc(size_t N) {
    int64_t n = nchar2nlong(N);
    return (char*)new_chunk(n);
}

char* stack_malloc_align(size_t N, int64_t k) {
    ulong d = ((ulong)avma) % k, e = ((ulong)N) % k;
    if (d) 
        (void)new_chunk(d / sizeof(int64_t));
    if (e) 
        N += k - e;
    return (char*)new_chunk(nchar2nlong(N));
}

char* stack_calloc(size_t N) {
    char* p = stack_malloc(N);
    memset(p, 0, N); 
    return p;
}

char* stack_calloc_align(size_t N, int64_t k) {
    ulong d = ((ulong)avma) % k, e = ((ulong)N) % k;
    if (d) 
        (void)new_chunk(d / sizeof(int64_t));
    if (e) 
        N += k - e;
    return stack_calloc(N);
}

/* cgetg(lg(x), typ(x)), set *lx. Implicit unsetisclone() */
GEN cgetg_copy(GEN x, int64_t* plx) {
    GEN y;
    *plx = lg(x); y = new_chunk((size_t)*plx);
    y[0] = x[0] & (TYPBITS | LGBITS); return y;
}

GEN cgetg_block(int64_t x, int64_t y) {
    GEN z = newblock((size_t)x);
    z[0] = CLONEBIT | evaltyp(y) | evallg(x);
    return z;
}

GEN cgetg(int64_t x, int64_t y) {
    GEN z = new_chunk((size_t)x);
    z[0] = evaltyp(y) | evallg(x);
    return z;
}

GEN cgeti(int64_t x) {
    GEN z = new_chunk((size_t)x);
    z[0] = evaltyp(t_INT) | evallg(x);
    return z;
}

GEN cgetipos(int64_t x) {
    GEN z = cgeti(x);
    z[1] = evalsigne(1) | evallgefint(x);
    return z;
}

GEN cgetineg(int64_t x) {
    GEN z = cgeti(x);
    z[1] = evalsigne(-1) | evallgefint(x);
    return z;
}

GEN cgetr_block(int64_t x) {
    GEN z = newblock((size_t)x);
    z[0] = CLONEBIT | evaltyp(t_REAL) | evallg(x);
    return z;
}

GEN cgetr(int64_t x) {
    GEN z = new_chunk((size_t)x);
    z[0] = evaltyp(t_REAL) | evallg(x);
    return z;
}

/*******************************************************************/
/*                                                                 */
/*                     COPY, NEGATION, ABSOLUTE VALUE              */
/*                                                                 */
/*******************************************************************/
/* cannot do memcpy because sometimes x and y overlap */
GEN leafcopy(GEN x) {
    int64_t lx = lg(x);
    GEN y = new_chunk(lx); /* can't use cgetg_copy, in case x,y overlap */
    while (--lx > 0) 
        y[lx] = x[lx];
    y[0] = x[0] & (TYPBITS | LGBITS); 
    return y;
}

GEN icopy(GEN x) {
    int64_t i = lgefint(x), lx = i;
    GEN y = new_chunk(lx); /* can't use cgeti, in case x,y overlap */
    while (--i > 0) 
        y[i] = x[i];
    y[0] = evaltyp(t_INT) | evallg(lx);
    return y;
}

GEN icopyspec(GEN x, int64_t nx) {
    int64_t i = nx + 2, lx = i;
    GEN y = new_chunk(lx); /* can't use cgeti, in case x,y overlap */
    x -= 2; 
    while (--i >= 2) 
        y[i] = x[i];
    y[1] = evalsigne(1) | evallgefint(lx);
    y[0] = evaltyp(t_INT) | evallg(lx);
    return y;
}

GEN rcopy(GEN x) { 
    return leafcopy(x); 
}

GEN mpcopy(GEN x) { 
    return leafcopy(x); 
}

GEN mpabs(GEN x) { 
    GEN y = leafcopy(x); 
    setabssign(y); 
    return y; 
}

GEN mpabs_shallow(GEN x) { 
    return signe(x) < 0 ? mpabs(x) : x; 
}

GEN absi(GEN x) { 
    return mpabs(x); 
}
GEN absi_shallow(GEN x) { 
    return signe(x) < 0 ? negi(x) : x; 
}

GEN absr(GEN x) { 
    return mpabs(x); 
}

GEN mpneg(GEN x) { 
    GEN y = leafcopy(x); 
    togglesign(y); 
    return y; 
}

GEN negi(GEN x) { 
    return mpneg(x); 
}

GEN negr(GEN x) { 
    return mpneg(x); 
}

/* negate in place */
void togglesign(GEN x) { 
    if (x[1] & SIGNBITS) { 
        x[1] ^= HIGHBIT; 
    } 
} 

void setabssign(GEN x) { 
    x[1] &= ~HIGHBIT; 
}

/* negate in place, except universal constants */
void togglesign_safe(GEN* px) {
    switch (*px - gen_1) /* gen_1, gen_2, gen_m1, gen_m2 */
    {
    case 0: *px = gen_m1; break;
    case 3: *px = gen_m2;  break;
    case 6: *px = gen_1; break;
    case 9: *px = gen_2;  break;

    default: togglesign(*px);
    }
}

/* setsigne(y, signe(x)) */
void affectsign(GEN x, GEN y) {
    y[1] = (x[1] & SIGNBITS) | (y[1] & ~SIGNBITS);
}

/* copies sign in place, except for universal constants */
void affectsign_safe(GEN x, GEN* py) {
    if (((*py)[1] ^ x[1]) & HIGHBIT) 
        togglesign_safe(py);
}

/*******************************************************************/
/*                                                                 */
/*                     GEN -> LONG, LONG -> GEN                    */
/*                                                                 */
/*******************************************************************/
/* assume x != 0, return -x as a t_INT */
GEN utoineg(ulong x) { 
    GEN y = cgetineg(3); 
    y[2] = x; 
    return y; 
}

/* assume x != 0, return utoi(x) */
GEN utoipos(ulong x) { 
    GEN y = cgetipos(3); 
    y[2] = x; 
    return y; 
}

GEN utoi(ulong x) { 
    return x ? utoipos(x) : gen_0; 
}

GEN stoi(int64_t x) {
    if (!x) 
        return gen_0;
    return x > 0 ? utoipos((ulong)x) : utoineg((ulong)-x);
}

/* x 2^BIL + y */
GEN uutoi(ulong x, ulong y) {
    GEN z;
    if (!x) 
        return utoi(y);
    z = cgetipos(4);
    *int_W_lg(z, 1, 4) = x;
    *int_W_lg(z, 0, 4) = y; 
    return z;
}

/* - (x 2^BIL + y) */
GEN uutoineg(ulong x, ulong y) {
    GEN z;
    if (!x) 
        return y ? utoineg(y) : gen_0;
    z = cgetineg(4);
    *int_W_lg(z, 1, 4) = x;
    *int_W_lg(z, 0, 4) = y; 
    return z;
}

int64_t itos(GEN x) {
    int64_t s = signe(x);
    int64_t u;

    if (!s) 
        return 0;
    u = x[2];
    if (lgefint(x) > 3 || u < 0)
        pari_err_OVERFLOW("t_INT-->int64_t assignment");
    return (s > 0) ? u : -u;
}

/* as itos, but return 0 if too large. Cf is_bigint */
int64_t itos_or_0(GEN x) {
    int64_t n;
    if (lgefint(x) != 3 || (n = x[2]) & HIGHBIT) 
        return 0;
    return signe(x) > 0 ? n : -n;
}

ulong itou(GEN x) {
    switch (lgefint(x)) {
    case 2: return 0;
    case 3: return x[2];
    default:
        pari_err_OVERFLOW("t_INT-->ulong assignment");
        return 0; /* LCOV_EXCL_LINE */
    }
}

/* as itou, but return 0 if too large. Cf is_bigint */
ulong itou_or_0(GEN x) {
    if (lgefint(x) != 3) 
        return 0;
    return (ulong)x[2];
}

ulong umuluu_or_0(ulong x, ulong y) {
    ulong z;
    LOCAL_HIREMAINDER;
    z = mulll(x, y, &hiremainder);
    return hiremainder ? 0 : z;
}

/* return x*y if <= n, else 0. Beware overflow */
ulong umuluu_le(ulong x, ulong y, ulong n) {
    ulong z;
    LOCAL_HIREMAINDER;
    z = mulll(x, y, &hiremainder);
    return (hiremainder || z > n) ? 0 : z;
}

GEN real_0_bit(int64_t bitprec) { 
    GEN x = cgetr(2); 
    x[1] = evalexpo(bitprec); 
    return x; 
}

GEN real_0(int64_t prec) { 
    return real_0_bit(-prec2nbits(prec)); 
}

GEN real_1_bit(int64_t bit) { 
    return real_1(nbits2prec(bit)); 
}

GEN real_1(int64_t prec) {
    GEN x = cgetr(prec);
    int64_t i;
    x[1] = evalsigne(1) | _evalexpo(0);
    x[2] = (int64_t)HIGHBIT; for (i = 3; i < prec; i++) x[i] = 0;
    return x;
}

GEN real_m1(int64_t prec) {
    GEN x = cgetr(prec);
    int64_t i;
    x[1] = evalsigne(-1) | _evalexpo(0);
    x[2] = (int64_t)HIGHBIT; for (i = 3; i < prec; i++) x[i] = 0;
    return x;
}

/* 2.^n */
GEN real2n(int64_t n, int64_t prec) { 
    GEN z = real_1(prec); 
    setexpo(z, n); 
    return z; 
}

GEN real_m2n(int64_t n, int64_t prec) { 
    GEN z = real_m1(prec); 
    setexpo(z, n); 
    return z; 
}

GEN stor(int64_t s, int64_t prec) { 
    GEN z = cgetr(prec); 
    affsr(s, z); 
    return z; 
}

GEN utor(ulong s, int64_t prec) { 
    GEN z = cgetr(prec); 
    affur(s, z); 
    return z; 
}

GEN itor(GEN x, int64_t prec) { 
    GEN z = cgetr(prec); 
    affir(x, z); 
    return z; 
}

GEN rtor(GEN x, int64_t prec) { 
    GEN z = cgetr(prec); 
    affrr(x, z); 
    return z; 
}

ulong int_bit(GEN x, int64_t n) {
    int64_t r, q = dvmdsBIL(n, &r);
    return q < lgefint(x) - 2 ? ((ulong)*int_W(x, q) >> r) & 1ULL : 0;
}

/*******************************************************************/
/*                                                                 */
/*                           COMPARISON                            */
/*                                                                 */
/*******************************************************************/
int cmpss(int64_t a, int64_t b) {
    return a > b ? 1 : (a < b ? -1 : 0);
}

int cmpuu(ulong a, ulong b) {
    return a > b ? 1 : (a < b ? -1 : 0);
}

int cmpir(GEN x, GEN y) {
    pari_sp av;
    GEN z;

    if (!signe(x)) return -signe(y);
    if (!signe(y))
    {
        if (expo(y) >= expi(x)) 
            return 0;
        return signe(x);
    }
    av = avma; 
    z = itor(x, realprec(y)); 
    set_avma(av);
    return cmprr(z, y); /* cmprr does no memory adjustment */
}

int cmpri(GEN x, GEN y) { 
    return -cmpir(y, x); 
}

int cmpsr(int64_t x, GEN y) {
    pari_sp av;
    GEN z;

    if (!x) 
        return -signe(y);
    av = avma; 
    z = stor(x, LOWDEFAULTPREC); 
    set_avma(av);
    return cmprr(z, y);
}

int cmprs(GEN x, int64_t y) { 
    return -cmpsr(y, x); 
}

/* compare x and y */
int cmpui(ulong x, GEN y) {
    ulong p;
    if (!x) 
        return -signe(y);
    if (signe(y) <= 0) 
        return 1;
    if (lgefint(y) > 3) 
        return -1;
    p = y[2]; 
    if (p == x) 
        return 0;
    return p < x ? 1 : -1;
}

int cmpiu(GEN x, ulong y) { 
    return -cmpui(y, x); 
}

/* compare x and |y| */
int abscmpui(ulong x, GEN y) {
    int64_t l = lgefint(y);
    ulong p;

    if (!x) 
        return (l > 2) ? -1 : 0;
    if (l == 2) 
        return 1;
    if (l > 3) 
        return -1;
    p = y[2]; 
    if (p == x) 
        return 0;
    return p < x ? 1 : -1;
}

int abscmpiu(GEN x, ulong y) { 
    return -abscmpui(y, x); 
}

int cmpsi(int64_t x, GEN y) {
    ulong p;

    if (!x) 
        return -signe(y);

    if (x > 0) {
        if (signe(y) <= 0) 
            return 1;
        if (lgefint(y) > 3) 
            return -1;
        p = y[2]; 
        if (p == (ulong)x) 
            return 0;
        return p < (ulong)x ? 1 : -1;
    }

    if (signe(y) >= 0) 
        return -1;
    if (lgefint(y) > 3) 
        return 1;
    p = y[2]; 
    if (p == (ulong)-x) 
        return 0;
    return p < (ulong)(-x) ? -1 : 1;
}

int cmpis(GEN x, int64_t y) { 
    return -cmpsi(y, x); 
}

int mpcmp(GEN x, GEN y) {
    if (typ(x) == t_INT)
        return (typ(y) == t_INT) ? cmpii(x, y) : cmpir(x, y);
    return (typ(y) == t_INT) ? -cmpir(y, x) : cmprr(x, y);
}

/* x == y ? */
int equalui(ulong x, GEN y) {
    if (!x) 
        return !signe(y);
    if (signe(y) <= 0 || lgefint(y) != 3) 
        return 0;
    return ((ulong)y[2] == (ulong)x);
}

/* x == y ? */
int equalsi(int64_t x, GEN y) {
    if (!x) 
        return !signe(y);
    if (x > 0)     {
        if (signe(y) <= 0 || lgefint(y) != 3) 
            return 0;
        return ((ulong)y[2] == (ulong)x);
    }
    if (signe(y) >= 0 || lgefint(y) != 3) 
        return 0;
    return ((ulong)y[2] == (ulong)-x);
}
/* x == |y| ? */
int absequalui(ulong x, GEN y) {
    if (!x) 
        return !signe(y);
    return (lgefint(y) == 3 && (ulong)y[2] == x);
}

int absequaliu(GEN x, ulong y) { 
    return absequalui(y, x); 
}

int equalis(GEN x, int64_t y) { 
    return equalsi(y, x); 
}

int equaliu(GEN x, ulong y) { 
    return equalui(y, x); 
}

/* assume x != 0, is |x| == 2^n ? */
int absrnz_equal2n(GEN x) {
    if ((ulong)x[2] == HIGHBIT)
    {
        int64_t i, lx = lg(x);
        for (i = 3; i < lx; i++)
            if (x[i]) return 0;
        return 1;
    }
    return 0;
}

/* assume x != 0, is |x| == 1 ? */
int absrnz_equal1(GEN x) { 
    return !expo(x) && absrnz_equal2n(x); 
}

int64_t maxss(int64_t x, int64_t y) { 
    return x > y ? x : y; 
}

int64_t minss(int64_t x, int64_t y) { 
    return x < y ? x : y; 
}

int64_t minuu(ulong x, ulong y) { 
    return x < y ? x : y; 
}

int64_t maxuu(ulong x, ulong y) { 
    return x > y ? x : y; 
}

double maxdd(double x, double y) { 
    return x > y ? x : y; 
}

double mindd(double x, double y) { 
    return x < y ? x : y; 
}

/*******************************************************************/
/*                                                                 */
/*                             ADD / SUB                           */
/*                                                                 */
/*******************************************************************/
GEN subuu(ulong x, ulong y) {
    ulong z;
    LOCAL_OVERFLOW;
    z = subll(x, y, &overflow);
    /*z = x - y;
    overflow = (ulong)y > (ulong)x;*/
    return overflow ? utoineg(-(int64_t)z) : utoi(z);
}

GEN adduu(ulong x, ulong y) { 
    ulong t = x + y; 
    return uutoi((t < x), t); 
}

GEN addss(int64_t x, int64_t y) {
    if (!x) 
        return stoi(y);
    if (!y) 
        return stoi(x);
    if (x > 0) 
        return y > 0 ? adduu(x, y) : subuu(x, -y);

    if (y > 0) 
        return subuu(y, -x);
    else { /* - adduu(-x, -y) */
        ulong t = (-x) + (-y); 
        return uutoineg((t < (ulong)(-x)), t);
    }
}
GEN subss(int64_t x, int64_t y) { 
    return addss(-y, x); 
}


GEN
subii(GEN x, GEN y) {
    if (x == y) 
        return gen_0; /* frequent with x = y = gen_0 */
    return addii_sign(x, signe(x), y, -signe(y));
}

GEN addii(GEN x, GEN y) { 
    return addii_sign(x, signe(x), y, signe(y)); 
}

GEN addrr(GEN x, GEN y) { 
    return addrr_sign(x, signe(x), y, signe(y)); 
}

GEN subrr(GEN x, GEN y) { 
    return addrr_sign(x, signe(x), y, -signe(y)); 
}

GEN addir(GEN x, GEN y) { 
    return addir_sign(x, signe(x), y, signe(y)); 
}

GEN subir(GEN x, GEN y) { 
    return addir_sign(x, signe(x), y, -signe(y)); 
}

GEN subri(GEN x, GEN y) { 
    return addir_sign(y, -signe(y), x, signe(x)); 
} 

GEN addsi(int64_t x, GEN y) { 
    return addsi_sign(x, y, signe(y)); 
}

GEN addui(ulong x, GEN y) { 
    return addui_sign(x, y, signe(y)); 
}

GEN subsi(int64_t x, GEN y) { 
    return addsi_sign(x, y, -signe(y)); 
}

GEN subui(ulong x, GEN y) { 
    return addui_sign(x, y, -signe(y)); 
}

/*******************************************************************/
/*                                                                 */
/*                           MOD, REM, DIV                         */
/*                                                                 */
/*******************************************************************/
ulong mod2BIL(GEN x) { return *int_LSW(x); }
int64_t mod64(GEN x) { return mod2BIL(x) & 63; }
int64_t mod32(GEN x) { return mod2BIL(x) & 31; }
int64_t mod16(GEN x) { return mod2BIL(x) & 15; }
int64_t mod8(GEN x) { return mod2BIL(x) & 7; }
int64_t mod4(GEN x) { return mod2BIL(x) & 3; }
int64_t mod2(GEN x) { return mod2BIL(x) & 1; }

int mpodd(GEN x) { 
    return signe(x) && mod2(x); 
}

/* x mod 2^n, n < BITS_IN_LONG */
ulong umodi2n(GEN x, int64_t n) {
    int64_t s = signe(x);
    const ulong _2n = 1ULL << n;
    ulong m;
    if (!s) return 0;
    m = *int_LSW(x) & (_2n - 1);
    if (s < 0 && m) m = _2n - m;
    return m;
}

ulong Mod64(GEN x) { return umodi2n(x, 6); }
ulong Mod32(GEN x) { return umodi2n(x, 5); }
ulong Mod16(GEN x) { return umodi2n(x, 4); }
ulong Mod8(GEN x) { return umodi2n(x, 3); }
ulong Mod4(GEN x) { return umodi2n(x, 2); }
ulong Mod2(GEN x) { return umodi2n(x, 1); }

GEN truedivii(GEN a, GEN b) { 
    return truedvmdii(a, b, NULL); 
}

GEN truedivis(GEN a, int64_t b) { 
    return truedvmdis(a, b, NULL); 
} 

GEN truedivsi(int64_t a, GEN b) { 
    return truedvmdsi(a, b, NULL); 
}

GEN divii(GEN a, GEN b) { 
    return dvmdii(a, b, NULL); 
}

GEN remii(GEN a, GEN b) { 
    return dvmdii(a, b, ONLY_REM); 
}

GEN divss(int64_t x, int64_t y) { 
    return stoi(x / y); 
}

GEN modss(int64_t x, int64_t y) { 
    return utoi(smodss(x, y)); 
}

GEN remss(int64_t x, int64_t y) { 
    return stoi(x % y); 
}

int64_t smodss(int64_t x, int64_t y) {
    int64_t r = x % y;
    return (r >= 0) ? r : labs(y) + r;
}

ulong umodsu(int64_t x, ulong y) {
    return x >= 0 ? x % y : Fl_neg((-x) % y, y);
}

int64_t sdivss_rem(int64_t x, int64_t y, int64_t* r) {
    int64_t q;
    LOCAL_HIREMAINDER;
    if (!y) 
        pari_err_INV("sdivss_rem", gen_0);
    hiremainder = 0; 
    q = divll((ulong)labs(x), (ulong)labs(y), &hiremainder);
    if (x < 0) { 
        hiremainder = -((int64_t)hiremainder); 
        q = -q; 
    }
    if (y < 0) 
        q = -q;
    *r = hiremainder; 
    return q;
}

GEN divss_rem(int64_t x, int64_t y, int64_t* r) { 
    return stoi(sdivss_rem(x, y, r)); 
}

ulong udivuu_rem(ulong x, ulong y, ulong* r) {
    if (!y) 
        pari_err_INV("udivuu_rem", gen_0);
    *r = x % y; 
    return x / y;
}

ulong ceildivuu(ulong a, ulong b) {
    ulong c = a / b;
    return (a % b) ? c + 1 : c;
}

ulong uabsdivui_rem(ulong x, GEN y, ulong* r) {
    int64_t q, s = signe(y);
    LOCAL_HIREMAINDER;

    if (!s) 
        pari_err_INV("uabsdivui_rem", gen_0);
    if (!x || lgefint(y) > 3) { 
        *r = x; 
        return 0; 
    }
    hiremainder = 0; 
    q = (int64_t)divll(x, (ulong)y[2], &hiremainder);
    if (s < 0) 
        q = -q;
    *r = hiremainder; return q;
}

/* assume d != 0 and |n| / d can be represented as an ulong.
 * Return |n|/d, set *r = |n| % d */
ulong uabsdiviu_rem(GEN n, ulong d, ulong* r) {
    switch (lgefint(n))
    {
    case 2: *r = 0; 
        return 0;
    case 3:
        {
            ulong nn = n[2];
            *r = nn % d; 
            return nn / d;
        }
    default: /* 4 */
        {
            ulong n1, n0, q;
            LOCAL_HIREMAINDER;
            n0 = *int_W(n, 0);
            n1 = *int_W(n, 1);
            hiremainder = n1;
            q = divll(n0, d, &hiremainder);
            *r = hiremainder; 
            return q;
        }
    }
}

int64_t sdivsi_rem(int64_t x, GEN y, int64_t* r) {
    int64_t q, s = signe(y);
    LOCAL_HIREMAINDER;

    if (!s) 
        pari_err_INV("sdivsi_rem", gen_0);
    if (!x || lgefint(y) > 3 || ((int64_t)y[2]) < 0) { 
        *r = x; 
        return 0; 
    }
    hiremainder = 0; 
    q = (int64_t)divll(labs(x), (ulong)y[2], &hiremainder);
    if (x < 0) { 
        hiremainder = -((int64_t)hiremainder); 
        q = -q; 
    }
    if (s < 0) 
        q = -q;
    *r = hiremainder; 
    return q;
}

GEN divsi_rem(int64_t s, GEN y, int64_t* r) { 
    return stoi(sdivsi_rem(s, y, r)); 
}

int64_t sdivsi(int64_t x, GEN y) {
    int64_t q, s = signe(y);

    if (!s) 
        pari_err_INV("sdivsi", gen_0);
    if (!x || lgefint(y) > 3 || ((int64_t)y[2]) < 0) 
        return 0;
    q = labs(x) / y[2];
    if (x < 0) 
        q = -q;
    if (s < 0) 
        q = -q;
    return q;
}

GEN dvmdss(int64_t x, int64_t y, GEN* z) {
    int64_t r;
    GEN q = divss_rem(x, y, &r);
    *z = stoi(r); return q;
}

int64_t dvmdsBIL(int64_t n, int64_t* r) { 
    *r = remsBIL(n); 
    return divsBIL(n); 
}

ulong dvmduBIL(ulong n, ulong* r) { 
    *r = remsBIL(n); 
    return divsBIL(n); 
}

GEN dvmdsi(int64_t x, GEN y, GEN* z) {
    int64_t r;
    GEN q = divsi_rem(x, y, &r);
    *z = stoi(r); return q;
}

GEN dvmdis(GEN x, int64_t y, GEN* z) {
    int64_t r;
    GEN q = divis_rem(x, y, &r);
    *z = stoi(r); return q;
}

int64_t smodis(GEN x, int64_t y) {
    pari_sp av = avma;
    int64_t r; (void)divis_rem(x, y, &r);
    return gc_long(av, (r >= 0) ? r : labs(y) + r);
}

GEN modis(GEN x, int64_t y) { 
    return stoi(smodis(x, y)); 
}

GEN modsi(int64_t x, GEN y) {
    int64_t r; 
    (void)sdivsi_rem(x, y, &r);
    return (r >= 0) ? stoi(r) : addsi_sign(r, y, 1);
}

ulong umodui(ulong x, GEN y) {
    if (!signe(y)) 
        pari_err_INV("umodui", gen_0);
    if (!x || lgefint(y) > 3) 
        return x;
    return x % (ulong)y[2];
}

ulong ugcdiu(GEN x, ulong y) { 
    return ugcd(umodiu(x, y), y); 
}

ulong ugcdui(ulong y, GEN x) { 
    return ugcd(umodiu(x, y), y); 
}

GEN remsi(int64_t x, GEN y) {
    int64_t r; 
    (void)sdivsi_rem(x, y, &r); 
    return stoi(r);
}

GEN remis(GEN x, int64_t y) {
    pari_sp av = avma;
    int64_t r;
    (void)divis_rem(x, y, &r); 
    set_avma(av); 
    return stoi(r);
}

GEN rdivis(GEN x, int64_t y, int64_t prec) {
    GEN z = cgetr(prec);
    pari_sp av = avma;
    affrr(divrs(itor(x, prec), y), z);
    set_avma(av); return z;
}

GEN rdivsi(int64_t x, GEN y, int64_t prec) {
    GEN z = cgetr(prec);
    pari_sp av = avma;
    affrr(divsr(x, itor(y, prec)), z);
    set_avma(av); return z;
}

GEN rdivss(int64_t x, int64_t y, int64_t prec) {
    GEN z = cgetr(prec);
    pari_sp av = avma;
    affrr(divrs(stor(x, prec), y), z);
    set_avma(av); return z;
}

void rdiviiz(GEN x, GEN y, GEN z) {
    int64_t prec = realprec(z), lx = lgefint(x), ly = lgefint(y);
    if (lx == 2) { 
        affur(0, z); 
        return; 
    }
    if (ly == 3) {
        affir(x, z); 
        if (signe(y) < 0) 
            togglesign(z);
        affrr(divru(z, y[2]), z);
    }
    else if (lx > prec + 1 || ly > prec + 1)  {
        affir(x, z); 
        affrr(divri(z, y), z);
    }
    else  {
        int64_t b = bit_accuracy(prec) + expi(y) - expi(x) + 1;
        GEN q = divii(b > 0 ? shifti(x, b) : x, y);
        affir(q, z); 
        if (b > 0) 
            shiftr_inplace(z, -b);
    }
    set_avma((ulong)z);
}

GEN rdivii(GEN x, GEN y, int64_t prec) {
    GEN z = cgetr(prec); 
    rdiviiz(x, y, z); 
    return z;
}

GEN fractor(GEN x, int64_t prec) {
    return rdivii(gel(x, 1), gel(x, 2), prec);
}

int dvdii(GEN x, GEN y) {
    pari_sp av = avma;
    GEN r;
    if (!signe(x)) 
        return 1;
    if (!signe(y)) 
        return 0;
    r = remii(x, y);
    return gc_bool(av, r == gen_0);
}

int dvdsi(int64_t x, GEN y) {
    if (x == 0) 
        return 1;
    if (!signe(y)) 
        return 0;
    if (lgefint(y) != 3) 
        return 0;
    return x % y[2] == 0;
}

int dvdui(ulong x, GEN y) {
    if (x == 0) 
        return 1;
    if (!signe(y)) 
        return 0;
    if (lgefint(y) != 3) 
        return 0;
    return x % y[2] == 0;
}

int dvdis(GEN x, int64_t y) {
    return y ? smodis(x, y) == 0 : signe(x) == 0;
}

int dvdiu(GEN x, ulong y) {
    return y ? umodiu(x, y) == 0 : signe(x) == 0;
}

int dvdisz(GEN x, int64_t y, GEN z) {
    const pari_sp av = avma;
    int64_t r;
    GEN p1 = divis_rem(x, y, &r);
    set_avma(av); if (r) return 0;
    affii(p1, z); return 1;
}

int dvdiuz(GEN x, ulong y, GEN z) {
    const pari_sp av = avma;
    ulong r;
    GEN p1 = absdiviu_rem(x, y, &r);
    set_avma(av); if (r) return 0;
    affii(p1, z); return 1;
}

int dvdiiz(GEN x, GEN y, GEN z) {
    const pari_sp av = avma;
    GEN p2, p1 = dvmdii(x, y, &p2);
    if (signe(p2)) return gc_bool(av, 0);
    affii(p1, z); return gc_bool(av, 1);
}

/* copied from divll_pre.h*/
ulong /* requires u1 <= n, n normalised */
remll_pre_normalized(ulong u1, ulong u0, ulong n, ulong ninv) {
    ulong q0, q1, r;
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;
    q0 = mulll(ninv, u1, &hiremainder);
    q1 = hiremainder;
    q0 = addll(q0, u0, &overflow);   
   /* overflow = (q0 + u0) < q0;
    q0 += u0;*/
    q1 = addllx(q1, u1, &overflow);  
   /* tmp = q1 + overflow;
    overflow = tmp < q1;
    overflow |= (tmp + u1) < tmp;
    q1 = tmp + u1;*/

    r = u0 - (q1 + 1) * n;
    if (r >= q0)
        r += n;
    return r < n ? r : r - n;
}

ulong /* reduce <a_hi, a_lo> mod n */
remll_pre(ulong a_hi, ulong a_lo, ulong n, ulong ninv) {
    int norm = bfffo(n);
    int bits = BITS_IN_LONG - norm;
    ulong sn = n << norm;
    if (a_hi >= n) /* reduce a_hi first */
    {
        const ulong u1 = norm ? a_hi >> bits : 0;
        const ulong u0 = a_hi << norm;
        a_hi = remll_pre_normalized(u1, u0, sn, ninv) >> norm;
    }
    /* now reduce <a_hi, a_lo> */
    {
        const ulong u1 = ((a_hi << norm) | (norm ? a_lo >> bits : 0));
        const ulong u0 = a_lo << norm;
        return remll_pre_normalized(u1, u0, sn, ninv) >> norm;
    }
}

ulong remlll_pre(ulong u2, ulong u1, ulong u0, ulong n, ulong ninv) {
    u1 = remll_pre(u2, u1, n, ninv);
    return remll_pre(u1, u0, n, ninv);
}

ulong Fl_sqr_pre(ulong a, ulong p, ulong pi) {
    ulong x;
    LOCAL_HIREMAINDER;
    x = mulll(a, a, &hiremainder);
    return remll_pre(hiremainder, x, p, pi);
}

ulong Fl_mul_pre(ulong a, ulong b, ulong p, ulong pi) {
    ulong x;
    LOCAL_HIREMAINDER;
    x = mulll(a, b, &hiremainder);
    return remll_pre(hiremainder, x, p, pi);
}

ulong Fl_addmul_pre(ulong y0, ulong x0, ulong x1, ulong p, ulong pi)
{
    ulong l0, h0;
    LOCAL_HIREMAINDER;
    hiremainder = y0;
    l0 = addmul(x0, x1, &hiremainder);
    h0 = hiremainder;
    return remll_pre(h0, l0, p, pi);
}

ulong Fl_addmulmul_pre(ulong x0, ulong y0, ulong x1, ulong y1, ulong p, ulong pi) {
    ulong l0, l1, h0, h1;
    ulong temp;
    LOCAL_OVERFLOW;
    LOCAL_HIREMAINDER;
    l0 = mulll(x0, y0, &hiremainder); 
    h0 = hiremainder;
    l1 = mulll(x1, y1, &hiremainder); 
    h1 = hiremainder;

    //l0 = addll(l0, l1);    /* expanded inline */
    overflow = (l0 + l1) < l1;
    l0 += l1;
    //h0 = addllx(h0, h1);   /* expanded inline */
    temp = h0 + overflow;
    overflow = temp < h0;
    overflow |= (temp + h1) < temp;
    h0 += h1;
    return overflow ? remlll_pre(1, h0, l0, p, pi) : remll_pre(h0, l0, p, pi);
}

ulong Fl_ellj_pre(ulong a4, ulong a6, ulong p, ulong pi)
{
    /* a43 = 4 a4^3 */
    ulong a43 = Fl_double(Fl_double(
        Fl_mul_pre(a4, Fl_sqr_pre(a4, p, pi), p, pi), p), p);
    /* a62 = 27 a6^2 */
    ulong a62 = Fl_mul_pre(Fl_sqr_pre(a6, p, pi), 27 % p, p, pi);
    ulong z1 = Fl_mul_pre(a43, 1728 % p, p, pi);
    ulong z2 = Fl_add(a43, a62, p);
    return Fl_div(z1, z2, p);
}

/*******************************************************************/
/*                                                                 */
/*                        MP (INT OR REAL)                         */
/*                                                                 */
/*******************************************************************/
GEN mptrunc(GEN x) { 
    return typ(x) == t_INT ? icopy(x) : truncr(x); 
}

GEN mpfloor(GEN x) { 
    return typ(x) == t_INT ? icopy(x) : floorr(x); 
}
GEN mpceil(GEN x) { 
    return typ(x) == t_INT ? icopy(x) : ceilr(x); 
}

GEN mpround(GEN x) { 
    return typ(x) == t_INT ? icopy(x) : roundr(x); 
}

int64_t mpexpo(GEN x) { 
    return typ(x) == t_INT ? expi(x) : expo(x); 
}

GEN mpadd(GEN x, GEN y) {
    if (typ(x) == t_INT)
        return (typ(y) == t_INT) ? addii(x, y) : addir(x, y);
    return (typ(y) == t_INT) ? addir(y, x) : addrr(x, y);
}

GEN mpsub(GEN x, GEN y) {
    if (typ(x) == t_INT)
        return (typ(y) == t_INT) ? subii(x, y) : subir(x, y);
    return (typ(y) == t_INT) ? subri(x, y) : subrr(x, y);
}

GEN mpmul(GEN x, GEN y) {
    if (typ(x) == t_INT)
        return (typ(y) == t_INT) ? mulii(x, y) : mulir(x, y);
    return (typ(y) == t_INT) ? mulir(y, x) : mulrr(x, y);
}

GEN mpsqr(GEN x) { 
    return (typ(x) == t_INT) ? sqri(x) : sqrr(x); 
}

GEN mpdiv(GEN x, GEN y) {
    if (typ(x) == t_INT)
        return (typ(y) == t_INT) ? divii(x, y) : divir(x, y);
    return (typ(y) == t_INT) ? divri(x, y) : divrr(x, y);
}

/*******************************************************************/
/*                                                                 */
/*                          Z/nZ, n ULONG                          */
/*                                                                 */
/*******************************************************************/
ulong Fl_double(ulong a, ulong p)
{
    ulong res = a << 1;
    return (res >= p || res < a) ? res - p : res;
}

ulong Fl_triple(ulong a, ulong p) {
    ulong res = a << 1;
    if (res >= p || res < a) res -= p;
    res += a;
    return (res >= p || res < a) ? res - p : res;
}

ulong Fl_halve(ulong a, ulong p) {
    ulong ap, ap2;
    if ((a & 1ULL) == 0) 
        return a >> 1;
    ap = a + p; ap2 = ap >> 1;
    return ap >= a ? ap2 : (ap2 | HIGHBIT);
}

ulong Fl_add(ulong a, ulong b, ulong p) {
    ulong res = a + b;
    return (res >= p || res < a) ? res - p : res;
}

ulong Fl_neg(ulong x, ulong p) { 
    return x ? p - x : 0; 
}

ulong Fl_sub(ulong a, ulong b, ulong p) {
    ulong res = a - b;
    return (res > a) ? res + p : res;
}

/* centerlift(u mod p) */
int64_t Fl_center(ulong u, ulong p, ulong ps2) { 
    return (int64_t)(u > ps2) ? u - p : u; 
}

/* return (a*b) mod p */
ulong Fl_mul(ulong a, ulong b, ulong p) {
    ulong x;
    LOCAL_HIREMAINDER;
    x = mulll(a, b, &hiremainder);
    if (!hiremainder) 
        return x % p;
    (void)divll(x, p, &hiremainder);
    return hiremainder;
}

/* return a^2 mod p */
ulong Fl_sqr(ulong a, ulong p) {
    ulong x;
    LOCAL_HIREMAINDER;
    x = mulll(a, a, &hiremainder);
    if (!hiremainder) 
        return x % p;
    (void)divll(x, p, &hiremainder);
    return hiremainder;
}

/* don't assume that p is prime: can't special case a = 0 */
ulong Fl_div(ulong a, ulong b, ulong p) {
    return Fl_mul(a, Fl_inv(b, p), p);
}

/*******************************************************************/
/*                                                                 */
/*        DEFINED FROM EXISTING ONE EXPLOITING COMMUTATIVITY       */
/*                                                                 */
/*******************************************************************/
GEN addri(GEN x, GEN y) { 
    return addir(y, x); 
}

GEN addis(GEN x, int64_t s) { 
    return addsi(s, x); 
}

GEN
addiu(GEN x, ulong s) { 
    return addui(s, x); 
}

GEN addrs(GEN x, int64_t s) { 
    return addsr(s, x); 
}

GEN subiu(GEN x, int64_t y) { 
    GEN z = subui(y, x); 
    togglesign(z); 
    return z; 
}

GEN subis(GEN x, int64_t y) { 
    return addsi(-y, x); 
}

GEN subrs(GEN x, int64_t y) { 
    return addsr(-y, x); 
}

GEN mulis(GEN x, int64_t s) { 
    return mulsi(s, x); 
}

GEN muliu(GEN x, ulong s) { 
    return mului(s, x); 
}

GEN mulru(GEN x, ulong s) { 
    return mulur(s, x); 
}

GEN mulri(GEN x, GEN s) { 
    return mulir(s, x); 
}

GEN mulrs(GEN x, int64_t s) { 
    return mulsr(s, x); 
}

/*******************************************************************/
/*                                                                 */
/*                  VALUATION, EXPONENT, SHIFTS                    */
/*                                                                 */
/*******************************************************************/
int64_t vali(GEN x) {
    int64_t i;
    GEN xp;

    if (!signe(x)) return -1;
    xp = int_LSW(x);
    for (i = 0; !*xp; i++) xp = int_nextW(xp);
    return vals(*xp) + i * BITS_IN_LONG;
}

/* assume x > 0 */
int64_t expu(ulong x) { 
    return (BITS_IN_LONG - 1) - (int64_t)bfffo(x); 
}

int64_t expi(GEN x) {
    const int64_t lx = lgefint(x);
    return lx == 2 ? -(int64_t)HIGHEXPOBIT 
                   : bit_accuracy(lx) - (int64_t)bfffo(*int_MSW(x)) - 1;
}

GEN shiftr(GEN x, int64_t n) {
    const int64_t e = evalexpo(expo(x) + n);
    const GEN y = rcopy(x);

    if (e & ~EXPOBITS) 
        pari_err_OVERFLOW("expo()");
    y[1] = (y[1] & ~EXPOBITS) | e; 
    return y;
}

GEN mpshift(GEN x, int64_t s) { 
    return (typ(x) == t_INT) ? shifti(x, s) : shiftr(x, s); 
}

/* FIXME: adapt/use mpn_[lr]shift instead */
/* z2[imin..imax] := z1[imin..imax].f shifted left sh bits
 * (feeding f from the right). Assume sh > 0 */
void shift_left(GEN z2, GEN z1, int64_t imin, int64_t imax, ulong f, ulong sh) {
    GEN sb = z1 + imin, se = z1 + imax, te = z2 + imax;
    ulong l, m = BITS_IN_LONG - sh, k = f >> m;
    while (se > sb) {
        l = *se--;
        *te-- = (l << sh) | k;
        k = l >> m;
    }
    *te = (((ulong)*se) << sh) | k;
}

/* z2[imin..imax] := f.z1[imin..imax-1] shifted right sh bits
 * (feeding f from the left). Assume sh > 0 */
void shift_right(GEN z2, GEN z1, int64_t imin, int64_t imax, ulong f, ulong sh) {
    GEN sb = z1 + imin, se = z1 + imax, tb = z2 + imin;
    ulong k, l = *sb++, m = BITS_IN_LONG - sh;
    *tb++ = (l >> sh) | (f << m);
    while (sb < se) {
        k = l << m;
        l = *sb++;
        *tb++ = (l >> sh) | k;
    }
}

/* Backward compatibility. Inefficient && unused */
ulong shiftl(ulong x, ulong y) {
    ulong hiremainder = x >> (BITS_IN_LONG - y); 
    return (x << y);
}

ulong shiftlr(ulong x, ulong y) {
    ulong hiremainder = x << (BITS_IN_LONG - y); 
    return (x >> y);
}

void shiftr_inplace(GEN z, int64_t d) {
    setexpo(z, expo(z) + d);
}

/*******************************************************************/
/*                                                                 */
/*                           ASSIGNMENT                            */
/*                                                                 */
/*******************************************************************/
void affii(GEN x, GEN y) {
    int64_t lx = lgefint(x);
    if (lg(y) < lx) pari_err_OVERFLOW("t_INT-->t_INT assignment");
    while (--lx) y[lx] = x[lx];
}

void affsi(int64_t s, GEN x) {
    if (!s) 
        x[1] = evalsigne(0) | evallgefint(2);
    else
    {
        if (s > 0) { 
            x[1] = evalsigne(1) | evallgefint(3); 
            x[2] = s; 
        }
        else { 
            x[1] = evalsigne(-1) | evallgefint(3); 
            x[2] = -s; }
    }
}

void affui(ulong u, GEN x) {
    if (!u) 
        x[1] = evalsigne(0) | evallgefint(2);
    else { 
        x[1] = evalsigne(1) | evallgefint(3); 
        x[2] = u; }
}

void affsr(int64_t x, GEN y) {
    int64_t sh, i, ly = lg(y);

    if (!x) {
        y[1] = evalexpo(-prec2nbits(ly));
        return;
    }
    if (x < 0) {
        x = -x; 
        sh = bfffo(x);
        y[1] = evalsigne(-1) | _evalexpo((BITS_IN_LONG - 1) - sh);
    }
    else {
        sh = bfffo(x);
        y[1] = evalsigne(1) | _evalexpo((BITS_IN_LONG - 1) - sh);
    }
    y[2] = ((ulong)x) << sh; 
    for (i = 3; i < ly; i++) 
        y[i] = 0;
}

void affur(ulong x, GEN y) {
    int64_t sh, i, ly = lg(y);

    if (!x) {
        y[1] = evalexpo(-prec2nbits(ly));
        return;
    }
    sh = bfffo(x);
    y[1] = evalsigne(1) | _evalexpo((BITS_IN_LONG - 1) - sh);
    y[2] = x << sh; 
    for (i = 3; i < ly; i++) 
        y[i] = 0;
}

void affiz(GEN x, GEN y) { 
    if (typ(y) == t_INT) 
        affii(x, y); 
    else 
        affir(x, y); 
}

void affsz(int64_t x, GEN y) { 
    if (typ(y) == t_INT) 
        affsi(x, y); 
    else 
        affsr(x, y); 
}

void mpaff(GEN x, GEN y) { 
    if (typ(x) == t_INT) 
        affiz(x, y); 
    else affrr(x, y); 
}

/*******************************************************************/
/*                                                                 */
/*                    OPERATION + ASSIGNMENT                       */
/*                                                                 */
/*******************************************************************/

void addiiz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(addii(x, y), z); 
    set_avma(av);
}

void addirz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(addir(x, y), z); 
    set_avma(av);
}

void addriz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(addri(x, y), z); 
    set_avma(av);
}

void addrrz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(addrr(x, y), z); 
    set_avma(av);
}

void addsiz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(addsi(s, y), z); 
    set_avma(av);
}

void addsrz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(addsr(s, y), z); 
    set_avma(av);
}

void addssz(int64_t s, int64_t y, GEN z) {
    pari_sp av = avma; 
    affii(addss(s, y), z); 
    set_avma(av);
}

void diviiz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(divii(x, y), z); 
    set_avma(av);
}

void divirz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(divir(x, y), z); 
    set_avma(av);
}

void divisz(GEN x, int64_t y, GEN z) {
    pari_sp av = avma; 
    affii(divis(x, y), z); 
    set_avma(av);
}

void divriz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(divri(x, y), z); 
    set_avma(av);
}

void divrrz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(divrr(x, y), z); 
    set_avma(av);
}

void divrsz(GEN y, int64_t s, GEN z) {
    pari_sp av = avma; 
    affrr(divrs(y, s), z); 
    set_avma(av);
}

void divsiz(int64_t x, GEN y, GEN z) {
    int64_t junk; 
    affsi(sdivsi_rem(x, y, &junk), z);
}

void divsrz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(divsr(s, y), z); 
    set_avma(av);
}

void divssz(int64_t x, int64_t y, GEN z) {
    affsi(x / y, z);
}

void modisz(GEN y, int64_t s, GEN z) {
    affsi(smodis(y, s), z);
}

void modsiz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(modsi(s, y), z); 
    set_avma(av);
}

void modssz(int64_t s, int64_t y, GEN z) {
    affsi(smodss(s, y), z);
}

void mpaddz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(mpadd(x, y), z); 
    set_avma(av);
}

void mpsubz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(mpsub(x, y), z); 
    set_avma(av);
}

void mpmulz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(mpmul(x, y), z); 
    set_avma(av);
}

void muliiz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(mulii(x, y), z); 
    set_avma(av);
}
void mulirz(GEN x, GEN y, GEN z)
{
    pari_sp av = avma; mpaff(mulir(x, y), z); set_avma(av);
}

void mulriz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(mulri(x, y), z); 
    set_avma(av);
}

void mulrrz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(mulrr(x, y), z); 
    set_avma(av);
}

void mulsiz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(mulsi(s, y), z); 
    set_avma(av);
}

void mulsrz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    mpaff(mulsr(s, y), z); 
    set_avma(av);
}

void mulssz(int64_t s, int64_t y, GEN z) {
    pari_sp av = avma; 
    affii(mulss(s, y), z); 
    set_avma(av);
}

void remiiz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(remii(x, y), z); 
    set_avma(av);
}

void remisz(GEN y, int64_t s, GEN z) {
    pari_sp av = avma; 
    affii(remis(y, s), z); 
    set_avma(av);
}

void remsiz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(remsi(s, y), z); 
    set_avma(av);
}

void remssz(int64_t s, int64_t y, GEN z) {
    pari_sp av = avma; 
    affii(remss(s, y), z); 
    set_avma(av);
}

void subiiz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(subii(x, y), z); 
    set_avma(av);
}

void subirz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(subir(x, y), z); 
    set_avma(av);
}

void subisz(GEN y, int64_t s, GEN z) {
    pari_sp av = avma; 
    affii(addsi(-s, y), z); 
    set_avma(av);
}

void subriz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(subri(x, y), z); 
    set_avma(av);
}

void subrrz(GEN x, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(subrr(x, y), z); 
    set_avma(av);
}

void subrsz(GEN y, int64_t s, GEN z) {
    pari_sp av = avma; 
    affrr(addsr(-s, y), z); 
    set_avma(av);
}

void subsiz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affii(subsi(s, y), z); 
    set_avma(av);
}

void subsrz(int64_t s, GEN y, GEN z) {
    pari_sp av = avma; 
    affrr(subsr(s, y), z); 
    set_avma(av);
}

void subssz(int64_t x, int64_t y, GEN z) { 
    addssz(x, -y, z); 
}

void dvmdssz(int64_t x, int64_t y, GEN z, GEN t) {
    pari_sp av = avma;
    int64_t r;
    affii(divss_rem(x, y, &r), z); set_avma(av); affsi(r, t);
}

void dvmdsiz(int64_t x, GEN y, GEN z, GEN t) {
    pari_sp av = avma;
    int64_t r;
    affii(divsi_rem(x, y, &r), z); set_avma(av); affsi(r, t);
}

void dvmdisz(GEN x, int64_t y, GEN z, GEN t) {
    pari_sp av = avma;
    int64_t r;
    affii(divis_rem(x, y, &r), z); set_avma(av); affsi(r, t);
}

void dvmdiiz(GEN x, GEN y, GEN z, GEN t) {
    pari_sp av = avma;
    GEN r;
    affii(dvmdii(x, y, &r), z); affii(r, t); set_avma(av);
}
