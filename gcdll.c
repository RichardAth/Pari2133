//#line 2 "../src/kernel/none/gcdll.c"
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
#include "pari.h"
#include "paripriv.h"
#include "int.h"

#define LOCAL_OVERFLOW ulong overflow = 0
#define LOCAL_HIREMAINDER ulong hiremainder=0

/***********************************************************************/
/**                                                                   **/
/**                          GCD                                      **/
/**                                                                   **/
/***********************************************************************/
/* Fast ulong gcd.  Called with y odd; x can be arbitrary (but will most of
 * the time be smaller than y). */

/* Gotos are Harmful, and Programming is a Science.  E.W.Dijkstra. */
INLINE ulong
gcduodd(ulong x, ulong y)         /* assume y&1==1, y > 1 */
{
  if (!x) return y;
  /* fix up x */
  while (!(x&1)) x>>=1;
  if (x==1) return 1;
  if (x==y) return y;
  else if (x>y) goto xislarger;/* will be rare, given how we'll use this */
  /* loop invariants: x,y odd and distinct. */
 yislarger:
  if ((x^y)&2)                 /* ...01, ...11 or vice versa */
    y=(x >> 2)+(y >> 2)+1;         /* ==(x+y) >> 2 except it can't overflow */
  else                         /* ...01,...01 or ...11,...11 */
    y=(y-x) >> 2;                /* now y!=0 in either case */
  while (!(y&1)) y>>=1;        /* kill any windfall-gained powers of 2 */
  if (y==1) return 1;          /* comparand == return value... */
  if (x==y) return y;          /* this and the next is just one comparison */
  else if (x<y) goto yislarger;/* else fall through to xislarger */

 xislarger:                    /* same as above, seen through a mirror */
  if ((x^y)&2)
    x=(x >> 2)+(y >> 2)+1;
  else
    x=(x-y) >> 2;                /* x!=0 */
  while (!(x&1)) x>>=1;
  if (x==1) return 1;
  if (x==y) return y;
  else if (x>y) goto xislarger;

  goto yislarger;
}
/* Gotos are useful, and Programming is an Art.  D.E.Knuth. */
/* PS: Of course written with Dijkstra's lessons firmly in mind... --GN */

/* at least one of a or b is odd, return gcd(a,b) */
INLINE ulong
mygcduodd(ulong a, ulong b)
{
  ulong c;
  if (b&1)
  {
    if (a==1 || b==1)
      c = 1;
    else
     c = gcduodd(a, b);
  }
  else
  {
    if (a==1)
      c = 1;
    else
      c = gcduodd(b, a);
  }
  return c;
}

/* modified right shift binary algorithm with at most one division */
ulong
ugcd(ulong a,ulong b)
{
  int64_t v;
  if (!b) return a;
  if (!a) return b;
  if (a>b) { a %= b; if (!a) return b; }
  else     { b %= a; if (!b) return a; }
  v = vals(a|b);
  return mygcduodd(a >> v, b >> v) << v;
}
int64_t
cgcd(int64_t a,int64_t b) { return (int64_t)ugcd(labs(a), labs(b)); }

/* For gcdii(): assume a>b>0, return gcd(a,b) as a GEN */
GEN
igcduu(ulong a, ulong b)
{
  int64_t v;
  a %= b; if (!a) return utoipos(b);
  v = vals(a|b);
  return utoipos( mygcduodd(a >> v, b >> v) << v );
}

/*Warning: overflows silently if lcm does not fit*/
ulong
ulcm(ulong a, ulong b)
{
  ulong d = ugcd(a,b);
  if (!d) return 0;
  return d == 1? a*b: a*(b/d);
}
int64_t
clcm(int64_t a,int64_t b) { return ulcm(labs(a), labs(b)); }

/********************************************************************/
/**                                                                **/
/**               INTEGER EXTENDED GCD  (AND INVMOD)               **/
/**                                                                **/
/********************************************************************/
/* Two basic ideas - (1) avoid many integer divisions, especially when the
 * quotient is 1 which happens ~ 40% of the time.  (2) Use Lehmer's trick as
 * modified by Jebelean of extracting a couple of words' worth of leading bits
 * from both operands, and compute partial quotients from them as long as we
 * can be sure of their values.  Jebelean's modifications consist in
 * inequalities from which we can quickly decide whether to carry on or to
 * return to the outer loop, and in re-shifting after the first word's worth of
 * bits has been used up.  All of this is described in R. Lercier's thesis
 * [pp148-153 & 163f.], except his outer loop isn't quite right: the catch-up
 * divisions needed when one partial quotient is larger than a word are missing.
 *
 * The API consists of invmod() and bezout() below; the single-word routines
 * xgcduu and xxgcduu may be called directly if desired; lgcdii() probably
 * doesn't make much sense out of context.
 *
 * The whole lot is a factor 6 .. 8 faster on word-sized operands, and asym-
 * ptotically about a factor 2.5 .. 3, depending on processor architecture,
 * than the naive continued-division code.  Unfortunately, thanks to the
 * unrolled loops and all, the code is lengthy. */

/*==================================
 * xgcduu(d,d1,f,v,v1,s)
 * xxgcduu(d,d1,f,u,u1,v,v1,s)
 * rgcduu(d,d1,vmax,u,u1,v,v1,s)
 *==================================*/
/*
 * Fast `final' extended gcd algorithm, acting on two ulongs.  Ideally this
 * should be replaced with assembler versions wherever possible.  The present
 * code essentially does `subtract, compare, and possibly divide' at each step,
 * which is reasonable when hardware division (a) exists, (b) is a bit slowish
 * and (c) does not depend a lot on the operand values (as on i486).  When
 * wordsize division is in fact an assembler routine based on subtraction,
 * this strategy may not be the most efficient one.
 *
 * xxgcduu() should be called with  d > d1 > 0, returns gcd(d,d1), and assigns
 * the usual signless cont.frac. recurrence matrix to [u, u1; v, v1]  (i.e.,
 * the product of all the [0, 1; 1 q_j] where the leftmost factor arises from
 * the quotient of the first division step),  and the information about the
 * implied signs to s  (-1 when an odd number of divisions has been done,
 * 1 otherwise).  xgcduu() is exactly the same except that u,u1 are not com-
 * puted (and not returned, of course).
 *
 * The input flag f should be set to 1 if we know in advance that gcd(d,d1)==1
 * (so we can stop the chain division one step early:  as soon as the remainder
 * equals 1).  Use this when you intend to use only what would be v, or only
 * what would be u and v, after that final division step, but not u1 and v1.
 * With the flag in force and thus without that final step, the interesting
 * quantity/ies will still sit in [u1 and] v1, of course.
 *
 * For computing the inverse of a single-word INTMOD known to exist, pass f=1
 * to xgcduu(), and obtain the result from s and v1.  (The routine does the
 * right thing when d1==1 already.)  For finishing a multiword modinv known
 * to exist, pass f=1 to xxgcduu(), and multiply the returned matrix  (with
 * properly adjusted signs)  onto the values v' and v1' previously obtained
 * from the multiword division steps.  Actually, just take the scalar product
 * of [v',v1'] with [u1,-v1], and change the sign if s==-1.  (If the final
 * step had been carried out, it would be [-u,v], and s would also change.)
 * For reducing a rational number to lowest terms, pass f=0 to xgcduu().
 * Finally, f=0 with xxgcduu() is useful for Bezout computations.
 * (It is safe for invmod() to call xgcduu() with f=1, because f&1 doesn't
 * make a difference when gcd(d,d1)>1.  The speedup is negligible.)
 *
 * In principle, when gcd(d,d1) is known to be 1, it is straightforward to
 * recover the final u,u1 given only v,v1 and s.  However, it probably isn't
 * worthwhile, as it trades a few multiplications for a division.
 *
 * rgcduu() is a variant of xxgcduu() which does not have f  (the effect is
 * that of f=0),  but instead has a ulong vmax parameter, for use in rational
 * reconstruction below.  It returns when v1 exceeds vmax;  v will never
 * exceed vmax.  (vmax=0 is taken as a synonym of ULONG_MAX i.e. unlimited,
 * in which case rgcduu behaves exactly like xxgcduu with f=0.)  The return
 * value of rgcduu() is typically meaningless;  the interesting part is the
 * matrix. */

ulong
xgcduu(ulong d, ulong d1, int f, ulong* v, ulong* v1, int64_t *s)
{
  ulong xv,xv1, xs, q,res;
  LOCAL_HIREMAINDER;

  /* The above blurb contained a lie.  The main loop always stops when d1
   * has become equal to 1.  If (d1 == 1 && !(f&1)) after the loop, we do
   * the final `division' of d by 1 `by hand' as it were.
   *
   * The loop has already been unrolled once.  Aggressive optimization could
   * well lead to a totally unrolled assembler version.
   *
   * On modern x86 architectures, this loop is a pig anyway.  The division
   * instruction always puts its result into the same pair of registers, and
   * we always want to use one of them straight away, so pipeline performance
   * will suck big time.  An assembler version should probably do a first loop
   * computing and storing all the quotients -- their number is bounded in
   * advance -- and then assembling the matrix in a second pass.  On other
   * architectures where we can cycle through four or so groups of registers
   * and exploit a fast ALU result-to-operand feedback path, this is much less
   * of an issue. */
  xs = res = 0;
  xv = 0ULL; xv1 = 1ULL;
  while (d1 > 1ULL)
  {
    d -= d1; /* no need to use subll */
    if (d >= d1)
    {
      hiremainder = 0; 
      q = 1 + divll(d,d1, &hiremainder); 
      d = hiremainder;
      xv += q * xv1;
    }
    else
      xv += xv1;
    if (d <= 1ULL) { 
        xs=1; 
        break; /* possible loop exit */
    } 
    /* repeat with inverted roles */
    d1 -= d;
    if (d1 >= d)     {
      hiremainder = 0; 
      q = 1 + divll(d1,d, &hiremainder); 
      d1 = hiremainder;
      xv1 += q * xv;
    }
    else
      xv1 += xv;
  }

  if (!(f&1))
  { /* division by 1 postprocessing if needed */
    if (xs && d==1) { 
        xv1 += d1 * xv; 
        xs = 0; 
        res = 1ULL; }
    else 
        if (!xs && d1==1) { 
            xv += d * xv1; 
            xs = 1; 
            res = 1ULL; 
        }
  }

  if (xs) {
    *s = -1; 
    *v = xv1; 
    *v1 = xv;
    return (res ? res : (d==1 ? 1ULL : d1));
  }
  else {
    *s = 1; 
    *v = xv; 
    *v1 = xv1;
    return (res ? res : (d1==1 ? 1ULL : d));
  }
}

ulong xxgcduu(ulong d, ulong d1, int f, ulong* u, ulong* u1,
        ulong* v, ulong* v1, int64_t *s)
{
  ulong xu,xu1, xv,xv1, xs, q,res;
  LOCAL_HIREMAINDER;

  xs = res = 0;
  xu = xv1 = 1ULL;
  xu1 = xv = 0ULL;
  while (d1 > 1ULL)
  {
    /* no need to use subll */
    d -= d1;
    if (d >= d1)
    {
      hiremainder = 0; 
      q = 1 + divll(d,d1, &hiremainder); 
      d = hiremainder;
      xv += q * xv1;
      xu += q * xu1;
    }
    else { 
        xv += xv1; 
        xu += xu1; 
    }
    if (d <= 1ULL) { 
        xs=1; 
        break; /* possible loop exit */
    } 
    /* repeat with inverted roles */
    d1 -= d;
    if (d1 >= d)
    {
      hiremainder = 0; q = 1 + divll(d1,d, &hiremainder); d1 = hiremainder;
      xv1 += q * xv;
      xu1 += q * xu;
    }
    else
    { xv1 += xv; xu1 += xu; }
  }

  if (!(f&1))
  { /* division by 1 postprocessing if needed */
    if (xs && d==1)
    {
      xv1 += d1 * xv;
      xu1 += d1 * xu;
      xs = 0; res = 1ULL;
    }
    else if (!xs && d1==1)
    {
      xv += d * xv1;
      xu += d * xu1;
      xs = 1; res = 1ULL;
    }
  }

  if (xs)
  {
    *s = -1; *u = xu1; *u1 = xu; *v = xv1; *v1 = xv;
    return (res ? res : (d==1 ? 1ULL : d1));
  }
  else
  {
    *s = 1; *u = xu; *u1 = xu1; *v = xv; *v1 = xv1;
    return (res ? res : (d1==1 ? 1ULL : d));
  }
}

ulong
rgcduu(ulong d, ulong d1, ulong vmax,
       ulong* u, ulong* u1, ulong* v, ulong* v1, int64_t *s)
{
  ulong xu,xu1, xv,xv1, xs, q, res=0;
  int f = 0;
  LOCAL_HIREMAINDER;

  if (vmax == 0) vmax = ULONG_MAX;
  xs = res = 0;
  xu = xv1 = 1ULL;
  xu1 = xv = 0ULL;
  while (d1 > 1ULL)
  {
    d -= d1; /* no need to use subll */
    if (d >= d1)
    {
      hiremainder = 0; 
      q = 1 + divll(d,d1, &hiremainder); 
      d = hiremainder;
      xv += q * xv1;
      xu += q * xu1;
    }
    else
    { xv += xv1; xu += xu1; }
    /* possible loop exit */
    if (xv > vmax) { 
        xs=f=1; 
        break; 
    }
    if (d <= 1ULL) { 
        xs=1; 
        break; 
    }
    /* repeat with inverted roles */
    d1 -= d;
    if (d1 >= d)
    {
      hiremainder = 0; 
      q = 1 + divll(d1,d, &hiremainder); 
      d1 = hiremainder;
      xv1 += q * xv;
      xu1 += q * xu;
    }
    else { 
        xv1 += xv; 
        xu1 += xu; 
    }
    /* possible loop exit */
    if (xv1 > vmax) { 
        f=1; 
        break; 
    }
  }

  if (!(f&1))
  { /* division by 1 postprocessing if needed */
    if (xs && d==1) {
      xv1 += d1 * xv;
      xu1 += d1 * xu;
      xs = 0; res = 1ULL;
    }
    else if (!xs && d1==1)  {
      xv += d * xv1;
      xu += d * xu1;
      xs = 1; 
      res = 1ULL;
    }
  }

  if (xs)  {
    *s = -1; 
    *u = xu1; 
    *u1 = xu; 
    *v = xv1; 
    *v1 = xv;
    return (res ? res : (d==1 ? 1ULL : d1));
  }
  else
  {
    *s = 1; *u = xu; *u1 = xu1; *v = xv; *v1 = xv1;
    return (res ? res : (d1==1 ? 1ULL : d));
  }
}

/*==================================
 * cbezout(a,b,uu,vv)
 *==================================
 * Same as bezout() but for C longs.
 *    Return g = gcd(a,b) >= 0, and assign longs u,v through pointers uu,vv
 *    such that g = u*a + v*b.
 * Special cases:
 *    a = b = 0 ==> pick u=1, v=0 (and return 1, surprisingly)
 *    a != 0 = b ==> keep v=0
 *    a = 0 != b ==> keep u=0
 *    |a| = |b| != 0 ==> keep u=0, set v=+-1
 * Assignments through uu,vv happen unconditionally. */
int64_t
cbezout(int64_t a,int64_t b,int64_t *uu,int64_t *vv)
{
  int64_t s,*t;
  ulong d = labs(a), d1 = labs(b);
  ulong r,u,u1,v,v1;

  if (!b)
  {
    *vv=0L;
    if (!a) { *uu=1L; return 0L; }
    *uu = a < 0 ? -1L : 1L;
    return (int64_t)d;
  }
  else if (!a || (d == d1))
  {
    *uu = 0L; *vv = b < 0 ? -1L : 1L;
    return (int64_t)d1;
  }
  else if (d == 1) /* frequently used by nfinit */
  {
    *uu = a; *vv = 0L;
    return 1L;
  }
  else if (d < d1)
  {
/* bug in gcc-2.95.3:
 * s = a; a = b; b = s; produces wrong result a = b. This is OK:  */
    { int64_t _x = a; a = b; b = _x; } /* in order to keep the right signs */
    r = d; d = d1; d1 = r;
    t = uu; uu = vv; vv = t;
  }
  /* d > d1 > 0 */
  r = xxgcduu(d, d1, 0, &u, &u1, &v, &v1, &s);
  if (s < 0)
  {
    *uu = a < 0 ? (int64_t)u : -(int64_t)u;
    *vv = b < 0 ? -(int64_t)v : (int64_t)v;
  }
  else
  {
    *uu = a < 0 ? -(int64_t)u : (int64_t)u;
    *vv = b < 0 ? (int64_t)v : -(int64_t)v;
  }
  return (int64_t)r;
}

/*==================================
 * lgcdii(d,d1,u,u1,v,v1,vmax)
 *==================================*/
/* Lehmer's partial extended gcd algorithm, acting on two t_INT GENs.
 *
 * Tries to determine, using the leading 2*BITS_IN_LONG significant bits of d
 * and a quantity of bits from d1 obtained by a shift of the same displacement,
 * as many partial quotients of d/d1 as possible, and assigns to [u,u1;v,v1]
 * the product of all the [0,1; 1,qj] thus obtained, where the leftmost
 * factor arises from the quotient of the first division step.
 *
 * For use in rational reconstruction, vmax can be given a nonzero value.
 * In this case, we will return early as soon as v1 > vmax (i.e. v <= vmax)
 *
 * MUST be called with  d > d1 > 0, and with  d occupying more than one
 * significant word.  Returns the number of reduction/swap steps carried out,
 * possibly zero, or under certain conditions minus that number.  When the
 * return value is nonzero, the caller should use the returned recurrence
 * matrix to update its own copies of d,d1.  When the return value is
 * nonpositive, and the latest remainder after updating turns out to be
 * nonzero, the caller should at once attempt a full division, rather than
 * trying lgcdii() again -- this typically happens when we are about to
 * encounter a quotient larger than half a word. (This is not detected
 * infallibly -- after a positive return value, it is possible that the next
 * stage will end up needing a full division.  After a negative return value,
 * however, this is certain, and should be acted upon.)
 *
 * The sign information, for which xgcduu() has its return argument s, is now
 * implicit in the LSB of our return value, and the caller may take advantage
 * of the fact that a return value of +-1 implies u==0,u1==v==1  [only v1 pro-
 * vides interesting information in this case].  One might also use the fact
 * that if the return value is +-2, then u==1, but this is rather marginal.
 *
 * If it was not possible to determine even the first quotient, either because
 * we're too close to an integer quotient or because the quotient would be
 * larger than one word  (if the `leading digit' of d1 after shifting is all
 * zeros), we return 0 and do not assign anything to the last four args.
 *
 * The division chain might even run to completion. It is up to the caller to
 * detect this case. This routine does not change d or d1; this is also up to
 * the caller */
int
lgcdii(ulong* d, ulong* d1, ulong* u, ulong* u1, ulong* v, ulong* v1,
       ulong vmax)
{
  /* Strategy:  (1) Extract/shift most significant bits.  We assume that d
   * has at least two significant words, but we can cope with a one-word d1.
   * Let dd,dd1 be the most significant dividend word and matching part of the
   * divisor.
   * (2) Check for overflow on the first division.  For our purposes, this
   * happens when the upper half of dd1 is zero.  (Actually this is detected
   * during extraction.)
   * (3) Get a fix on the first quotient.  We compute q = floor(dd/dd1), which
   * is an upper bound for floor(d/d1), and which gives the true value of the
   * latter if (and-almost-only-if) the remainder dd' = dd-q*dd1 is >= q.
   * (If it isn't, we give up.  This is annoying because the subsequent full
   * division will repeat some work already done, but it happens infrequently.
   * Doing the extra-bit-fetch in this case would be awkward.)
   * (4) Finish initializations.
   *
   * The remainder of the action is comparatively boring... The main loop has
   * been unrolled once (so we don't swap things and we can apply Jebelean's
   * termination conditions which alternatingly take two different forms during
   * successive iterations).  When we first run out of sufficient bits to form
   * a quotient, and have an extra word of each operand, we pull out two whole
   * word's worth of dividend bits, and divisor bits of matching significance;
   * to these we apply our partial matrix (disregarding overflow because the
   * result mod 2^(2*BITS_IN_LONG) will in fact give the correct values), and
   * re-extract one word's worth of the current dividend and a matching amount
   * of divisor bits.  The affair will normally terminate with matrix entries
   * just short of a whole word.  (We terminate the inner loop before these can
   * possibly overflow.) */
  ulong dd,dd1,ddlo,dd1lo, sh,shc; /* `digits', shift count */
  ulong xu,xu1, xv,xv1, q; /* recurrences, partial quotient, count */
  int res;
  ulong tmp0,tmp1,tmp2,tmpd,tmpu,tmpv; /* temps */
  ulong dm1, d1m1;
  int64_t ld, ld1, lz;
  int skip = 0;
  LOCAL_OVERFLOW;
  LOCAL_HIREMAINDER;

  /* following is just for convenience: vmax==0 means no bound */
  if (vmax == 0) vmax = ULONG_MAX;
  ld = lgefint(d); ld1 = lgefint(d1); lz = ld - ld1; /* >= 0 */
  if (lz > 1) return 0; /* rare */

  d = int_MSW(d);  dm1  = *int_precW(d);
  d1 = int_MSW(d1);d1m1 = *int_precW(d1);
  dd1lo = 0; /* unless we find something better */
  sh = bfffo(*d);

  if (sh)
  { /* do the shifting */
    shc = BITS_IN_LONG - sh;
    if (lz)
    { /* dividend longer than divisor */
      dd1 = (*d1 >> shc);
      if (!(HIGHMASK & dd1)) return 0; /* overflow detected */
      if (ld1 > 3)
        dd1lo = (*d1 << sh) + (d1m1 >> shc);
      else
        dd1lo = (*d1 << sh);
    }
    else
    { /* dividend and divisor have the same length */
      dd1 = (*d1 << sh);
      if (!(HIGHMASK & dd1)) return 0;
      if (ld1 > 3)
      {
        dd1 += (d1m1 >> shc);
        if (ld1 > 4)
          dd1lo = (d1m1 << sh) + (*int_precW(int_precW(d1)) >> shc);
        else
          dd1lo = (d1m1 << sh);
      }
    }
    /* following lines assume d to have 2 or more significant words */
    dd = (*d << sh) + (dm1 >> shc);
    if (ld > 4)
      ddlo = (dm1 << sh) + (*int_precW(int_precW(d)) >> shc);
    else
      ddlo = (dm1 << sh);
  }
  else
  { /* no shift needed */
    if (lz) return 0; /* dividend longer than divisor: overflow */
    dd1 = *d1;
    if (!(HIGHMASK & dd1)) return 0;
    if(ld1 > 3) dd1lo = d1m1;
    /* assume again that d has another significant word */
    dd = *d; ddlo = dm1;
  }
  /* First subtraction/division stage.  (If a subtraction initially suffices,
   * we don't divide at all.)  If a Jebelean condition is violated, and we
   * can't fix it even by looking at the low-order bits in ddlo,dd1lo, we
   * give up and ask for a full division.  Otherwise we commit the result,
   * possibly deciding to re-shift immediately afterwards. */
  dd -= dd1;
  if (dd < dd1)
  { /* first quotient known to be == 1 */
    xv1 = 1ULL;
    if (!dd) /* !(Jebelean condition), extraspecial case */
    { /* This actually happens. Now q==1 is known, but we underflow already.
       * OTOH we've just shortened d by a whole word. Thus we are happy and
       * return. */
      *u = 0; *v = *u1 = *v1 = 1ULL;
      return -1; /* Next step will be a full division. */
    }
  }
  else
  { /* division indicated */
    hiremainder = 0;
    xv1 = 1 + divll(dd, dd1, &hiremainder); /* xv1: alternative spelling of `q', here ;) */
    dd = hiremainder;
    if (dd < xv1) /* !(Jebelean cond'), nonextra special case */
    { /* Attempt to complete the division using the less significant bits,
       * before skipping right past the 1st loop to the reshift stage. */
      ddlo = subll(ddlo, mulll(xv1, dd1lo, &hiremainder), &overflow);
      dd = subllx(dd, hiremainder, &overflow);

      /* If we now have an overflow, q was too large. Thanks to our decision
       * not to get here unless the original dd1 had bits set in the upper half
       * of the word, we now know that the correct quotient is in fact q-1. */
      if (overflow)
      {
        xv1--;
        ddlo = addll(ddlo,dd1lo, &overflow);
        dd = addllx(dd,dd1, &overflow); /* overflows again which cancels the borrow */
        /* ...and fall through to skip=1 below */
      }
      else
      /* Test Jebelean condition anew, at this point using _all_ the extracted
       * bits we have.  This is clutching at straws; we have a more or less
       * even chance of succeeding this time.  Note that if we fail, we really
       * do not know whether the correct quotient would have been q or some
       * smaller value. */
        if (!dd && ddlo < xv1) return 0;

      /* Otherwise, we now know that q is correct, but we cannot go into the
       * 1st loop.  Raise a flag so we'll remember to skip past the loop.
       * Get here also after the q-1 adjustment case. */
      skip = 1;
    } /* if !(Jebelean) then */
  }
  res = 1;
  if (xv1 > vmax)
  { /* gone past the bound already */
    *u = 0ULL; *u1 = 1ULL; *v = 1ULL; *v1 = xv1;
    return res;
  }
  xu = 0ULL; xv = xu1 = 1ULL;

  /* Some invariants from here across the first loop:
   *
   * At this point, and again after we are finished with the first loop and
   * subsequent conditional, a division and the attached update of the
   * recurrence matrix have just been carried out completely.  The matrix
   * xu,xu1;xv,xv1 has been initialized (or updated, possibly with permuted
   * columns), and the current remainder == next divisor (dd at the moment)
   * is nonzero (it might be zero here, but then skip will have been set).
   *
   * After the first loop, or when skip is set already, it will also be the
   * case that there aren't sufficiently many bits to continue without re-
   * shifting.  If the divisor after reshifting is zero, or indeed if it
   * doesn't have more than half a word's worth of bits, we will have to
   * return at that point.  Otherwise, we proceed into the second loop.
   *
   * Furthermore, when we reach the re-shift stage, dd:ddlo and dd1:dd1lo will
   * already reflect the result of applying the current matrix to the old
   * ddorig:ddlo and dd1orig:dd1lo.  (For the first iteration above, this
   * was easy to achieve, and we didn't even need to peek into the (now
   * no longer existent!) saved words.  After the loop, we'll stop for a
   * moment to merge in the ddlo,dd1lo contributions.)
   *
   * Note that after the 1st division, even an a priori quotient of 1 cannot be
   * trusted until we've checked Jebelean's condition: it might be too small */
  if (!skip)
  {
    for(;;)
    { /* First half of loop divides dd into dd1, and leaves the recurrence
       * matrix xu,...,xv1 groomed the wrong way round (xu,xv will be the newer
       * entries) when successful. */
      tmpd = dd1 - dd;
      if (tmpd < dd)
      { /* quotient suspected to be 1 */
        tmpu = xu + xu1; /* cannot overflow -- everything bounded by
                          * the original dd during first loop */
        tmpv = xv + xv1;
      }
      else
      { /* division indicated */
        hiremainder = 0;
        q = 1 + divll(tmpd, dd, &hiremainder);
        tmpd = hiremainder;
        tmpu = xu + q*xu1; /* can't overflow, but may need to be undone */
        tmpv = xv + q*xv1;
      }

      tmp0 = addll(tmpv, xv1, &overflow);
      if ((tmpd < tmpu) || overflow ||
          (dd - tmpd < tmp0)) /* !(Jebelean cond.) */
        break; /* skip ahead to reshift stage */
      else
      { /* commit dd1, xu, xv */
        res++;
        dd1 = tmpd; xu = tmpu; xv = tmpv;
        if (xv > vmax) { 
            *u = xu1; 
            *u1 = xu; 
            *v = xv1; 
            *v1 = xv; 
            return res; 
        }
      }

      /* Second half of loop divides dd1 into dd, and the matrix returns to its
       * normal arrangement. */
      tmpd = dd - dd1;
      if (tmpd < dd1)
      { /* quotient suspected to be 1 */
        tmpu = xu1 + xu; /* cannot overflow */
        tmpv = xv1 + xv;
      }
      else
      { /* division indicated */
        hiremainder = 0;
        q = 1 + divll(tmpd, dd1, &hiremainder);
        tmpd = hiremainder;
        tmpu = xu1 + q*xu;
        tmpv = xv1 + q*xv;
      }

      tmp0 = addll(tmpu, xu, &overflow);
      if ((tmpd < tmpv) || overflow ||
          (dd1 - tmpd < tmp0)) /* !(Jebelean cond.) */
        break; /* skip ahead to reshift stage */
      else
      { /* commit dd, xu1, xv1 */
        res++;
        dd = tmpd; xu1 = tmpu; xv1 = tmpv;
        if (xv1 > vmax) { 
            *u = xu; 
            *u1 = xu1; 
            *v = xv; 
            *v1 = xv1; 
            return res; 
        }
      }
    } /* end of first loop */

    /* Intermezzo: update dd:ddlo, dd1:dd1lo.  (But not if skip is set.) */
    if (res&1)
    { /* after failed division in 1st half of loop:
       * [dd1:dd1lo,dd:ddlo] = [ddorig:ddlo,dd1orig:dd1lo]
       *                       * [ -xu, xu1 ; xv, -xv1 ]
       * Actually, we only multiply [ddlo,dd1lo] onto the matrix and add the
       * high-order remainders + overflows onto [dd1,dd] */
      tmp1 = mulll(ddlo, xu, &hiremainder); 
      tmp0 = hiremainder;
      tmp1 = subll(mulll(dd1lo,xv, &hiremainder), tmp1, &overflow);
      dd1 += subllx(hiremainder, tmp0, &overflow);
      tmp2 = mulll(ddlo, xu1, &hiremainder); 
      tmp0 = hiremainder;
      ddlo = subll(tmp2, mulll(dd1lo,xv1, &hiremainder), &overflow);
      dd += subllx(tmp0, hiremainder, &overflow);
      dd1lo = tmp1;
    }
    else
    { /* after failed division in 2nd half of loop:
       * [dd:ddlo,dd1:dd1lo] = [ddorig:ddlo,dd1orig:dd1lo]
       *                       * [ xu1, -xu ; -xv1, xv ]
       * Actually, we only multiply [ddlo,dd1lo] onto the matrix and add the
       * high-order remainders + overflows onto [dd,dd1] */
      tmp1 = mulll(ddlo, xu1, &hiremainder); 
      tmp0 = hiremainder;
      tmp1 = subll(tmp1, mulll(dd1lo,xv1, &hiremainder), &overflow);
      dd += subllx(tmp0, hiremainder, &overflow);
      tmp2 = mulll(ddlo, xu, &hiremainder); 
      tmp0 = hiremainder;
      dd1lo = subll(mulll(dd1lo,xv, &hiremainder), tmp2, &overflow);
      dd1 += subllx(hiremainder, tmp0, &overflow);
      ddlo = tmp1;
    }
  } /* end of skip-pable section:  get here also, with res==1, when there
     * was a problem immediately after the very first division. */

  /* Re-shift.  Note: the shift count _can_ be zero, viz. under the following
   * precise conditions:  The original dd1 had its topmost bit set, so the 1st
   * q was 1, and after subtraction, dd had its topmost bit unset.  If now dd=0,
   * we'd have taken the return exit already, so we couldn't have got here.
   * If not, then it must have been the second division which has gone amiss
   * (because dd1 was very close to an exact multiple of the remainder dd value,
   * so this will be very rare).  At this point, we'd have a fairly slim chance
   * of fixing things by re-examining dd1:dd1lo vs. dd:ddlo, but this is not
   * guaranteed to work. Instead of trying, we return at once and let caller
   * see to it that the initial subtraction is re-done usingall the bits of
   * both operands, which already helps, and the next round will either be a
   * full division  (if dd occupied a halfword or less), or another llgcdii()
   * first step.  In the latter case, since we try a little harder during our
   * first step, we may actually be able to fix the problem, and get here again
   * with improved low-order bits and with another step under our belt.
   * Otherwise we'll have given up above and forced a full division.
   *
   * If res is even, the shift count _cannot_ be zero.  (The first step forces
   * a zero into the remainder's MSB, and all subsequent remainders will have
   * inherited it.)
   *
   * The re-shift stage exits if the next divisor has at most half a word's
   * worth of bits.
   *
   * For didactic reasons, the second loop will be arranged in the same way
   * as the first -- beginning with the division of dd into dd1, as if res
   * was odd.  To cater for this, if res is actually even, we swap things
   * around during reshifting.  (During the second loop, the parity of res
   * does not matter; we know in which half of the loop we are when we decide
   * to return.) */
  if (res&1)
  { /* after odd number of division(s) */
    if (dd1 && (sh = bfffo(dd1)))
    {
      shc = BITS_IN_LONG - sh;
      dd = (ddlo >> shc) + (dd << sh);
      if (!(HIGHMASK & dd))
      {
        *u = xu; *u1 = xu1; *v = xv; *v1 = xv1;
        return -(int)res; /* full division asked for */
      }
      dd1 = (dd1lo >> shc) + (dd1 << sh);
    }
    else
    { /* time to return: <= 1 word left, or sh==0 */
      *u = xu; *u1 = xu1; *v = xv; *v1 = xv1;
      return res;
    }
  }
  else
  { /* after even number of divisions */
    if (dd)
    {
      sh = bfffo(dd); /* > 0 */
      shc = BITS_IN_LONG - sh;
      /* dd:ddlo will become the new dd1, and v.v. */
      tmpd = (ddlo >> shc) + (dd << sh);
      dd = (dd1lo >> shc) + (dd1 << sh);
      dd1 = tmpd;
      /* This completes the swap; now dd is again the current divisor */
      if (HIGHMASK & dd)
      { /* recurrence matrix is the wrong way round; swap it */
        tmp0 = xu; xu = xu1; xu1 = tmp0;
        tmp0 = xv; xv = xv1; xv1 = tmp0;
      }
      else
      { /* recurrence matrix is the wrong way round; fix this */
        *u = xu1; *u1 = xu; *v = xv1; *v1 = xv;
        return -(int)res;                /* full division asked for */
      }
    }
    else
    { /* time to return: <= 1 word left */
      *u = xu1; *u1 = xu; *v = xv1; *v1 = xv;
      return res;
    }
  } /* end reshift */

  /* The Second Loop.  Rip-off of the first, but we now check for overflow
   * in the recurrences.  Returns instead of breaking when we cannot fix the
   * quotient any longer. */
  for(;;)
  { /* First half of loop divides dd into dd1, and leaves the recurrence
     * matrix xu,...,xv1 groomed the wrong way round (xu,xv will be the newer
     * entries) when successful */
    tmpd = dd1 - dd;
    if (tmpd < dd)
    { /* quotient suspected to be 1 */
      tmpu = xu + xu1;
      tmpv = addll(xv, xv1, &overflow); /* xv,xv1 will overflow first */
      tmp1 = overflow;
    }
    else
    { /* division indicated */
      hiremainder = 0;
      q = 1 + divll(tmpd, dd, &hiremainder);
      tmpd = hiremainder;
      tmpu = xu + q*xu1;
      tmpv = addll(xv, mulll(q,xv1, &hiremainder), &overflow);
      tmp1 = overflow | hiremainder;
    }

    tmp0 = addll(tmpv, xv1, &overflow);
    if ((tmpd < tmpu) || overflow || tmp1 ||
        (dd - tmpd < tmp0)) /* !(Jebelean cond.) */
    { /* The recurrence matrix has not yet been warped... */
      *u = xu; *u1 = xu1; *v = xv; *v1 = xv1;
      break;
    }
    /* commit dd1, xu, xv */
    res++;
    dd1 = tmpd; xu = tmpu; xv = tmpv;
    if (xv > vmax) { 
        *u = xu1; 
        *u1 = xu; 
        *v = xv1; 
        *v1 = xv; 
        return res; 
    }

    /* Second half of loop divides dd1 into dd, and the matrix returns to its
     * normal arrangement */
    tmpd = dd - dd1;
    if (tmpd < dd1)
    { /* quotient suspected to be 1 */
      tmpu = xu1 + xu;
      tmpv = addll(xv1, xv, &overflow);
      tmp1 = overflow;
    }
    else
    { /* division indicated */
      hiremainder = 0;
      q = 1 + divll(tmpd, dd1, &hiremainder);
      tmpd = hiremainder;
      tmpu = xu1 + q*xu;
      tmpv = addll(xv1, mulll(q, xv, &hiremainder), &overflow);
      tmp1 = overflow | hiremainder;
    }

    tmp0 = addll(tmpu, xu, &overflow);
    if ((tmpd < tmpv) || overflow || tmp1 ||
        (dd1 - tmpd < tmp0)) /* !(Jebelean cond.) */
    { /* The recurrence matrix has not yet been unwarped, so it is
       * the wrong way round;  fix this. */
      *u = xu1; *u1 = xu; *v = xv1; *v1 = xv;
      break;
    }

    res++; /* commit dd, xu1, xv1 */
    dd = tmpd; xu1 = tmpu; xv1 = tmpv;
    if (xv1 > vmax) { 
        *u = xu; 
        *u1 = xu1; 
        *v = xv; 
        *v1 = xv1; 
        return res; 
    }
  } /* end of second loop */

  return res;
}

/* 1 / Mod(x,p). Assume x < p */
ulong
Fl_invsafe(ulong x, ulong p)
{
  int64_t s;
  ulong xv, xv1, g = xgcduu(p, x, 1, &xv, &xv1, &s);
  if (g != 1ULL) return 0ULL;
  xv = xv1 % p; if (s < 0) xv = p - xv;
  return xv;
}

/* 1 / Mod(x,p). Assume x < p */
ulong
Fl_inv(ulong x, ulong p)
{
  ulong xv  = Fl_invsafe(x, p);
  if (!xv && p!=1ULL) 
      pari_err_INV("Fl_inv", mkintmod(utoi(x), utoi(p)));
  return xv;
}
