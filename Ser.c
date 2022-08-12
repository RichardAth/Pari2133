/* Copyright (C) 2000, 2012  The PARI group.

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
/*******************************************************************/
/*                                                                 */
/*                     Conversion --> t_SER                        */
/*                                                                 */
/*******************************************************************/
static GEN
RgX_to_ser_i(GEN x, int64_t l, int64_t v, int copy)
{
  int64_t i = 2, lx = lg(x), vx = varn(x);
  GEN y;
  if (lx == 2) return zeroser(vx, minss(l - 2, v));
  if (l <= 2)
  {
    if (l == 2 && v != LONG_MAX) return zeroser(vx, v);
    pari_err_BUG("RgX_to_ser (l < 2)");
  }
  y = cgetg(l,t_SER);
  if (v == LONG_MAX) { v = 1; lx = 3; } /* e.g. Mod(0,3) * x^0 */
  else if (v)
  {
    int64_t w = v;
    while (isrationalzero(gel(x,2))) { x++; w--; }
    lx -= v;
    if (w)
    { /* keep type information, e.g. Mod(0,3) + x */
      GEN z = gel(x,2); /* = 0 */
      x += w; gel(y,2) = gadd(gel(x,2), z); i++;
    }
  }
  y[1] = evalvarn(vx) | evalvalp(v); /* must come here because of LONG_MAX */
  if (lx > l) lx = l;
  /* 3 <= lx <= l */
  if (copy)
    for (; i <lx; i++) gel(y,i) = gcopy(gel(x,i));
  else
    for (; i <lx; i++) gel(y,i) = gel(x,i);
  for (     ; i < l; i++) gel(y,i) = gen_0;
  return normalize(y);
}
/* enlarge/truncate t_POL x to a t_SER with lg l */
GEN
RgX_to_ser(GEN x, int64_t l) { return RgX_to_ser_i(x, l, RgX_val(x), 0); }
GEN
RgX_to_ser_inexact(GEN x, int64_t l)
{
  int64_t i, lx = lg(x);
  int first = 1;
  for (i = 2; i < lx && gequal0(gel(x,i)); i++) /* ~ RgX_valrem + normalize */
    if (first && !isexactzero(gel(x,i)))
    {
      pari_warn(warner,"normalizing a series with 0 leading term");
      first = 0;
    }
  return RgX_to_ser_i(x, l, i - 2, 0);
}
GEN
rfrac_to_ser(GEN x, int64_t l)
{
  GEN d = gel(x,2);
  if (l == 2)
  {
    int64_t v = varn(d);
    return zeroser(varn(d), gvaluation(x, pol_x(v)));
  }
  return gdiv(gel(x,1), RgX_to_ser(d, l));
}

static GEN
RgV_to_ser_i(GEN x, int64_t v, int64_t l, int copy)
{
  int64_t j, lx = minss(lg(x), l-1);
  GEN y;
  if (lx == 1) return zeroser(v, l-2);
  y = cgetg(l, t_SER); y[1] = evalvarn(v)|evalvalp(0);
  x--;
  if (copy)
    for (j = 2; j <= lx; j++) gel(y,j) = gcopy(gel(x,j));
  else
    for (j = 2; j <= lx; j++) gel(y,j) = gel(x,j);
  for (     ; j < l;   j++) gel(y,j) = gen_0;
  return normalize(y);
}
GEN
RgV_to_ser(GEN x, int64_t v, int64_t l) { return RgV_to_ser_i(x, v, l, 0); }

/* x a t_SER, prec >= 0 */
GEN
sertoser(GEN x, int64_t prec)
{
  int64_t i, lx = lg(x), l;
  GEN y;
  if (lx == 2) return zeroser(varn(x), prec);
  l = prec+2; lx = minss(lx, l);
  y = cgetg(l,t_SER); y[1] = x[1];
  for (i = 2; i < lx; i++) gel(y,i) = gel(x,i);
  for (     ; i < l;  i++) gel(y,i) = gen_0;
  return y;
}

/* R(1/x) = x^v * n/d, val(n) = val(d) = 0 */
int64_t
rfracrecip(GEN *pn, GEN *pd)
{
  int64_t v = degpol(*pd);
  if (typ(*pn) == t_POL && varn(*pn) == varn(*pd))
  {
    v -= degpol(*pn);
    (void)RgX_valrem(*pn, pn); *pn = RgX_recip(*pn);
  }
  (void)RgX_valrem(*pd, pd); *pd = RgX_recip(*pd);
  return v;
}

/* R(1/x) + O(x^N) */
GEN
rfracrecip_to_ser_absolute(GEN R, int64_t N)
{
  GEN n = gel(R,1), d = gel(R,2);
  int64_t v = rfracrecip(&n, &d); /* R(1/x) = x^v * n/d, val(n) = val(d) = 0 */
  if (N <= v) return zeroser(varn(d), N);
  R = gdiv(n, RgX_to_ser(d, N-v+2));
  setvalp(R, v); return R;
}

/* assume prec >= 0 */
GEN
scalarser(GEN x, int64_t v, int64_t prec)
{
  int64_t i, l;
  GEN y;

  if (gequal0(x))
  {
    if (isrationalzero(x)) return zeroser(v, prec);
    if (!isexactzero(x)) prec--;
    y = cgetg(3, t_SER);
    y[1] = evalsigne(0) | _evalvalp(prec) | evalvarn(v);
    gel(y,2) = gcopy(x); return y;
  }
  l = prec + 2; y = cgetg(l, t_SER);
  y[1] = evalsigne(1) | _evalvalp(0) | evalvarn(v);
  gel(y,2) = gcopy(x); for (i=3; i<l; i++) gel(y,i) = gen_0;
  return y;
}

/* assume x a t_[SER|POL], apply gtoser to all coeffs */
static GEN
coefstoser(GEN x, int64_t v, int64_t prec)
{
  int64_t i, lx;
  GEN y = cgetg_copy(x, &lx); y[1] = x[1];
  for (i=2; i<lx; i++) gel(y,i) = gtoser(gel(x,i), v, prec);
  return y;
}

static void
err_ser_priority(GEN x, int64_t v) { pari_err_PRIORITY("Ser", x, "<", v); }
/* x a t_POL */
static GEN
poltoser(GEN x, int64_t v, int64_t prec)
{
  int64_t s = varncmp(varn(x), v);
  if (s < 0) err_ser_priority(x,v);
  if (s > 0) return scalarser(x, v, prec);
  return RgX_to_ser_i(x, prec+2, RgX_val(x), 1);
}
/* x a t_RFRAC */
static GEN
rfractoser(GEN x, int64_t v, int64_t prec)
{
  int64_t s = varncmp(varn(gel(x,2)), v);
  pari_sp av;
  if (s < 0) err_ser_priority(x,v);
  if (s > 0) return scalarser(x, v, prec);
  av = avma; return gerepileupto(av, rfrac_to_ser(x, prec+2));
}
GEN
toser_i(GEN x)
{
  switch(typ(x))
  {
    case t_SER: return x;
    case t_POL: return RgX_to_ser(x, precdl+2);
    case t_RFRAC: return rfrac_to_ser(x, precdl+2);
  }
  return NULL;
}

/* conversion: prec ignored if t_VEC or t_SER in variable v */
GEN
gtoser(GEN x, int64_t v, int64_t prec)
{
  int64_t tx = typ(x);

  if (v < 0) v = 0;
  if (prec < 0) pari_err_DOMAIN("Ser", "precision", "<", gen_0, stoi(prec));
  if (tx == t_SER)
  {
    int64_t s = varncmp(varn(x), v);
    if      (s < 0) return coefstoser(x, v, prec);
    else if (s > 0) return scalarser(x, v, prec);
    return gcopy(x);
  }
  if (is_scalar_t(tx)) return scalarser(x,v,prec);
  switch(tx)
  {
    case t_POL: return poltoser(x, v, prec);
    case t_RFRAC: return rfractoser(x, v, prec);
    case t_QFR:
    case t_QFI: return RgV_to_ser_i(x, v, 4+1, 1);
    case t_VECSMALL: x = zv_to_ZV(x);/*fall through*/
    case t_VEC: case t_COL:
      if (varncmp(gvar(x), v) <= 0) pari_err_PRIORITY("Ser", x, "<=", v);
      return RgV_to_ser_i(x, v, lg(x)+1, 1);
  }
  pari_err_TYPE("Ser",x);
  return NULL; /* LCOV_EXCL_LINE */
}
/* impose prec */
GEN
gtoser_prec(GEN x, int64_t v, int64_t prec)
{
  pari_sp av = avma;
  if (v < 0) v = 0;
  if (prec < 0) pari_err_DOMAIN("Ser", "precision", "<", gen_0, stoi(prec));
  switch(typ(x))
  {
    case t_SER: if (varn(x) != v) break;
                return gerepilecopy(av, sertoser(x, prec));
    case t_QFR:
    case t_QFI:
      x = RgV_to_ser_i(mkvec3(gel(x,1),gel(x,2),gel(x,3)), v, prec+2, 1);
      return gerepileupto(av, x);
    case t_VECSMALL: x = zv_to_ZV(x);/*fall through*/
    case t_VEC: case t_COL:
      if (varncmp(gvar(x), v) <= 0) pari_err_PRIORITY("Ser", x, "<=", v);
      return RgV_to_ser_i(x, v, prec+2, 1);
  }
  return gtoser(x, v, prec);
}
GEN
Ser0(GEN x, int64_t v, GEN d, int64_t prec)
{
  if (!d) return gtoser(x, v, prec);
  if (typ(d) != t_INT)
  {
    d = gceil(d);
    if (typ(d) != t_INT) pari_err_TYPE("Ser [precision]",d);
  }
  return gtoser_prec(x, v, itos(d));
}
