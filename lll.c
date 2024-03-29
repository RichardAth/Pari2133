/* Copyright (C) 2008  The PARI group.

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

/* default quality ratio for LLL */
static const double LLLDFT = 0.99;

/* assume flag & (LLL_KER|LLL_IM|LLL_ALL). LLL_INPLACE implies LLL_IM */
static GEN
lll_trivial(GEN x, int64_t flag)
{
  if (lg(x) == 1)
  { /* dim x = 0 */
    if (! (flag & LLL_ALL)) return cgetg(1,t_MAT);
    retmkvec2(cgetg(1,t_MAT), cgetg(1,t_MAT));
  }
  /* dim x = 1 */
  if (gequal0(gel(x,1)))
  {
    if (flag & LLL_KER) return matid(1);
    if (flag & (LLL_IM|LLL_INPLACE)) return cgetg(1,t_MAT);
    retmkvec2(matid(1), cgetg(1,t_MAT));
  }
  if (flag & LLL_INPLACE) return gcopy(x);
  if (flag & LLL_KER) return cgetg(1,t_MAT);
  if (flag & LLL_IM)  return matid(1);
  retmkvec2(cgetg(1,t_MAT), (flag & LLL_GRAM)? gcopy(x): matid(1));
}

/* vecslice(x,#x-k,#x) in place. Works for t_MAT, t_VEC/t_COL */
static GEN
vectail_inplace(GEN x, int64_t k)
{
  if (!k) return x;
  x[k] = ((ulong)x[0] & ~LGBITS) | evallg(lg(x) - k);
  return x + k;
}

/* k = dim Kernel */
static GEN
lll_finish(GEN h, int64_t k, int64_t flag)
{
  GEN g;
  if (!(flag & (LLL_IM|LLL_KER|LLL_ALL|LLL_INPLACE))) return h;
  if (flag & (LLL_IM|LLL_INPLACE)) return vectail_inplace(h, k);
  if (flag & LLL_KER) { setlg(h,k+1); return h; }
  g = vecslice(h,1,k); /* done first: vectail_inplace kills h */
  return mkvec2(g, vectail_inplace(h, k));
}

INLINE GEN
mulshift(GEN y, GEN z, int64_t e)
{
  int64_t ly = lgefint(y), lz;
  pari_sp av;
  GEN t;
  if (ly == 2) return gen_0;
  lz = lgefint(z);
  av = avma; (void)new_chunk(ly+lz+nbits2lg(e)); /* HACK */
  t = mulii(z, y);
  set_avma(av); return shifti(t, e);
}

INLINE GEN
submulshift(GEN x, GEN y, GEN z, int64_t e)
{
  int64_t lx = lgefint(x), ly, lz;
  pari_sp av;
  GEN t;
  if (!e) return submulii(x, y, z);
  if (lx == 2) { t = mulshift(y, z, e); togglesign(t); return t; }
  ly = lgefint(y);
  if (ly == 2) return icopy(x);
  lz = lgefint(z);
  av = avma; (void)new_chunk(lx+ly+lz+nbits2lg(e)); /* HACK */
  t = shifti(mulii(z, y), e);
  set_avma(av); return subii(x, t);
}

/********************************************************************/
/**                                                                **/
/**                   FPLLL (adapted from D. Stehle's code)        **/
/**                                                                **/
/********************************************************************/
/* n < 10; gerepileall supporting &NULL arguments. Maybe rename and export ? */
INLINE void
gc_lll(pari_sp av, int n, ...)
{
  int i, j;
  GEN *gptr[10];
  va_list a; va_start(a, n);
  for (i=j=0; i<n; i++)
  {
    GEN *x = va_arg(a,GEN*);
    if (*x) { gptr[j++] = x; *x = (GEN)copy_bin(*x); }
  }
  set_avma(av);
  for (--j; j>=0; j--) *gptr[j] = bin_copy((GENbin*)*gptr[j]);
  va_end(a);
}
/* Babai() and fplll() are a conversion to libpari API and data types
   of the file proved.c in fplll-1.3 by Damien Stehle'.

  Copyright 2005, 2006 Damien Stehle'.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program implements ideas from the paper "Floating-point LLL Revisited",
  by Phong Nguyen and Damien Stehle', in the Proceedings of Eurocrypt'2005,
  Springer-Verlag; and was partly inspired by Shoup's NTL library:
  http://www.shoup.net/ntl/
*/

/* x t_REAL, |x| >= 1/2. Test whether |x| <= 1.5 */
static int
absrsmall(GEN x)
{
  int64_t e = expo(x);
  return (e < 0) || (e == 0 && (ulong)x[2] <= (3ULL << (BITS_IN_LONG-2)));
}

/***********************************************/
/* Babai's Nearest Plane algorithm (iterative) */
/***********************************************/
/* Size-reduces b_kappa using mu_{i,j} and r_{i,j} for j<=i <kappa
Updates B (kappa); computes mu_{kappa,j}, r_{kappa,j} for j<=kappa, and s(kappa)
mu, r, s updated in place (affrr).
*/
static int64_t
Babai(pari_sp av, int64_t kappa, GEN *pG, GEN *pB, GEN *pU, GEN mu, GEN r, GEN s,
      int64_t a, int64_t zeros, int64_t maxG, GEN eta, int64_t prec)
{
  pari_sp av0 = avma;
  GEN G = *pG, B = *pB, U = *pU, tmp, rtmp, ztmp;
  int64_t k, d, n, aa = a > zeros? a: zeros+1;
  GEN maxmu = gen_0, max2mu = gen_0;

  /* N.B: we set d = 0 (resp. n = 0) to avoid updating U (resp. B) */
  d = U? lg(U)-1: 0;
  n = B? nbrows(B): 0;
  if (gc_needed(av,2))
  {
    if(DEBUGMEM>1) pari_warn(warnmem,"Babai[0], a=%ld", aa);
    gc_lll(av,3,&G,&B,&U);
  }
  for (;;) {
    int go_on = 0;
    GEN max3mu;
    int64_t i, j;

    if (gc_needed(av0,2))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Babai[1], a=%ld", aa);
      gc_lll(av,5,&G,&B,&U,&maxmu,&max2mu);
    }
    /* Step2: compute the GSO for stage kappa */
    max3mu = max2mu;
    max2mu = maxmu;
    maxmu = real_0(prec);
    for (j=aa; j<kappa; j++)
    {
      pari_sp btop = avma;
      k = zeros+1;
      if (j > k)
      {
        tmp  = mulrr(gmael(mu,j,k), gmael(r,kappa,k));
        rtmp = subir(gmael(G,kappa,j), tmp);
        for (k++; k<j; k++)
        {
          tmp  = mulrr(gmael(mu,j,k), gmael(r,kappa,k));
          rtmp = subrr(rtmp,tmp);
        }
        affrr(rtmp, gmael(r,kappa,j));
      }
      else
        affir(gmael(G,kappa,j), gmael(r,kappa,j));
      affrr(divrr(gmael(r,kappa,j), gmael(r,j,j)), gmael(mu,kappa,j));
      if (abscmprr(maxmu, gmael(mu,kappa,j))<0)
        maxmu = gmael(mu,kappa,j);
      set_avma(btop);
    }
    maxmu = absr(maxmu); /* copy needed: we 'affrr' to mu later */
    if (typ(max3mu)==t_REAL && abscmprr(max3mu, shiftr(max2mu, 5)) <= 0)
    { /* precision too low */
      *pG = G; *pB = B; *pU = U; return kappa;
    }

    /* Step3--5: compute the X_j's  */
    for (j=kappa-1; j>zeros; j--)
    {
      tmp = gmael(mu,kappa,j);
      if (abscmprr(tmp, eta) <= 0) continue; /* (essentially) size-reduced */

      if (gc_needed(av0,2))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"Babai[2], a=%ld, j=%ld", aa,j);
        gc_lll(av,5,&G,&B,&U,&maxmu,&max2mu);
      }
      go_on = 1;
      /* we consider separately the case |X| = 1 */
      if (absrsmall(tmp))
      {
        if (signe(tmp) > 0) { /* in this case, X = 1 */
          pari_sp btop = avma;
          for (k=zeros+1; k<j; k++)
            affrr(subrr(gmael(mu,kappa,k), gmael(mu,j,k)), gmael(mu,kappa,k));
          set_avma(btop);

          for (i=1; i<=n; i++)
            gmael(B,kappa,i) = subii(gmael(B,kappa,i), gmael(B,j,i));
          for (i=1; i<=d; i++)
            gmael(U,kappa,i) = subii(gmael(U,kappa,i), gmael(U,j,i));
          btop = avma;
          ztmp = subii(gmael(G,j,j), shifti(gmael(G,kappa,j), 1));
          ztmp = addii(gmael(G,kappa,kappa), ztmp);
          gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
          for (i=1; i<=j; i++)
            gmael(G,kappa,i) = subii(gmael(G,kappa,i), gmael(G,j,i));
          for (i=j+1; i<kappa; i++)
            gmael(G,kappa,i) = subii(gmael(G,kappa,i), gmael(G,i,j));
          for (i=kappa+1; i<=maxG; i++)
            gmael(G,i,kappa) = subii(gmael(G,i,kappa), gmael(G,i,j));
        } else { /* otherwise X = -1 */
          pari_sp btop = avma;
          for (k=zeros+1; k<j; k++)
            affrr(addrr(gmael(mu,kappa,k), gmael(mu,j,k)), gmael(mu,kappa,k));
          set_avma(btop);

          for (i=1; i<=n; i++)
            gmael(B,kappa,i) = addii(gmael(B,kappa,i),gmael(B,j,i));
          for (i=1; i<=d; i++)
            gmael(U,kappa,i) = addii(gmael(U,kappa,i),gmael(U,j,i));
          btop = avma;
          ztmp = addii(gmael(G,j,j), shifti(gmael(G,kappa,j), 1));
          ztmp = addii(gmael(G,kappa,kappa), ztmp);
          gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
          for (i=1; i<=j; i++)
            gmael(G,kappa,i) = addii(gmael(G,kappa,i), gmael(G,j,i));
          for (i=j+1; i<kappa; i++)
            gmael(G,kappa,i) = addii(gmael(G,kappa,i), gmael(G,i,j));
          for (i=kappa+1; i<=maxG; i++)
            gmael(G,i,kappa) = addii(gmael(G,i,kappa), gmael(G,i,j));
        }
        continue;
      }
      /* we have |X| >= 2 */
      ztmp = roundr_safe(tmp);
      if (lgefint(ztmp) == 3)
      {
        pari_sp btop = avma;
        ulong xx = ztmp[2]; /* X fits in an ulong */
        if (signe(ztmp) > 0) /* = xx */
        {
          for (k=zeros+1; k<j; k++)
          {
            rtmp = subrr(gmael(mu,kappa,k), mulur(xx, gmael(mu,j,k)));
            affrr(rtmp, gmael(mu,kappa,k));
          }
          set_avma(btop);
          for (i=1; i<=n; i++)
            gmael(B,kappa,i) = submuliu_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
          for (i=1; i<=d; i++)
            gmael(U,kappa,i) = submuliu_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
          btop = avma;
          ztmp = shifti(muliu(gmael(G,kappa,j), xx), 1);
          ztmp = subii(mulii(gmael(G,j,j), sqru(xx)), ztmp);
          ztmp = addii(gmael(G,kappa,kappa), ztmp);
          gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
          for (i=1; i<=j; i++)
            gmael(G,kappa,i) = submuliu_inplace(gmael(G,kappa,i), gmael(G,j,i), xx);
          for (i=j+1; i<kappa; i++)
            gmael(G,kappa,i) = submuliu_inplace(gmael(G,kappa,i), gmael(G,i,j), xx);
          for (i=kappa+1; i<=maxG; i++)
            gmael(G,i,kappa) = submuliu_inplace(gmael(G,i,kappa), gmael(G,i,j), xx);
        }
        else /* = -xx */
        {
          for (k=zeros+1; k<j; k++)
          {
            rtmp = addrr(gmael(mu,kappa,k), mulur(xx, gmael(mu,j,k)));
            affrr(rtmp, gmael(mu,kappa,k));
          }
          set_avma(btop);
          for (i=1; i<=n; i++)
            gmael(B,kappa,i) = addmuliu_inplace(gmael(B,kappa,i), gmael(B,j,i), xx);
          for (i=1; i<=d; i++)
            gmael(U,kappa,i) = addmuliu_inplace(gmael(U,kappa,i), gmael(U,j,i), xx);
          btop = avma;
          ztmp = shifti(muliu(gmael(G,kappa,j), xx), 1);
          ztmp = addii(mulii(gmael(G,j,j), sqru(xx)), ztmp);
          ztmp = addii(gmael(G,kappa,kappa), ztmp);
          gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
          for (i=1; i<=j; i++)
            gmael(G,kappa,i) = addmuliu_inplace(gmael(G,kappa,i), gmael(G,j,i), xx);
          for (i=j+1; i<kappa; i++)
            gmael(G,kappa,i) = addmuliu_inplace(gmael(G,kappa,i), gmael(G,i,j), xx);
          for (i=kappa+1; i<=maxG; i++)
            gmael(G,i,kappa) = addmuliu_inplace(gmael(G,i,kappa), gmael(G,i,j), xx);
        }
      }
      else
      {
        pari_sp btop;
        int64_t e = expi(ztmp) - prec2nbits(prec);
        GEN X = ztmp;

        if (e <= 0) e = 0; else X = shifti(X, -e);
        btop = avma;
        for (k=zeros+1; k<j; k++)
        {
          rtmp = subrr(gmael(mu,kappa,k), mulir(ztmp, gmael(mu,j,k)));
          affrr(rtmp, gmael(mu,kappa,k));
        }
        set_avma(btop);
        for (i=1; i<=n; i++)
          gmael(B,kappa,i) = submulshift(gmael(B,kappa,i), gmael(B,j,i), X, e);
        for (i=1; i<=d; i++)
          gmael(U,kappa,i) = submulshift(gmael(U,kappa,i), gmael(U,j,i), X, e);
        btop = avma;
        ztmp = shifti(mulii(gmael(G,kappa,j), X), e+1);
        ztmp = subii(shifti(mulii(gmael(G,j,j), sqri(X)), 2*e), ztmp);
        ztmp = addii(gmael(G,kappa,kappa), ztmp);
        gmael(G,kappa,kappa) = gerepileuptoint(btop, ztmp);
        for (i=1; i<=j; i++)
          gmael(G,kappa,i) = submulshift(gmael(G,kappa,i), gmael(G,j,i), X, e);
        for (   ; i<kappa; i++)
          gmael(G,kappa,i) = submulshift(gmael(G,kappa,i), gmael(G,i,j), X, e);
        for (i=kappa+1; i<=maxG; i++)
          gmael(G,i,kappa) = submulshift(gmael(G,i,kappa), gmael(G,i,j), X, e);
      }
    }
    if (!go_on) break; /* Anything happened? */
    aa = zeros+1;
  }

  affir(gmael(G,kappa,kappa), gel(s,zeros+1));
  /* the last s[kappa-1]=r[kappa][kappa] is computed only if kappa increases */
  av = avma;
  for (k=zeros+1; k<=kappa-2; k++)
  {
    tmp = subrr(gel(s,k), mulrr(gmael(mu,kappa,k), gmael(r,kappa,k)));
    affrr(tmp, gel(s,k+1));
  }
  *pG = G; *pB = B; *pU = U; return gc_long(av, 0);
}

static void
rotate(GEN A, int64_t k2, int64_t k)
{
  int64_t i;
  GEN B = gel(A,k2);
  for (i = k2; i > k; i--) gel(A,i) = gel(A,i-1);
  gel(A,k) = B;
}

/* ****************** */
/* The LLL Algorithm  */
/* ****************** */
/* G integral Gram matrix, LLL-reduces (G,B,U) in place [apply base change
 * transforms to B and U]. If (keepfirst), never swap with first vector.
 * If G = NULL, we compute the Gram matrix incrementally.
 * Return -1 on failure, else zeros = dim Kernel (>= 0) */
static int64_t
fplll(GEN *pG, GEN *pB, GEN *pU, GEN *pr, double DELTA, double ETA,
      int64_t keepfirst, int64_t prec)
{
  pari_sp av, av2;
  GEN mu, r, s, tmp, SPtmp, alpha, G = *pG, B = *pB, U = *pU;
  GEN delta = dbltor(DELTA), eta = dbltor(ETA);
  int64_t d, maxG, kappa, kappa2, i, j, zeros, kappamax, incgram = !G, cnt = 0;

  if (incgram)
  { /* incremental Gram matrix */
    maxG = 2; d = lg(B)-1;
    G = zeromatcopy(d, d);
  }
  else
    maxG = d = lg(G)-1;

  mu = cgetg(d+1, t_MAT);
  r  = cgetg(d+1, t_MAT);
  s  = cgetg(d+1, t_VEC);
  for (j = 1; j <= d; j++)
  {
    GEN M = cgetg(d+1, t_COL), R = cgetg(d+1, t_COL);
    gel(mu,j)= M;
    gel(r,j) = R;
    gel(s,j) = cgetr(prec);
    for (i = 1; i <= d; i++)
    {
      gel(R,i) = cgetr(prec);
      gel(M,i) = cgetr(prec);
    }
  }
  SPtmp = cgetg(d+1, t_VEC);
  alpha = cgetg(d+1, t_VECSMALL);
  av = avma;

  /* Step2: Initializing the main loop */
  kappamax = 1;
  i = 1;
  do {
    if (incgram) gmael(G,i,i) = ZV_dotsquare(gel(B,i));
    affir(gmael(G,i,i), gmael(r,i,i));
  } while (!signe(gmael(G,i,i)) && ++i <= d);
  zeros = i-1; /* all basis vectors b_i with i <= zeros are zero vectors */
  kappa = i;
  if (zeros < d) affir(gmael(G,zeros+1,zeros+1), gmael(r,zeros+1,zeros+1));
  for (i=zeros+1; i<=d; i++) alpha[i]=1;

  while (++kappa <= d)
  {
    if (kappa > kappamax)
    {
      if (DEBUGLEVEL>=4) err_printf("K%ld ",kappa);
      kappamax = kappa;
      if (incgram)
      {
        for (i=zeros+1; i<=kappa; i++)
          gmael(G,kappa,i) = ZV_dotproduct(gel(B,kappa), gel(B,i));
        maxG = kappamax;
      }
    }
    /* Step3: Call to the Babai algorithm, mu,r,s updated in place */
    if (Babai(av, kappa, &G,&B,&U, mu,r,s, alpha[kappa], zeros, maxG, eta,prec))
    {
      if (incgram) G = NULL;
      *pG = G; *pB = B; *pU = U; return -1;
    }
    av2 = avma;
    if ((keepfirst && kappa == 2) ||
        cmprr(mulrr(gmael(r,kappa-1,kappa-1), delta), gel(s,kappa-1)) <= 0)
    { /* Step4: Success of Lovasz's condition */
      alpha[kappa] = kappa;
      tmp = mulrr(gmael(mu,kappa,kappa-1), gmael(r,kappa,kappa-1));
      affrr(subrr(gel(s,kappa-1), tmp), gmael(r,kappa,kappa));
      set_avma(av2); continue;
    }
    /* Step5: Find the right insertion index kappa, kappa2 = initial kappa */
    if (DEBUGLEVEL>=4 && kappa==kappamax && signe(gel(s,kappa-1)))
      if (++cnt > 20) { cnt = 0; err_printf("(%ld) ", expo(gel(s,1))); }
    kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros+2 + (keepfirst ? 1: 0)) break;
      tmp = mulrr(gmael(r,kappa-1,kappa-1), delta);
    } while (cmprr(gel(s,kappa-1), tmp) <= 0);
    set_avma(av2);

    for (i=kappa; i<kappa2; i++)
      if (kappa <= alpha[i]) alpha[i] = kappa;
    for (i=kappa2; i>kappa; i--) alpha[i] = alpha[i-1];
    for (i=kappa2+1; i<=kappamax; i++)
      if (kappa < alpha[i]) alpha[i] = kappa;
    alpha[kappa] = kappa;

    /* Step6: Update the mu's and r's */
    rotate(mu, kappa2, kappa);
    rotate(r, kappa2, kappa);
    affrr(gel(s,kappa), gmael(r,kappa,kappa));

    /* Step7: Update G, B, U */
    if (U) rotate(U, kappa2, kappa);
    if (B) rotate(B, kappa2, kappa);
    for (i=1; i<=kappa2; i++) gel(SPtmp,i) = gmael(G,kappa2,i);
    for (   ; i<=maxG;   i++) gel(SPtmp,i) = gmael(G,i,kappa2);
    for (i=kappa2; i>kappa; i--)
    {
      for (j=1; j<kappa; j++) gmael(G,i,j) = gmael(G,i-1,j);
      gmael(G,i,kappa) = gel(SPtmp,i-1);
      for (j=kappa+1; j<=i; j++) gmael(G,i,j) = gmael(G,i-1,j-1);
      for (j=kappa2+1; j<=maxG; j++) gmael(G,j,i) = gmael(G,j,i-1);
    }
    for (i=1; i<kappa; i++) gmael(G,kappa,i) = gel(SPtmp,i);
    gmael(G,kappa,kappa) = gel(SPtmp,kappa2);
    for (i=kappa2+1; i<=maxG; i++) gmael(G,i,kappa) = gel(SPtmp,i);

    /* Step8: Prepare the next loop iteration */
    if (kappa == zeros+1 && !signe(gmael(G,kappa,kappa)))
    {
      zeros++; kappa++;
      affir(gmael(G,kappa,kappa), gmael(r,kappa,kappa));
    }
  }
  if (pr) *pr = RgM_diagonal_shallow(r);
  *pG = G; *pB = B; *pU = U; return zeros; /* success */
}

/* Assume x a ZM, if pN != NULL, set it to Gram-Schmidt (squared) norms
 * The following modes are supported:
 * - flag & LLL_INPLACE: x a lattice basis, return x*U
 * - flag & LLL_GRAM: x a Gram matrix / else x a lattice basis; return
 *     LLL base change matrix U [LLL_IM]
 *     kernel basis [LLL_KER, nonreduced]
 *     both [LLL_ALL] */
GEN
ZM_lll_norms(GEN x, double DELTA, int64_t flag, GEN *pN)
{
  pari_sp av = avma;
  const double ETA = 0.51;
  int64_t p, zeros, n = lg(x)-1;
  GEN G, B, U;
  pari_timer T;

  if (n <= 1) return lll_trivial(x, flag);
  x = RgM_shallowcopy(x);
  if (flag & LLL_GRAM)
  { G = x; B = NULL; U = matid(n); }
  else
  { G = NULL; B = x; U = (flag & LLL_INPLACE)? NULL: matid(n); }
  for (p = DEFAULTPREC;; p += EXTRAPREC64)
  {
    if (DEBUGLEVEL>=4)
    {
      err_printf("Entering L^2: LLL-parameters (%.3f,%.3f), prec = %d\n",
                 DELTA,ETA, p);
      timer_start(&T);
    }
    zeros = fplll(&G, &B, &U, pN, DELTA, ETA, flag & LLL_KEEP_FIRST, p);
    if (zeros >= 0) break;
    gc_lll(av, 3, &G, &B, &U);
  }
  if (DEBUGLEVEL>=4) timer_printf(&T,"LLL");
  return lll_finish(U? U: B, zeros, flag);
}

/********************************************************************/
/**                                                                **/
/**                        LLL OVER K[X]                           **/
/**                                                                **/
/********************************************************************/
static int
pslg(GEN x)
{
  int64_t tx;
  if (gequal0(x)) return 2;
  tx = typ(x); return is_scalar_t(tx)? 3: lg(x);
}

static int
REDgen(int64_t k, int64_t l, GEN h, GEN L, GEN B)
{
  GEN q, u = gcoeff(L,k,l);
  int64_t i;

  if (pslg(u) < pslg(B)) return 0;

  q = gneg(gdeuc(u,B));
  gel(h,k) = gadd(gel(h,k), gmul(q,gel(h,l)));
  for (i=1; i<l; i++) gcoeff(L,k,i) = gadd(gcoeff(L,k,i), gmul(q,gcoeff(L,l,i)));
  gcoeff(L,k,l) = gadd(gcoeff(L,k,l), gmul(q,B)); return 1;
}

static int
do_SWAPgen(GEN h, GEN L, GEN B, int64_t k, GEN fl, int *flc)
{
  GEN p1, la, la2, Bk;
  int64_t ps1, ps2, i, j, lx;

  if (!fl[k-1]) return 0;

  la = gcoeff(L,k,k-1); la2 = gsqr(la);
  Bk = gel(B,k);
  if (fl[k])
  {
    GEN q = gadd(la2, gmul(gel(B,k-1),gel(B,k+1)));
    ps1 = pslg(gsqr(Bk));
    ps2 = pslg(q);
    if (ps1 <= ps2 && (ps1 < ps2 || !*flc)) return 0;
    *flc = (ps1 != ps2);
    gel(B,k) = gdiv(q, Bk);
  }

  swap(gel(h,k-1), gel(h,k)); lx = lg(L);
  for (j=1; j<k-1; j++) swap(gcoeff(L,k-1,j), gcoeff(L,k,j));
  if (fl[k])
  {
    for (i=k+1; i<lx; i++)
    {
      GEN t = gcoeff(L,i,k);
      p1 = gsub(gmul(gel(B,k+1),gcoeff(L,i,k-1)), gmul(la,t));
      gcoeff(L,i,k) = gdiv(p1, Bk);
      p1 = gadd(gmul(la,gcoeff(L,i,k-1)), gmul(gel(B,k-1),t));
      gcoeff(L,i,k-1) = gdiv(p1, Bk);
    }
  }
  else if (!gequal0(la))
  {
    p1 = gdiv(la2, Bk);
    gel(B,k+1) = gel(B,k) = p1;
    for (i=k+2; i<=lx; i++) gel(B,i) = gdiv(gmul(p1,gel(B,i)),Bk);
    for (i=k+1; i<lx; i++)
      gcoeff(L,i,k-1) = gdiv(gmul(la,gcoeff(L,i,k-1)), Bk);
    for (j=k+1; j<lx-1; j++)
      for (i=j+1; i<lx; i++)
        gcoeff(L,i,j) = gdiv(gmul(p1,gcoeff(L,i,j)), Bk);
  }
  else
  {
    gcoeff(L,k,k-1) = gen_0;
    for (i=k+1; i<lx; i++)
    {
      gcoeff(L,i,k) = gcoeff(L,i,k-1);
      gcoeff(L,i,k-1) = gen_0;
    }
    gel(B,k) = gel(B,k-1); fl[k] = 1; fl[k-1] = 0;
  }
  return 1;
}

static void
incrementalGSgen(GEN x, GEN L, GEN B, int64_t k, GEN fl)
{
  GEN u = NULL; /* gcc -Wall */
  int64_t i, j;
  for (j = 1; j <= k; j++)
    if (j==k || fl[j])
    {
      u = gcoeff(x,k,j);
      if (!is_extscalar_t(typ(u))) pari_err_TYPE("incrementalGSgen",u);
      for (i=1; i<j; i++)
        if (fl[i])
        {
          u = gsub(gmul(gel(B,i+1),u), gmul(gcoeff(L,k,i),gcoeff(L,j,i)));
          u = gdiv(u, gel(B,i));
        }
      gcoeff(L,k,j) = u;
    }
  if (gequal0(u)) gel(B,k+1) = gel(B,k);
  else
  {
    gel(B,k+1) = gcoeff(L,k,k); gcoeff(L,k,k) = gen_1; fl[k] = 1;
  }
}

static GEN
lllgramallgen(GEN x, int64_t flag)
{
  int64_t lx = lg(x), i, j, k, l, n;
  pari_sp av;
  GEN B, L, h, fl;
  int flc;

  n = lx-1; if (n<=1) return lll_trivial(x,flag);
  if (lgcols(x) != lx) pari_err_DIM("lllgramallgen");

  fl = cgetg(lx, t_VECSMALL);

  av = avma;
  B = scalarcol_shallow(gen_1, lx);
  L = cgetg(lx,t_MAT);
  for (j=1; j<lx; j++) { gel(L,j) = zerocol(n); fl[j] = 0; }

  h = matid(n);
  for (i=1; i<lx; i++)
    incrementalGSgen(x, L, B, i, fl);
  flc = 0;
  for(k=2;;)
  {
    if (REDgen(k, k-1, h, L, gel(B,k))) flc = 1;
    if (do_SWAPgen(h, L, B, k, fl, &flc)) { if (k > 2) k--; }
    else
    {
      for (l=k-2; l>=1; l--)
        if (REDgen(k, l, h, L, gel(B,l+1))) flc = 1;
      if (++k > n) break;
    }
    if (gc_needed(av,1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lllgramallgen");
      gerepileall(av,3,&B,&L,&h);
    }
  }
  k=1; while (k<lx && !fl[k]) k++;
  return lll_finish(h,k-1,flag);
}

static int
RgM_square(GEN x) { int64_t l = lg(x); return l == 1 || l == lgcols(x); }
static GEN
lllallgen(GEN x, int64_t flag)
{
  pari_sp av = avma;
  if (!(flag & LLL_GRAM)) x = gram_matrix(x);
  else if (!RgM_square(x)) pari_err_DIM("qflllgram");
  return gerepilecopy(av, lllgramallgen(x, flag));
}
GEN
lllgen(GEN x) { return lllallgen(x, LLL_IM); }
GEN
lllkerimgen(GEN x) { return lllallgen(x, LLL_ALL); }
GEN
lllgramgen(GEN x)  { return lllallgen(x, LLL_IM|LLL_GRAM); }
GEN
lllgramkerimgen(GEN x)  { return lllallgen(x, LLL_ALL|LLL_GRAM); }

static GEN
lllall(GEN x, int64_t flag)
{ pari_sp av = avma; return gerepilecopy(av, ZM_lll(x, LLLDFT, flag)); }
GEN
lllint(GEN x) { return lllall(x, LLL_IM); }
GEN
lllkerim(GEN x) { return lllall(x, LLL_ALL); }
GEN
lllgramint(GEN x)
{ if (!RgM_square(x)) pari_err_DIM("qflllgram");
  return lllall(x, LLL_IM | LLL_GRAM); }
GEN
lllgramkerim(GEN x)
{ if (!RgM_square(x)) pari_err_DIM("qflllgram");
  return lllall(x, LLL_ALL | LLL_GRAM); }

GEN
lllfp(GEN x, double D, int64_t flag)
{
  int64_t n = lg(x)-1;
  pari_sp av = avma;
  GEN h;
  if (n <= 1) return lll_trivial(x,flag);
  if ((flag & LLL_GRAM) && !RgM_square(x)) pari_err_DIM("qflllgram");
  h = ZM_lll(RgM_rescale_to_int(x), D, flag);
  return gerepilecopy(av, h);
}

GEN
lllgram(GEN x) { return lllfp(x,LLLDFT,LLL_GRAM|LLL_IM); }
GEN
lll(GEN x) { return lllfp(x,LLLDFT,LLL_IM); }

GEN
qflll0(GEN x, int64_t flag)
{
  if (typ(x) != t_MAT) pari_err_TYPE("qflll",x);
  switch(flag)
  {
    case 0: return lll(x);
    case 1: RgM_check_ZM(x,"qflll"); return lllint(x);
    case 2: RgM_check_ZM(x,"qflll"); return lllintpartial(x);
    case 3: RgM_check_ZM(x,"qflll"); return lllall(x, LLL_INPLACE);
    case 4: RgM_check_ZM(x,"qflll"); return lllkerim(x);
    case 5: return lllkerimgen(x);
    case 8: return lllgen(x);
    default: pari_err_FLAG("qflll");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

GEN
qflllgram0(GEN x, int64_t flag)
{
  if (typ(x) != t_MAT) pari_err_TYPE("qflllgram",x);
  switch(flag)
  {
    case 0: return lllgram(x);
    case 1: RgM_check_ZM(x,"qflllgram"); return lllgramint(x);
    case 4: RgM_check_ZM(x,"qflllgram"); return lllgramkerim(x);
    case 5: return lllgramkerimgen(x);
    case 8: return lllgramgen(x);
    default: pari_err_FLAG("qflllgram");
  }
  return NULL; /* LCOV_EXCL_LINE */
}

/********************************************************************/
/**                                                                **/
/**                   INTEGRAL KERNEL (LLL REDUCED)                **/
/**                                                                **/
/********************************************************************/
static GEN
kerint0(GEN M)
{
  /* return ZM_lll(M, LLLDFT, LLL_KER); */
  GEN U, H = ZM_hnflll(M,&U,1);
  int64_t d = lg(M)-lg(H);
  if (!d) return cgetg(1, t_MAT);
  return ZM_lll(vecslice(U,1,d), LLLDFT, LLL_INPLACE);
}
GEN
kerint(GEN M)
{
  pari_sp av = avma;
  return gerepilecopy(av, kerint0(M));
}
/* OBSOLETE: use kerint */
GEN
matkerint0(GEN M, int64_t flag)
{
  pari_sp av = avma;
  if (typ(M) != t_MAT) pari_err_TYPE("matkerint",M);
  M = Q_primpart(M);
  RgM_check_ZM(M, "kerint");
  switch(flag)
  {
    case 0:
    case 1: return gerepilecopy(av, kerint0(M));
    default: pari_err_FLAG("matkerint");
  }
  return NULL; /* LCOV_EXCL_LINE */
}
