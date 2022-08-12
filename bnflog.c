/* Copyright (C) 2016  The PARI group.

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

/*******************************************************************/
/*                  LOGARITHMIC CLASS GROUP                        */
/*******************************************************************/
/* min(v, v(Log_p Norm_{F_\p/Q_p}(x))) */
static int64_t
vlognorm(GEN nf, GEN T, GEN x, GEN p, int64_t v)
{
  GEN a = nf_to_scalar_or_alg(nf, x);
  GEN N = RgXQ_norm(a, T);
  if (typ(N) != t_PADIC) N = cvtop(N, p, v);
  return minss(v, valp( Qp_log(N) ));
}
/* K number field, pr a maximal ideal, let K_pr be the attached local
 * field, K_pr = Q_p[X] / (T), T irreducible. Return \tilde{e}(K_pr/Q_p) */
static int64_t
etilde(GEN nf, GEN pr, GEN T)
{
  GEN gp = pr_get_p(pr);
  ulong e = pr_get_e(pr);
  int64_t v, voo, vmin, p, k;

  if (!u_pval(e, gp))
  {
    v = u_pval(pr_get_f(pr), gp);
    return itou( mului(e, powiu(gp, v)) );
  }
  nf = checknf(nf);
  p = itou(gp);
  k = e / (p-1) + 1;
  /* log Norm_{F_P/Q_p} (1 + P^k) = Tr(P^k) = p^[(k + v(Diff))/ e] Z_p */
  voo = (k + idealval(nf, nf_get_diff(nf), pr)) / e;
  vmin = vlognorm(nf, T, pr_get_gen(pr), gp, voo);
  if (k > 1)
  {
    GEN U = idealprincipalunits(nf, pr, k);
    GEN gen = abgrp_get_gen(U), cyc = abgrp_get_cyc(U);
    int64_t i, l = lg(cyc);
    for (i = 1; i < l; i++)
    {
      if (voo - Z_lval(gel(cyc,i), p) >= vmin) break;
      vmin = vlognorm(nf, T, gel(gen,i), gp, vmin);
    }
  }
  v = u_lval(degpol(T), p) + (p == 2UL? 2 : 1) - vmin;
  (void)u_lvalrem(e, p, &e);
  return e * upowuu(p,v);
}
static int64_t
ftilde_from_e(GEN pr, int64_t e) { return pr_get_e(pr) * pr_get_f(pr) / e; }
static int64_t
ftilde(GEN K, GEN pr, GEN T) { return ftilde_from_e(pr, etilde(K,pr, T)); }

static int64_t
get_ZpX_index(GEN K, GEN pr, GEN T)
{
  GEN p, pi;
  int64_t j, l = lg(T);
  if (l == 2) return 1;
  p = pr_get_p(pr); pi = nf_to_scalar_or_alg(K, pr_get_gen(pr));
  for (j = 1; j < l; j++)
  {
    GEN t = gel(T,j);
    if (t && gvaluation(RgXQ_norm(pi, t), p)) return j;
  }
  return 0;
}

/* Given a number field K and a prime p, return
 * S = places of K above p [primedec]
 * R = corresponding p-adic factors of K.pol (mod p^k), in the same order */
static GEN
padicfact(GEN K, GEN S, int64_t k)
{
    GEN R, p = pr_get_p(gel(S,1));
  GEN T = gel(factorpadic(nf_get_pol(K), p, k), 1);
  int64_t l, i;
  S = idealprimedec(K, p);
  R = cgetg_copy(S, &l);
  for (i = 1; i < l; i++)
  {
    int64_t j = get_ZpX_index(K, gel(S,i), T);
    gel(R,i) = gel(T,j);
    gel(T,j) = NULL;
  }
  return R;
}

/* K a bnf, compute Cl'(K) = ell-Sylow of Cl(K) / (places above ell).
 * Return [D, u, R0, U0, ordS]
 * - D: cyclic factors for Cl'(K)
 * - u: generators of cyclic factors (all coprime to ell)
 * - R0: subgroup isprincipal(<S>) (divides K.cyc)
 * - U0: generators of R0 are of the form S . U0
 * - ordS[i] = order of S[i] in CL(K)  */
static GEN
CL_prime(GEN K, GEN ell, GEN Sell)
{
  GEN g, ordS, R0, U0, U, D, u, cyc = bnf_get_cyc(K);
  int64_t i, l, lD, lS = lg(Sell);

  g = leafcopy(bnf_get_gen(K));
  l = lg(g);
  for (i = 1; i < l; i++)
  {
    GEN A = gel(g,i), a = gcoeff(A,1,1);
    int64_t v = Z_pvalrem(a, ell, &a);
    if (v) gel(g,i) = hnfmodid(A, a); /* make coprime to ell */
  }
  R0 = cgetg(lS, t_MAT);
  ordS = cgetg(lS, t_VEC);
  for (i = 1; i < lS; i++)
  {
    gel(R0,i) = isprincipal(K, gel(Sell,i));
    gel(ordS,i) = charorder(cyc, gel(R0,i)); /* order of Sell[i] */
  }
  R0 = shallowconcat(R0, diagonal_shallow(cyc));
  /* R0 = subgroup generated by S in Cl(K) [ divides diagonal(K.cyc) ]*/
  R0 = ZM_hnfall(R0, &U0, 2); /* [S | cyc] * U0 = R0 in HNF */
  D = ZM_snfall(R0, &U,NULL);
  D = RgM_diagonal_shallow(D);
  lD = lg(D);
  u = ZM_inv(U, NULL); settyp(u, t_VEC);
  for (i = 1; i < lD; i++) gel(u,i) = idealfactorback(K,g,gel(u,i),1);
  setlg(U0, l);
  U0 = rowslice(U0,1,lS-1); /* restrict to 'S' part */
  return mkvec5(D, u, R0, U0, ordS);
}

static GEN
ell1(GEN ell) { return equaliu(ell,2)? utoipos(5): addiu(ell,1); }

/* log N_{F_P/Q_p}(x) */
static GEN
vtilde_i(GEN K, GEN x, GEN T, GEN ell, int64_t prec)
{
  GEN N, cx;
  if (typ(x) != t_POL) x = nf_to_scalar_or_alg(K, x);
  if (typ(x) != t_POL) { cx = x; N = NULL; }
  else
  {
    x = Q_primitive_part(x,&cx);
    N = resultant(RgX_rem(x,T), T);
    N = cvtop(N,ell,prec);
  }
  if (cx)
  {
    (void)Q_pvalrem(cx, ell, &cx);
    if (!isint1(cx))
    {
      cx = gpowgs(cvtop(cx,ell,prec), degpol(T));
      N = N? gmul(N, cx): cx;
    }
  }
  return N? Qp_log(N): gen_0;
}
static GEN
vecvtilde_i(GEN K, GEN x, GEN T, GEN ell, int64_t prec)
{ pari_APPLY_same(vtilde_i(K, gel(x,i), T, ell, prec)); }
static GEN
vtilde(GEN K, GEN x, GEN T, GEN deg, GEN ell, int64_t prec)
{
  pari_sp av;
  GEN v, G, E;
  if (typ(x) != t_MAT) return gdiv(vtilde_i(K,x,T,ell,prec), deg);
  G = gel(x,1);
  E = gel(x,2); av = avma; v = vecvtilde_i(K,G,T,ell,prec);
  return gerepileupto(av, gdiv(RgV_dotproduct(E, v), deg));
}

/* v[i] = deg S[i] mod p^prec */
static GEN
get_vdegS(GEN Ftilde, GEN ell, int64_t prec)
{
  int64_t i, l = lg(Ftilde);
  GEN v = cgetg(l, t_VEC), degell = Qp_log( cvtop(ell1(ell), ell, prec) );
  for (i = 1; i < l; i++) gel(v,i) = gmulsg(Ftilde[i], degell);
  return v;
}
/* K a bnf. Compute kernel \tilde{Cl}_K(ell); return cyclic factors.
 * Set *pM to (vtilde_S[i](US[j]))_{i,j} */
static GEN
CL_tilde(GEN K, GEN US, GEN ell, GEN T, int64_t imin, GEN vdegS,
         GEN *pM, int64_t prec)
{
  int64_t i, j, k, lD, l = lg(T), lU = lg(US);
  GEN D, M, ellk;

  /* p = P^e: \tilde{Cl}(l) = (1) */
  if (l == 2) { *pM = cgetg(1, t_MAT); return cgetg(1, t_VEC); }
  M = cgetg(lU, t_MAT);
  for (j = 1; j < lU; j++)
  {
    GEN c = cgetg(l, t_COL), a = gel(US,j);
    for (i = 1; i < l; i++)
      gel(c,i) = vtilde(K, a, gel(T,i), gel(vdegS,i), ell, prec);
    gel(M,j) = c;
  }
  k = padicprec(M, ell); ellk = powiu(ell, k);
  *pM = M = gmod(M, ellk);
  M = ZM_hnfmodid(rowsplice(M, imin), ellk);
  D = matsnf0(M, 4); lD = lg(D);
  if (lD > 1 && Z_pval(gel(D,1), ell) >= k) return NULL;
  return D;
}

/* [L:K] = ell^k; return 1 if L/K is locally cyclotomic at ell, 0 otherwise */
int64_t
rnfislocalcyclo(GEN rnf)
{
  pari_sp av = avma;
  GEN K, L, S, SK, TK, SLs, SL2, TL, ell;
  ulong ll;
  int64_t i, j, k, lk, lSK;
  checkrnf(rnf);
  lk = rnf_get_degree(rnf);
  if (lk == 1) return 1;
  k = uisprimepower(lk, &ll);
  if (!k) pari_err_IMPL("rnfislocalcyclo for non-l-extensions");
  ell = utoi(ll);
  K = rnf_get_nf(rnf);
  L = rnf_build_nfabs(rnf, nf_get_prec(K));
  S = rnfidealprimedec(rnf, ell);
  SK  = gel(S,1);
  SLs = gel(S,2);
  SL2 = shallowconcat1(SLs);
  TK = padicfact(K, SK, 100); lSK = lg(SK);
  TL = padicfact(L, SL2, 100);
  for (i = 1; i < lSK; i++)
  {
    int64_t eK = etilde(K, gel(SK,i), gel(TK,i));
    GEN SL = gel(SLs,i);
    int64_t lSL = lg(SL);
    for (j = 1; j < lSL; j++)
    {
      int64_t iS = gen_search(SL2, gel(SL,j), 0, (void*)&cmp_prime_over_p,
                &cmp_nodata);
      int64_t eL = etilde(L, gel(SL,j), gel(TL,iS));
      if (dvdui(eL/eK, ell)) return gc_long(av,0);
    }
  };
  return gc_long(av,1);
}

#if 0
/* Return 1 if L/Q is locally cyclotomic at ell */
static int
islocalcycloQ(GEN L, GEN ell)
{
  GEN SL = idealprimedec(L,ell), TL;
  int64_t i, lSL = lg(SL);
  TL = padicfact(L,  SL, 100);
  for (i = 1; i < lSL; i++)
  {
    int64_t eL = etilde(L, gel(SL,i), gel(TL,i));
    if (dvdui(eL,ell)) return 0;
  }
  return 1;
}
#endif

/* true nf, pr a prid */
static int64_t
nfislocalpower_i(GEN nf, GEN pr, GEN a, GEN n)
{
  int64_t v, e, t;
  GEN p, G, L;
  a = nf_to_scalar_or_basis(nf,a);
  if (!signe(n)) return isint1(a);
  v = nfvalrem(nf, a, pr, &a); if (!dvdsi(v, n)) return 0;
  p = pr_get_p(pr);
  v = Z_pvalrem(n, p, &n);
  if (!equali1(n))
  {
    GEN T, modpr = zk_to_Fq_init(nf, &pr, &T, &p);
    GEN ap = nf_to_Fq(nf, a, modpr);
    if (!Fq_ispower(ap, n, T, p)) return 0;
  }
  if (!v) return 1;
  e = pr_get_e(pr);
  if (v == 1) /* optimal formula */
    t = itos( divii(mului(e,p), subiu(p,1)) ) + 1;
  else /* straight Hensel */
    t = 2 * e * v + 1;
  G = Idealstarprk(nf, pr, t, nf_INIT);
  L = ideallogmod(nf, a, G, powiu(p, v));
  return ZV_equal0(L) || ZV_pval(L, p) >= v;
}
int64_t
nfislocalpower(GEN nf, GEN pr, GEN a, GEN n)
{
  pari_sp av = avma;
  if (typ(n) != t_INT) pari_err_TYPE("nfislocalpower",n);
  nf = checknf(nf); checkprid(pr);
  return gc_long(av, nfislocalpower_i(nf, pr, a, n));
}

/* v_ell(  exponent(D) ) */
static int64_t
ellexpo(GEN D, GEN ell) { return lg(D) == 1? 0: Z_pval(gel(D,1), ell); }

static GEN
ellsylow(GEN cyc, GEN ell)
{
  int64_t i, l;
  GEN d = cgetg_copy(cyc, &l);
  for (i = 1; i < l; i++)
  {
    GEN c = gel(cyc,i), a;
    if (!Z_pvalrem(c, ell, &a)) break;
    gel(d,i) = diviiexact(c, a);
  }
  setlg(d, i); return d;
}

static int64_t
vnorm_x(GEN nf, GEN x, GEN ell)
{
  x = nf_to_scalar_or_alg(nf,x);
  if (typ(x) != t_POL) return 0;
  x = Q_primpart(x);
  return Q_pval(nfnorm(nf,x), ell);
}
static int64_t
vtilde_prec_x(GEN nf, GEN x, GEN ell)
{
  int64_t i, l, v;
  GEN G;
  if (typ(x) != t_MAT) return vnorm_x(nf,x,ell);
  G = gel(x,1); l = lg(G); v = 0;
  for (i = 1; i < l; i++) v = maxss(v, vnorm_x(nf,gel(G,i),ell));
  return v;
}
/* upper bound for \delta(vec): estimate loss of accuracy when evaluating
 * \tilde{v} on the vec[i] */
static int64_t
vtilde_prec(GEN nf, GEN vec, GEN ell)
{
  int64_t v0 = 0, i, l = lg(vec);
  for (i = 1; i < l; i++)
    v0 = maxss(v0, vtilde_prec_x(nf, gel(vec,i), ell));
  return 3 + v0 + z_pval(nf_get_degree(nf), ell);
}
static GEN
get_Ftilde(GEN nf, GEN S, GEN T, GEN ell, int64_t *pimin)
{
  int64_t j, lS = lg(S), vmin = lS;
  GEN Ftilde = cgetg(lS, t_VECSMALL);
  *pimin = 1;
  for (j = 1; j < lS; j++)
  {
    int64_t f = ftilde(nf, gel(S,j), gel(T,j)), v = z_pval(f, ell);
    Ftilde[j] = f; if (v < vmin) { vmin = v; *pimin = j; }
  }
  return Ftilde;
}
static GEN
bnflog_i(GEN bnf, GEN ell)
{
  int64_t prec0, prec;
  GEN nf, US, vdegS, S, T, M, CLp, CLt, Ftilde, vtG, ellk;
  GEN D, Ap, cycAp, fu;
  int64_t imin, i, j, lvAp;

  bnf = checkbnf(bnf); nf = bnf_get_nf(bnf);
  S = idealprimedec(nf, ell);
  US = sunits_mod_units(bnf, S);
  prec0 = maxss(30, vtilde_prec(nf, US, ell));
  if (!(fu = bnf_build_cheapfu(bnf)) && !(fu = bnf_compactfu(bnf)))
    bnf_build_units(bnf);
  US = shallowconcat(fu, US);
  settyp(US, t_COL);
  T = padicfact(nf, S, prec0);
  Ftilde = get_Ftilde(nf, S, T, ell, &imin);
  CLp = CL_prime(bnf, ell, S);
  cycAp = gel(CLp,1);
  Ap = gel(CLp,2);
  for(;;)
  {
    vdegS = get_vdegS(Ftilde, ell, prec0);
    CLt = CL_tilde(nf, US, ell, T, imin, vdegS, &vtG, prec0);
    if (CLt) break;
    prec0 <<= 1;
    T = padicfact(nf, S, prec0);
  }
  prec = ellexpo(cycAp, ell) + ellexpo(CLt,ell) + 1;
  if (prec == 1) return mkvec3(cgetg(1,t_VEC), cgetg(1,t_VEC), cgetg(1,t_VEC));

  ellk = powiu(ell, prec);
  lvAp = lg(Ap);
  if (lvAp > 1)
  {
    int64_t lS = lg(S);
    GEN Kcyc = bnf_get_cyc(bnf);
    GEN C = zeromatcopy(lvAp-1, lS-1);
    GEN Rell = gel(CLp,3), Uell = gel(CLp,4), ordS = gel(CLp,5);
    for (i = 1; i < lvAp; i++)
    {
      GEN a, b, bi, A = gel(Ap,i), d = gel(cycAp,i);
      bi = isprincipal(bnf, A);
      a = vecmodii(ZC_Z_mul(bi,d), Kcyc);
      /* a in subgroup generated by S = Rell; hence b integral */
      b = hnf_invimage(Rell, a);
      b = vecmodii(ZM_ZC_mul(Uell, ZC_neg(b)), ordS);
      A = mkvec2(A, trivial_fact());
      A = idealpowred(nf, A, d);
      /* find a principal representative of A_i^cycA_i up to elements of S */
      a = isprincipalfact(bnf,gel(A,1),S,b,nf_GENMAT|nf_FORCE);
      if (!gequal0(gel(a,1))) pari_err_BUG("bnflog");
      a = famat_mul_shallow(gel(A,2), gel(a,2)); /* principal part */
      if (lg(a) == 1) continue;
      for (j = 1; j < lS; j++)
        gcoeff(C,i,j) = vtilde(nf, a, gel(T,j), gel(vdegS,j), ell, prec0);
    }
    C = gmod(gneg(C),ellk);
    C = shallowtrans(C);
    M = mkmat2(mkcol2(diagonal_shallow(cycAp), C), mkcol2(gen_0, vtG));
    M = shallowmatconcat(M); /* relation matrix */
  }
  else
    M = vtG;
  M = ZM_hnfmodid(M, ellk);
  D = matsnf0(M, 4);
  if (lg(D) == 1 || !dvdii(gel(D,1), ellk))
    pari_err_BUG("bnflog [missing Z_l component]");
  D = vecslice(D,2,lg(D)-1);
  return mkvec3(D, CLt, ellsylow(cycAp, ell));
}
GEN
bnflog(GEN bnf, GEN ell)
{
  pari_sp av = avma;
  return gerepilecopy(av, bnflog_i(bnf, ell));
}

GEN
bnflogef(GEN nf, GEN pr)
{
  pari_sp av = avma;
  int64_t e, f, ef;
  GEN p;
  checkprid(pr); p = pr_get_p(pr);
  nf = checknf(nf);
  e = pr_get_e(pr);
  f = pr_get_f(pr); ef = e*f;
  if (u_pval(ef, p))
  {
    GEN T = gel(factorpadic(nf_get_pol(nf), p, 100), 1);
    int64_t j = get_ZpX_index(nf, pr, T);
    e = etilde(nf, pr, gel(T,j));
    f = ef / e;
  }
  set_avma(av); return mkvec2s(e,f);
}

GEN
bnflogdegree(GEN nf, GEN A, GEN ell)
{
  pari_sp av = avma;
  GEN AZ, A0Z, NA0;
  int64_t vAZ;

  if (typ(ell) != t_INT) pari_err_TYPE("bnflogdegree", ell);
  nf = checknf(nf);
  A = idealhnf(nf, A);
  AZ = gcoeff(A,1,1);
  vAZ = Z_pvalrem(AZ, ell, &A0Z);
  if (is_pm1(A0Z))
    NA0 = gen_1;
  else
    (void)Z_pvalrem(idealnorm(nf,A), ell, &NA0);
  if (vAZ)
  {
    GEN Aell = ZM_hnfmodid(A, powiu(ell,vAZ));
    GEN S = idealprimedec(nf, ell), T;
    int64_t l, i, s = 0;
    T = padicfact(nf, S, 100);
    l = lg(S);
    for (i = 1; i < l; i++)
    {
      GEN P = gel(S,i);
      int64_t v = idealval(nf, Aell, P);
      if (v) s += v * ftilde(nf, P, gel(T,i));
    }
    if (s) NA0 = gmul(NA0, gpowgs(ell1(ell), s));
  }
  return gerepileupto(av, NA0);
}
