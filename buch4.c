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

/*******************************************************************/
/*                                                                 */
/*               S-CLASS GROUP AND NORM SYMBOLS                    */
/*          (Denis Simon, desimon@math.u-bordeaux.fr)              */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "int.h"

/* p > 2, T ZX, p prime, x t_INT */
static int64_t
lemma6(GEN T, GEN p, int64_t nu, GEN x)
{
  int64_t la, mu;
  pari_sp av = avma;
  GEN gpx, gx = poleval(T, x);

  if (Zp_issquare(gx, p)) return gc_long(av,1);

  la = Z_pval(gx, p);
  gpx = poleval(ZX_deriv(T), x);
  mu = signe(gpx)? Z_pval(gpx,p)
                 : la+nu+1; /* mu = +oo */
  set_avma(av);
  if (la > mu<<1) return 1;
  if (la >= nu<<1 && mu >= nu) return 0;
  return -1;
}
/* p = 2, T ZX, x t_INT: return 1 = yes, -1 = no, 0 = inconclusive */
static int64_t
lemma7(GEN T, int64_t nu, GEN x)
{
  int64_t odd4, la, mu;
  pari_sp av = avma;
  GEN gpx, oddgx, gx = poleval(T, x);

  if (Zp_issquare(gx,gen_2)) return 1;

  gpx = poleval(ZX_deriv(T), x);
  la = Z_lvalrem(gx, 2, &oddgx);
  odd4 = umodiu(oddgx,4); set_avma(av);

  mu = vali(gpx);
  if (mu < 0) mu = la+nu+1; /* mu = +oo */

  if (la > mu<<1) return 1;
  if (nu > mu)
  {
    int64_t mnl = mu+nu-la;
    if (odd(la)) return -1;
    if (mnl==1) return 1;
    if (mnl==2 && odd4==1) return 1;
  }
  else
  {
    int64_t nu2 = nu << 1;
    if (la >= nu2) return 0;
    if (la == nu2 - 2 && odd4==1) return 0;
  }
  return -1;
}

/* T a ZX, p a prime, pnu = p^nu, x0 t_INT */
static int64_t
zpsol(GEN T, GEN p, int64_t nu, GEN pnu, GEN x0)
{
  int64_t i, res;
  pari_sp av = avma, btop;
  GEN x, pnup;

  res = absequaliu(p,2)? lemma7(T,nu,x0): lemma6(T,p,nu,x0);
  if (res== 1) return 1;
  if (res==-1) return 0;
  x = x0; pnup = mulii(pnu,p);
  btop = avma;
  for (i=0; i < itos(p); i++)
  {
    x = addii(x,pnu);
    if (zpsol(T,p,nu+1,pnup,x)) return gc_long(av,1);
    if (gc_needed(btop, 2))
    {
      x = gerepileupto(btop, x);
      if (DEBUGMEM > 1)
        pari_warn(warnmem, "hyperell_locally_soluble: %ld/%Ps",i,p);
    }
  }
  return gc_long(av,0);
}

/* return 1 if equation y^2=T(x) has a rational p-adic solution (possibly
 * infinite), 0 otherwise. */
int64_t
hyperell_locally_soluble(GEN T,GEN p)
{
  pari_sp av = avma;
  int64_t res;
  if (typ(T)!=t_POL) pari_err_TYPE("hyperell_locally_soluble",T);
  if (typ(p)!=t_INT) pari_err_TYPE("hyperell_locally_soluble",p);
  RgX_check_ZX(T, "hyperell_locally_soluble");
  res = zpsol(T,p,0,gen_1,gen_0) || zpsol(RgX_recip_shallow(T), p, 1, p, gen_0);
  return gc_long(av, res);
}

/* is t a square in (O_K/pr) ? Assume v_pr(t) = 0 */
static int64_t
quad_char(GEN nf, GEN t, GEN pr)
{
  GEN ord, ordp, T, p, modpr = zk_to_Fq_init(nf, &pr,&T,&p);
  t = nf_to_Fq(nf,t,modpr);
  if (T)
  {
    ord = subiu( pr_norm(pr), 1 ); /* |(O_K / pr)^*| */
    ordp= subiu( p, 1);            /* |F_p^*|        */
    t = Fq_pow(t, diviiexact(ord, ordp), T,p); /* in F_p^* */
    if (typ(t) == t_POL)
    {
      if (degpol(t)) pari_err_BUG("nfhilbertp");
      t = gel(t,2);
    }
  }
  return kronecker(t, p);
}
/* quad_char(x), x in Z, nonzero mod p */
static int64_t
Z_quad_char(GEN x, GEN pr)
{
  int64_t f = pr_get_f(pr);
  if (!odd(f)) return 1;
  return kronecker(x, pr_get_p(pr));
}

/* (pr,2) = 1. return 1 if x in Z_K is a square in Z_{K_pr}, 0 otherwise.
 * modpr = zkmodprinit(nf,pr) */
static int64_t
psquarenf(GEN nf,GEN x,GEN pr,GEN modpr)
{
  pari_sp av = avma;
  GEN p = pr_get_p(pr);
  int64_t v;

  x = nf_to_scalar_or_basis(nf, x);
  if (typ(x) == t_INT) {
    if (!signe(x)) return 1;
    v = Z_pvalrem(x, p, &x) * pr_get_e(pr);
    if (v&1) return 0;
    v = (Z_quad_char(x, pr) == 1);
  } else {
    v = ZC_nfvalrem(x, pr, &x);
    if (v&1) return 0;
    v = (quad_char(nf, x, modpr) == 1);
  }
  return gc_long(av,v);
}

static int
ZV_iseven(GEN zlog)
{
  int64_t i, l = lg(zlog);
  for (i = 1; i < l; i++)
    if (mpodd(gel(zlog,i))) return 0;
  return 1;
}

/* pr | 2, project to principal units (trivializes later discrete log) */
static GEN
to_principal_unit(GEN nf, GEN x, GEN pr, GEN sprk)
{
  if (pr_get_f(pr) != 1)
  {
    GEN prk = gel(sprk,3);
    x = nfpowmodideal(nf, x, gmael(sprk,5,1), prk);
  }
  return x;
}
/* pr | 2. Return 1 if x in Z_K is square in Z_{K_pr}, 0 otherwise */
static int
psquare2nf(GEN nf, GEN x, GEN pr, GEN sprk)
{
  int64_t v = nfvalrem(nf, x, pr, &x);
  if (v == LONG_MAX) return 1; /* x = 0 */
  /* (x,pr) = 1 */
  if (odd(v)) return 0;
  x = to_principal_unit(nf, x, pr, sprk); /* = 1 mod pr */
  return ZV_iseven(sprk_log_prk1(nf, x, sprk));
}

/*
For z in nf, z != 0.
quadratic characters modulo the prime ideal pr in nf.
pr output by nfmodprinit
pstar output by idealstar (only for p | 2).
For p odd, the output is a vector [v,c]*Mod(1,2), where
v = valuation(z,pr)
c = !issquare( z/pr^v mod pr)
For p | 2, the output is similar, with a longer sequence of 0,1 for c.
*/

GEN
nf_quadchar_modpr(GEN nf, GEN z, GEN modpr, GEN pstar)
{
  pari_sp av = avma;
  GEN pr = modpr_get_pr(modpr);
  GEN v = stoi(nfvalrem(nf, z, pr, &z));
  if( equaliu(pr_get_p(pr),2))
  {
    GEN c = ideallog(nf, z, pstar);
    return gerepilecopy(av, shallowconcat(v, shallowtrans(c)));
  }
  else
  {
    GEN c = quad_char(nf, z, modpr)==1? gen_0: gen_1;
    return gerepilecopy(av, mkvec2(v,c));
  }
}

/* pr above an odd prime */
static int64_t
lemma6nf(GEN nf, GEN T, GEN pr, int64_t nu, GEN x, GEN modpr)
{
  pari_sp av = avma;
  int64_t la, mu;
  GEN gpx, gx = nfpoleval(nf, T, x);

  if (psquarenf(nf,gx,pr,modpr)) return 1;

  la = nfval(nf,gx,pr);
  gpx = nfpoleval(nf, RgX_deriv(T), x);
  mu = gequal0(gpx)? la+nu+1 /* +oo */: nfval(nf,gpx,pr);
  set_avma(av);
  if (la > (mu<<1)) return 1;
  if (la >= (nu<<1) && mu >= nu) return 0;
  return -1;
}
/* pr above 2 */
static int64_t
lemma7nf(GEN nf, GEN T, GEN pr, int64_t nu, GEN x, GEN sprk)
{
  int64_t i, res, la, mu, q, e, v;
  GEN M, y, gpx, loggx = NULL, gx = nfpoleval(nf, T, x);

  la = nfvalrem(nf, gx, pr, &gx); /* gx /= pi^la, pi a pr-uniformizer */
  if (la == LONG_MAX) return 1;
  if (!odd(la))
  {
    gx = to_principal_unit(nf, gx, pr, sprk); /* now 1 mod pr */
    loggx = sprk_log_prk1(nf, gx, sprk); /* cheap */
    if (ZV_iseven(loggx)) return 1;
  }
  gpx = nfpoleval(nf, RgX_deriv(T), x);
  mu = gequal0(gpx)? la+nu+1 /* oo */: nfval(nf,gpx,pr);

  if (la > (mu << 1)) return 1;
  if (nu > mu)
  {
    if (odd(la)) return -1;
    q = mu+nu-la; res = 1;
  }
  else
  {
    q = (nu << 1) - la;
    if (q <= 0) return 0;
    if (odd(la)) return -1;
    res = 0;
  }
  /* la even */
  e = pr_get_e(pr);
  if (q > e << 1)  return -1;
  if (q == 1) return res;

  /* gx = 1 mod pr; square mod pi^q ? */
  v = nfvalrem(nf, nfadd(nf, gx, gen_m1), pr, &y);
  if (v >= q) return res;
  /* 1 + pi^v y = (1 + pi^vz z)^2 mod pr^q ? v < q <= 2e => vz < e => vz = vy/2
   * => y = x^2 mod pr^(min(q-v, e+v/2)), (y,pr) = 1 */
  if (odd(v)) return -1;
  /* e > 1 */
  M = cgetg(2*e+1 - q + 1, t_VEC);
  for (i = q+1; i <= 2*e+1; i++) gel(M, i-q) = sprk_log_gen_pr(nf, sprk, i);
  M = ZM_hnfmodid(shallowconcat1(M), gen_2);
  return hnf_solve(M,loggx)? res: -1;
}
/* zinit either a sprk (pr | 2) or a modpr structure (pr | p odd).
   pnu = pi^nu, pi a uniformizer */
static int64_t
zpsolnf(GEN nf,GEN T,GEN pr,int64_t nu,GEN pnu,GEN x0,GEN repr,GEN zinit)
{
  int64_t i, res;
  pari_sp av = avma;
  GEN pnup;

  res = typ(zinit) == t_VEC? lemma7nf(nf,T,pr,nu,x0,zinit)
                           : lemma6nf(nf,T,pr,nu,x0,zinit);
  set_avma(av);
  if (res== 1) return 1;
  if (res==-1) return 0;
  pnup = nfmul(nf, pnu, pr_get_gen(pr));
  nu++;
  for (i=1; i<lg(repr); i++)
  {
    GEN x = nfadd(nf, x0, nfmul(nf,pnu,gel(repr,i)));
    if (zpsolnf(nf,T,pr,nu,pnup,x,repr,zinit)) return gc_long(av,1);
  }
  return gc_long(av,0);
}

/* Let y = copy(x); y[k] := j; return y */
static GEN
ZC_add_coeff(GEN x, int64_t k, int64_t j)
{ GEN y = shallowcopy(x); gel(y, k) = utoipos(j); return y; }

/* system of representatives for Zk/pr */
static GEN
repres(GEN nf, GEN pr)
{
  int64_t f = pr_get_f(pr), N = nf_get_degree(nf), p = itos(pr_get_p(pr));
  int64_t i, j, k, pi, pf = upowuu(p, f);
  GEN perm = pr_basis_perm(nf, pr), rep = cgetg(pf+1,t_VEC);

  gel(rep,1) = zerocol(N);
  for (pi=i=1; i<=f; i++,pi*=p)
  {
    int64_t t = perm[i];
    for (j=1; j<p; j++)
      for (k=1; k<=pi; k++) gel(rep, j*pi+k) = ZC_add_coeff(gel(rep,k), t, j);
  }
  return rep;
}

/* = 1 if equation y^2 = z^deg(T) * T(x/z) has a pr-adic rational solution
 * (possibly (1,y,0) = oo), 0 otherwise.
 * coeffs of T are algebraic integers in nf */
static int64_t
locally_soluble(GEN nf,GEN T,GEN pr)
{
  GEN repr, zinit;

  if (typ(T)!=t_POL) pari_err_TYPE("nf_hyperell_locally_soluble",T);
  if (gequal0(T)) return 1;
  checkprid(pr); nf = checknf(nf);
  if (absequaliu(pr_get_p(pr), 2))
  { /* tough case */
    zinit = log_prk_init(nf, pr, 1+2*pr_get_e(pr), NULL);
    if (psquare2nf(nf,constant_coeff(T),pr,zinit)) return 1;
    if (psquare2nf(nf, leading_coeff(T),pr,zinit)) return 1;
  }
  else
  {
    zinit = zkmodprinit(nf, pr);
    if (psquarenf(nf,constant_coeff(T),pr,zinit)) return 1;
    if (psquarenf(nf, leading_coeff(T),pr,zinit)) return 1;
  }
  repr = repres(nf,pr); /* FIXME: inefficient if Npr is large */
  return zpsolnf(nf, T, pr, 0, gen_1, gen_0, repr, zinit) ||
         zpsolnf(nf, RgX_recip_shallow(T), pr, 1, pr_get_gen(pr),
                 gen_0, repr, zinit);
}
int64_t
nf_hyperell_locally_soluble(GEN nf,GEN T,GEN pr)
{
  pari_sp av = avma;
  return gc_long(av, locally_soluble(nf, T, pr));
}

/* return a * denom(a)^2, as an 'liftalg' */
static GEN
den_remove(GEN nf, GEN a)
{
  GEN da;
  a = nf_to_scalar_or_basis(nf, a);
  switch(typ(a))
  {
    case t_INT: return a;
    case t_FRAC: return mulii(gel(a,1), gel(a,2));
    case t_COL:
      a = Q_remove_denom(a, &da);
      if (da) a = ZC_Z_mul(a, da);
      return nf_to_scalar_or_alg(nf, a);
    default: pari_err_TYPE("nfhilbert",a);
      return NULL;/*LCOV_EXCL_LINE*/
  }
}

static int64_t
hilb2nf(GEN nf,GEN a,GEN b,GEN p)
{
  pari_sp av = avma;
  GEN pol;
  a = den_remove(nf, a);
  b = den_remove(nf, b);
  pol = mkpoln(3, a, gen_0, b);
  /* varn(nf.pol) = 0, pol is not a valid GEN  [as in Pol([x,x], x)].
   * But it is only used as a placeholder, hence it is not a problem */
  return gc_long(av, nf_hyperell_locally_soluble(nf,pol,p)? 1: -1);
}

/* local quadratic Hilbert symbol (a,b)_pr, for a,b (nonzero) in nf */
static int64_t
nfhilbertp(GEN nf, GEN a, GEN b, GEN pr)
{
  GEN t;
  int64_t va, vb;
  pari_sp av = avma;

  if (absequaliu(pr_get_p(pr), 2)) return hilb2nf(nf,a,b,pr);

  /* pr not above 2, compute t = tame symbol */
  va = nfval(nf,a,pr);
  vb = nfval(nf,b,pr);
  if (!odd(va) && !odd(vb)) return gc_long(av,1);
  /* Trick: pretend the exponent is 2, result is OK up to squares ! */
  t = famat_makecoprime(nf, mkvec2(a,b), mkvec2s(vb, -va),
                        pr, pr_hnf(nf, pr), gen_2);
  /* quad. symbol is image of t = (-1)^(v(a)v(b)) a^v(b) b^(-v(a))
   * by the quadratic character  */
  switch(typ(t))
  {
    default: /* t_COL */
      if (!ZV_isscalar(t)) break;
      t = gel(t,1); /* fall through */
    case t_INT:
      if (odd(va) && odd(vb)) t = negi(t);
      return gc_long(av,  Z_quad_char(t, pr));
  }
  if (odd(va) && odd(vb)) t = ZC_neg(t);
  return gc_long(av, quad_char(nf, t, pr));
}

/* Global quadratic Hilbert symbol (a,b):
 *  =  1 if X^2 - aY^2 - bZ^2 has a point in projective plane
 *  = -1 otherwise
 * a, b should be nonzero */
int64_t
nfhilbert(GEN nf, GEN a, GEN b)
{
  pari_sp av = avma;
  int64_t i, l;
  GEN S, S2, Sa, Sb, sa, sb;

  nf = checknf(nf);
  a = nf_to_scalar_or_basis(nf, a);
  b = nf_to_scalar_or_basis(nf, b);
  /* local solutions in real completions ? [ error in nfsign if arg is 0 ]*/
  sa = nfsign(nf, a);
  sb = nfsign(nf, b); l = lg(sa);
  for (i=1; i<l; i++)
    if (sa[i] && sb[i])
    {
      if (DEBUGLEVEL>3)
        err_printf("nfhilbert not soluble at real place %ld\n",i);
      return gc_long(av,-1);
    }

  /* local solutions in finite completions ? (pr | 2ab)
   * primes above 2 are toughest. Try the others first */
  Sa = idealfactor(nf, a);
  Sb = idealfactor(nf, b);
  S2 = idealfactor(nf, gen_2);
  S = merge_factor(Sa, Sb, (void*)&cmp_prime_ideal, &cmp_nodata);
  S = merge_factor(S,  S2, (void*)&cmp_prime_ideal, &cmp_nodata);
  S = gel(S,1);
  /* product of all hilbertp is 1 ==> remove one prime (above 2!) */
  for (i=lg(S)-1; i>1; i--)
    if (nfhilbertp(nf,a,b,gel(S,i)) < 0)
    {
      if (DEBUGLEVEL>3)
        err_printf("nfhilbert not soluble at finite place %Ps\n",S[i]);
      return gc_long(av,-1);
    }
  return gc_long(av,1);
}

int64_t
nfhilbert0(GEN nf,GEN a,GEN b,GEN p)
{
  nf = checknf(nf);
  if (p) {
    checkprid(p);
    if (gequal0(a)) pari_err_DOMAIN("nfhilbert", "a", "=", gen_0, a);
    if (gequal0(b)) pari_err_DOMAIN("nfhilbert", "b", "=", gen_0, b);
    return nfhilbertp(nf,a,b,p);
  }
  return nfhilbert(nf,a,b);
}

static void
p_append(GEN p, hashtable *H, hashtable *H2)
{
  ulong h = H->hash(p);
  hashentry *e = hash_search2(H, (void*)p, h);
  if (!e)
  {
    hash_insert2(H, (void*)p, NULL, h);
    if (H2) hash_insert2(H2, (void*)p, NULL, h);
  }
}

/* N a t_INT */
static void
Zfa_append(GEN N, hashtable *H, hashtable *H2)
{
  if (!is_pm1(N))
  {
    GEN v = gel(absZ_factor(N),1);
    int64_t i, l = lg(v);
    for (i=1; i<l; i++) p_append(gel(v,i), H, H2);
  }
}
/* N a t_INT or t_FRAC or ideal in HNF*/
static void
fa_append(GEN N, hashtable *H, hashtable *H2)
{
  switch(typ(N))
  {
    case t_INT:
      Zfa_append(N,H,H2);
      break;
    case t_FRAC:
      Zfa_append(gel(N,1),H,H2);
      Zfa_append(gel(N,2),H,H2);
      break;
    default: /*t_MAT*/
      Zfa_append(gcoeff(N,1,1),H,H2);
      break;
  }
}

/* apply lift(rnfeltup) to all coeffs, without rnf structure */
static GEN
nfX_eltup(GEN nf, GEN rnfeq, GEN x)
{
  int64_t i, l;
  GEN y = cgetg_copy(x, &l), zknf = nf_nfzk(nf, rnfeq);
  for (i=2; i<l; i++) gel(y,i) = nfeltup(nf, gel(x,i), zknf);
  y[1] = x[1]; return y;
}

static hashtable *
hash_create_INT(ulong s)
{ return hash_create(s, (ulong(*)(void*))&hash_GEN,
                        (int(*)(void*,void*))&equalii, 1); }
GEN
rnfisnorminit(GEN T, GEN R, int galois)
{
  pari_sp av = avma;
  int64_t i, l, dR;
  GEN S, gen, cyc, bnf, nf, nfabs, rnfeq, bnfabs, k, polabs;
  GEN y = cgetg(9, t_VEC);
  hashtable *H = hash_create_INT(100UL);

  if (galois < 0 || galois > 2) pari_err_FLAG("rnfisnorminit");
  T = get_bnfpol(T, &bnf, &nf);
  if (!bnf) bnf = Buchall(nf? nf: T, nf_FORCE, DEFAULTPREC);
  if (!nf) nf = bnf_get_nf(bnf);

  R = get_bnfpol(R, &bnfabs, &nfabs);
  if (!gequal1(leading_coeff(R))) pari_err_IMPL("non monic relative equation");
  dR = degpol(R);
  if (dR <= 2) galois = 1;

  R = RgX_nffix("rnfisnorminit", T, R, 1);
  if (nf_get_degree(nf) == 1) /* over Q */
    rnfeq = mkvec5(R,gen_0,gen_0,T,R);
  else if (galois == 2) /* needs eltup+abstorel */
    rnfeq = nf_rnfeq(nf, R);
  else /* needs abstorel */
    rnfeq = nf_rnfeqsimple(nf, R);
  polabs = gel(rnfeq,1);
  k = gel(rnfeq,3);
  if (!bnfabs || !gequal0(k))
    bnfabs = Buchall(polabs, nf_FORCE, nf_get_prec(nf));
  if (!nfabs) nfabs = bnf_get_nf(bnfabs);

  if (galois == 2)
  {
    GEN P = polabs==R? leafcopy(R): nfX_eltup(nf, rnfeq, R);
    setvarn(P, fetch_var_higher());
    galois = !!nfroots_if_split(&nfabs, P);
    (void)delete_var();
  }

  cyc = bnf_get_cyc(bnfabs);
  gen = bnf_get_gen(bnfabs); l = lg(cyc);
  for(i=1; i<l; i++)
  {
    GEN g = gel(gen,i);
    if (ugcdiu(gel(cyc,i), dR) == 1) break;
    Zfa_append(gcoeff(g,1,1), H, NULL);
  }
  if (!galois)
  {
    GEN Ndiscrel = diviiexact(nf_get_disc(nfabs), powiu(nf_get_disc(nf), dR));
    Zfa_append(Ndiscrel, H, NULL);
  }
  S = hash_keys(H); settyp(S,t_VEC);
  gel(y,1) = bnf;
  gel(y,2) = bnfabs;
  gel(y,3) = R;
  gel(y,4) = rnfeq;
  gel(y,5) = S;
  gel(y,6) = nf_pV_to_prV(nf, S);
  gel(y,7) = nf_pV_to_prV(nfabs, S);
  gel(y,8) = stoi(galois); return gerepilecopy(av, y);
}

/* T as output by rnfisnorminit
 * if flag=0 assume extension is Galois (==> answer is unconditional)
 * if flag>0 add to S all primes dividing p <= flag
 * if flag<0 add to S all primes dividing abs(flag)

 * answer is a vector v = [a,b] such that
 * x = N(a)*b and x is a norm iff b = 1  [assuming S large enough] */
GEN
rnfisnorm(GEN T, GEN x, int64_t flag)
{
  pari_sp av = avma;
  GEN bnf, rel, R, rnfeq, nfpol;
  GEN nf, aux, H, U, Y, M, A, bnfS, sunitrel, futu, S, S1, S2, Sx;
  int64_t L, i, itu;
  hashtable *H0, *H2;
  if (typ(T) != t_VEC || lg(T) != 9)
    pari_err_TYPE("rnfisnorm [please apply rnfisnorminit()]", T);
  bnf = gel(T,1);
  rel = gel(T,2);
  bnf = checkbnf(bnf);
  rel = checkbnf(rel);
  nf = bnf_get_nf(bnf);
  x = nf_to_scalar_or_alg(nf,x);
  if (gequal0(x)) { set_avma(av); return mkvec2(gen_0, gen_1); }
  if (gequal1(x)) { set_avma(av); return mkvec2(gen_1, gen_1); }
  R = gel(T,3);
  rnfeq = gel(T,4);
  if (gequalm1(x) && odd(degpol(R)))
  { set_avma(av); return mkvec2(gen_m1, gen_1); }

  /* build set T of ideals involved in the solutions */
  nfpol = nf_get_pol(nf);
  S = gel(T,5);
  H0 = hash_create_INT(100UL);
  H2 = hash_create_INT(100UL);
  L = lg(S);
  for (i = 1; i < L; i++) p_append(gel(S,i),H0,NULL);
  S1 = gel(T,6);
  S2 = gel(T,7);
  if (flag > 0)
  {
    forprime_t T;
    ulong p;
    u_forprime_init(&T, 2, flag);
    while ((p = u_forprime_next(&T))) p_append(utoipos(p), H0,H2);
  }
  else if (flag < 0)
    Zfa_append(utoipos(-flag),H0,H2);
  /* overkill: prime ideals dividing x would be enough */
  A = idealnumden(nf, x);
  fa_append(gel(A,1), H0,H2);
  fa_append(gel(A,2), H0,H2);
  Sx = hash_keys(H2); L = lg(Sx);
  if (L > 1)
  { /* new primes */
    settyp(Sx, t_VEC);
    S1 = shallowconcat(S1, nf_pV_to_prV(nf, Sx));
    S2 = shallowconcat(S2, nf_pV_to_prV(rel, Sx));
  }

  /* computation on T-units */
  futu = shallowconcat(bnf_get_fu(rel), bnf_get_tuU(rel));
  bnfS = bnfsunit(bnf,S1,LOWDEFAULTPREC);
  sunitrel = shallowconcat(futu, gel(bnfsunit(rel,S2,LOWDEFAULTPREC), 1));

  A = lift_shallow(bnfissunit(bnf,bnfS,x));
  L = lg(sunitrel);
  itu = lg(nf_get_roots(nf))-1; /* index of torsion unit in bnfsunit(nf) output */
  M = cgetg(L+1,t_MAT);
  for (i=1; i<L; i++)
  {
    GEN u = eltabstorel(rnfeq, gel(sunitrel,i));
    gel(sunitrel,i) = u;
    u = bnfissunit(bnf,bnfS, gnorm(u));
    if (lg(u) == 1) pari_err_BUG("rnfisnorm");
    gel(u,itu) = lift_shallow(gel(u,itu)); /* lift root of 1 part */
    gel(M,i) = u;
  }
  aux = zerocol(lg(A)-1); gel(aux,itu) = utoipos( bnf_get_tuN(rel) );
  gel(M,L) = aux;
  H = ZM_hnfall(M, &U, 2);
  Y = RgM_RgC_mul(U, inverseimage(H,A));
  /* Y: sols of MY = A over Q */
  setlg(Y, L);
  aux = factorback2(sunitrel, gfloor(Y));
  x = mkpolmod(x,nfpol);
  if (!gequal1(aux)) x = gdiv(x, gnorm(aux));
  x = lift_if_rational(x);
  if (typ(aux) == t_POLMOD && degpol(nfpol) == 1)
    gel(aux,2) = lift_if_rational(gel(aux,2));
  return gerepilecopy(av, mkvec2(aux, x));
}

GEN
bnfisnorm(GEN bnf, GEN x, int64_t flag)
{
  pari_sp av = avma;
  GEN T = rnfisnorminit(pol_x(fetch_var()), bnf, flag == 0? 1: 2);
  GEN r = rnfisnorm(T, x, flag == 1? 0: flag);
  (void)delete_var();
  return gerepileupto(av,r);
}
