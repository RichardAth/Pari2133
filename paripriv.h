#pragma once
/* Copyright (C) 2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

#ifdef __cplusplus
BEGINEXTERN
#endif
/* Note regarding exporting names: 
This header file is required by gp.c but is not normally included for other
programs linking to parilib, so only functions called from gp.c or whatnow.c 
should need to be exported. This is done with the macro PARILIB_API which
expands to __declspec(dllexport) or __declspec(dllimport) when needed. */

/* for qsort */
typedef int (*QSCOMP)(const void *, const void *);

#define uel(a,i)            (((ulong*)(a))[i])
#define ucoeff(a,i,j)       (((ulong**)(a))[j][i])
#define umael(a,i,j)        (((ulong**)(a))[i][j])
#define umael2(a,i,j)       (((ulong**)(a))[i][j])
#define umael3(a,i,j,k)     (((ulong***)(a))[i][j][k])
#define umael4(a,i,j,k,l)   (((ulong****)(a))[i][j][k][l])
#define umael5(a,i,j,k,l,m) (((ulong*****)(a))[i][j][k][l][m])

#define numberof(x) (sizeof(x) / sizeof((x)[0]))

/* to manipulate 'blocs' */
#define BL_HEAD 8
#define bl_base(x) (void*)((x) - BL_HEAD)
#define bl_height(x) (((GEN)x)[-8])
#define bl_left(x)   (((GEN*)x)[-7])
#define bl_right(x)  (((GEN*)x)[-6])
#define bl_size(x)   (((GEN)x)[-5])
#define bl_refc(x)   (((GEN)x)[-4])
#define bl_next(x)   (((GEN*)x)[-3])
#define bl_prev(x)   (((GEN*)x)[-2])
#define bl_num(x)    (((GEN)x)[-1])

void clone_lock(GEN C);
void clone_unlock(GEN C);
void clone_unlock_deep(GEN C);

/* swap */
#define lswap(x,y) {int64_t _z=x; x=y; y=_z;}
#define pswap(x,y) {GEN *_z=x; x=y; y=_z;}
#define swap(x,y)  {GEN  _z=x; x=y; y=_z;}
#define dswap(x,y) { double _t=x; x=y; y=_t; }
#define pdswap(x,y) { double* _t=x; x=y; y=_t; }
#define swapspec(x,y, nx,ny) {swap(x,y); lswap(nx,ny);}

/* loops */
GEN incloop(GEN a);
GEN resetloop(GEN a, GEN b);
GEN setloop(GEN a);

/* parser */

/* GP control structures */
#define EXPR_WRAP(code, call) \
{ GEN z; GEN __E = code; \
  push_lex(gen_0, __E); z = call; pop_lex(1); return z; }
#define EXPRVOID_WRAP(code, call) \
{ GEN __E = code; \
  push_lex(gen_0, __E); call; pop_lex(1); }
#define EXPR_ARG __E, &gp_eval
#define EXPR_ARGPREC __E, &gp_evalprec
#define EXPR_ARGUPTO __E, &gp_evalupto
#define EXPR_ARGBOOL __E, &gp_evalbool
#define EXPR_ARGVOID __E, &gp_evalvoid

/* sumiter.c */
GEN  iferrpari(GEN a, GEN b, GEN c);
void forfactored(GEN a, GEN b, GEN code);
void forpari(GEN a, GEN b, GEN node);
void foreachpari(GEN a, GEN node);
void forsquarefree(GEN a, GEN b, GEN code);
void untilpari(GEN a, GEN b);
void whilepari(GEN a, GEN b);

GEN  geval_gp(GEN x, GEN t);

/* anal.c */
GEN  ifpari(GEN g, GEN a, GEN b);
GEN  andpari(GEN a, GEN b);
GEN  orpari(GEN a, GEN b);
void ifpari_void(GEN g, GEN a, GEN b);
GEN  ifpari_multi(GEN g, GEN a);
GEN  gadde(GEN *x, GEN y);
GEN  gadd1e(GEN *x);
GEN  gdive(GEN *x, GEN y);
GEN  gdivente(GEN *x, GEN y);
GEN  gdivrounde(GEN *x, GEN y);
GEN  gmode(GEN *x, GEN y);
GEN  gmule(GEN *x, GEN y);
GEN  gshiftle(GEN *x, int64_t n);
GEN  gshiftre(GEN *x, int64_t n);
GEN  gsube(GEN *x, GEN y);
GEN  gsub1e(GEN *x);
GEN  gshift_right(GEN x, int64_t n);

/* sumiter.c */
GEN  asympnum0(GEN u, GEN alpha, int64_t prec);
GEN  asympnumraw0(GEN u, int64_t LIM, GEN alpha, int64_t prec);
GEN  derivnum0(GEN a, GEN code, GEN ind, int64_t prec);
GEN  derivfun0(GEN args, GEN def, GEN code, int64_t k, int64_t prec);
GEN  direuler0(GEN a, GEN b, GEN code, GEN c);

GEN  direuler_bad(void *E, GEN (*eval)(void *, GEN, int64_t), GEN a, GEN b, GEN c, GEN Sbad);
void forcomposite(GEN a, GEN b, GEN code);
void fordiv(GEN a, GEN code);
void fordivfactored(GEN a, GEN code);
void forell0(int64_t a, int64_t b, GEN code, int64_t flag);
void forperm0(GEN k, GEN code);
void forprime(GEN a, GEN b, GEN code);
void forprimestep(GEN a, GEN b, GEN q, GEN code);
void forstep(GEN a, GEN b, GEN s, GEN code);
void forsubgroup0(GEN cyc, GEN bound, GEN code);
void forsubset0(GEN nk, GEN code);
void forvec(GEN x, GEN code, int64_t flag);
void forpart0(GEN k, GEN code , GEN nbound, GEN abound);
GEN  intcirc0(GEN a, GEN R, GEN code, GEN tab, int64_t prec);
GEN  intfuncinit0(GEN a, GEN b, GEN code, int64_t m, int64_t prec);
GEN  intnum0(GEN a, GEN b, GEN code, GEN tab, int64_t prec);
GEN  intnumgauss0(GEN a, GEN b, GEN code, GEN tab, int64_t prec);
GEN  intnumromb0_bitprec(GEN a, GEN b, GEN code, int64_t flag, int64_t bit);
GEN  laurentseries0(GEN f, int64_t M, int64_t v, int64_t prec);
GEN  limitnum0(GEN u, GEN alpha, int64_t prec);
GEN  matrice(GEN nlig, GEN ncol, GEN code);
void pariplot0(GEN a, GEN b, GEN code, GEN ysmlu, GEN ybigu, int64_t prec);
GEN  prodeuler0(GEN a, GEN b, GEN code, int64_t prec);
GEN  prodinf0(GEN a, GEN code, int64_t flag, int64_t prec);
GEN  produit(GEN a, GEN b, GEN code, GEN x);
GEN  somme(GEN a, GEN b, GEN code, GEN x);
GEN  sumalt0(GEN a, GEN code,int64_t flag, int64_t prec);
GEN  sumdivexpr(GEN num, GEN code);
GEN  sumdivmultexpr0(GEN num, GEN code);
GEN  suminf0_bitprec(GEN a, GEN code, int64_t bit);
GEN  sumnum0(GEN a, GEN code, GEN tab, int64_t prec);
GEN  sumnumap0(GEN a, GEN code, GEN tab, int64_t prec);
GEN  sumnumlagrange0(GEN a, GEN code, GEN tab, int64_t prec);
GEN  sumnummonien0(GEN a, GEN code, GEN tab, int64_t prec);
GEN  sumpos0(GEN a, GEN code, int64_t flag,int64_t prec);
GEN  vecexpr0(GEN nmax, GEN code, GEN pred);
GEN  vecexpr1(GEN nmax, GEN code, GEN pred);
GEN  vecteursmall(GEN nmax, GEN code);
GEN  vecteur(GEN nmax, GEN n);
GEN  vvecteur(GEN nmax, GEN n);
GEN  zbrent0(GEN a, GEN b, GEN code, int64_t prec);
GEN  solvestep0(GEN a, GEN b, GEN step, GEN code, int64_t flag, int64_t prec);

GEN  ploth0(GEN a, GEN b, GEN code, int64_t flag, int64_t n, int64_t prec);
GEN  plothexport0(GEN fmt, GEN a, GEN b, GEN code, int64_t flags, int64_t n, int64_t prec);
GEN  psploth0(GEN a,GEN b,GEN code,int64_t flag,int64_t n,int64_t prec);
GEN  plotrecth0(int64_t ne,GEN a,GEN b,GEN code,ulong flags,int64_t n,int64_t prec);

GEN  listcreate_gp(int64_t n);

/* mt */
void mt_sigint(void);
void mt_err_recover(int er);
void mt_export_add(const char *str, GEN val);
void mt_export_del(const char *str);
void mt_init_stack(size_t s);
int  mt_is_thread(void);

GEN  eisker_worker(GEN Ei, GEN M, GEN D, GEN co, GEN CD);
GEN  pareval_worker(GEN code);
GEN  parselect_worker(GEN d, GEN code);
void parfor0(GEN a, GEN b, GEN code, GEN code2);
GEN  parfor_worker(GEN i, GEN C);
void parforeach0(GEN x, GEN code, GEN code2);
void parforprime0(GEN a, GEN b, GEN code, GEN code2);
void parforprimestep0(GEN a, GEN b, GEN q, GEN code, GEN code2);
void parforvec0(GEN a, GEN code, GEN code2, int64_t flag);
GEN  parvector_worker(GEN i, GEN C);
GEN  polmodular_worker(GEN pt, ulong L, GEN hilb, GEN factu,
       GEN vne, GEN vinfo, int64_t compute_derivs, GEN j_powers, GEN fdb);
GEN  nmV_polint_center_tree_worker(GEN Va, GEN T, GEN R, GEN xa, GEN m2);
GEN  nmV_chinese_center_tree_seq(GEN A, GEN P, GEN T, GEN R);
GEN  nxMV_polint_center_tree_worker(GEN Va, GEN T, GEN R, GEN xa, GEN m2);
GEN  nxMV_chinese_center_tree_seq(GEN A, GEN P, GEN T, GEN R);
GEN  F2xq_log_Coppersmith_worker(GEN u, int64_t i, GEN V, GEN R);
GEN  Flxq_log_Coppersmith_worker(GEN u, int64_t i, GEN V, GEN R);
GEN  Fp_log_sieve_worker(int64_t a, int64_t prmax, GEN C, GEN c, GEN Ci, GEN ci, GEN pr, GEN sz);
GEN  QM_charpoly_ZX_worker(GEN P, GEN M, GEN dM);
GEN  QXQ_div_worker(GEN P, GEN A, GEN B, GEN C);
GEN  QXQ_inv_worker(GEN P, GEN A, GEN B);
GEN  ZX_resultant_worker(GEN P, GEN A, GEN B, GEN dB);
GEN  ZXQX_resultant_worker(GEN P, GEN A, GEN B, GEN T, GEN dB);
GEN  ZX_ZXY_resultant_worker(GEN P, GEN A, GEN B, GEN dB, GEN v);
GEN  ZX_direct_compositum_worker(GEN P, GEN A, GEN B);
GEN  ZXQX_direct_compositum_worker(GEN P, GEN A, GEN B, GEN C);
GEN  ZX_gcd_worker(GEN P, GEN A, GEN B, GEN g);
GEN  ZXQ_minpoly_worker(GEN P, GEN A, GEN B, int64_t d);
GEN  ZM_det_worker(GEN P, GEN A);
GEN  ZM_inv_worker(GEN P, GEN A);
GEN  ZM_ker_worker(GEN P, GEN A);
GEN  ZM_mul_worker(GEN P, GEN A, GEN B);
GEN  ZabM_inv_worker(GEN P, GEN A, GEN Q);
GEN  aprcl_step4_worker(ulong q, GEN pC, GEN N, GEN v);
GEN  aprcl_step6_worker(GEN r, int64_t t, GEN N, GEN N1, GEN et);
GEN  ecpp_sqrt_worker(GEN g, GEN N, GEN p);
GEN  ecpp_ispsp_worker(GEN N);
GEN  ecpp_step2_worker(GEN S, GEN HD, GEN primelist);
GEN  primecertisvalid_ecpp_worker(GEN certi);
GEN  lfuninit_worker(int64_t r, GEN K, GEN L, GEN peh2d, GEN vroots, GEN dr, GEN di, GEN an, GEN bn);
GEN  lfuninit_theta2_worker(int64_t r, GEN L, GEN qk, GEN a, GEN di, GEN an, GEN bn);
GEN  gen_parapply(GEN worker, GEN D);
GEN  parapply_slice_worker(GEN worker, GEN D);
GEN  gen_parapply_slice(GEN worker, GEN D, int64_t mmin);
GEN  gen_crt(const char *str, GEN worker, forprime_t *S, GEN dB, ulong bound, int64_t mmin, GEN *pt_mod,
             GEN crt(GEN, GEN, GEN*), GEN center(GEN, GEN, GEN));
void gen_inccrt(const char *str, GEN worker, GEN dB, int64_t n, int64_t mmin,
           forprime_t *S, GEN *pt_H, GEN *pt_mod, GEN crt(GEN, GEN, GEN*),
           GEN center(GEN, GEN, GEN));
void gen_inccrt_i(const char *str, GEN worker, GEN dB, int64_t n, int64_t mmin,
           forprime_t *S, GEN *pH, GEN *pmod, GEN crt(GEN, GEN, GEN*),
           GEN center(GEN, GEN, GEN));
GEN  direllnf_worker(GEN P, ulong X, GEN E);
GEN  dirartin_worker(GEN P, ulong X, GEN nf, GEN G, GEN V, GEN aut);
GEN  direllsympow_worker(GEN P, ulong X, GEN E, ulong m);
GEN  dirgenus2_worker(GEN P, ulong X, GEN Q);
GEN  pardireuler(GEN worker, GEN a, GEN b, GEN c, GEN Sbad);
GEN  FpM_ratlift_worker(GEN A, GEN mod, GEN B);
GEN  ellQ_factorback_worker(GEN A, GEN P, GEN L, GEN c4);
GEN  chinese_unit_worker(GEN P, GEN A, GEN U, GEN B, GEN D, GEN C);

/* Relative number fields */
enum { rnf_NFABS = 1, rnf_MAPS };

/* Finite fields */
enum { t_FF_FpXQ = 0, t_FF_Flxq = 1, t_FF_F2xq = 2 };
GEN FF_ellinit(GEN E, GEN fg);
GEN FF_elldata(GEN E, GEN fg);

/* L functions */
enum { t_LFUN_GENERIC, t_LFUN_ZETA, t_LFUN_NF, t_LFUN_ELL, t_LFUN_KRONECKER,
       t_LFUN_CHIZ, t_LFUN_CHIGEN, t_LFUN_ETA,
       t_LFUN_DIV, t_LFUN_MUL, t_LFUN_CONJ,
       t_LFUN_SYMPOW_ELL, t_LFUN_QF, t_LFUN_ARTIN, t_LFUN_MFCLOS,
       t_LFUN_GENUS2, t_LFUN_TWIST, t_LFUN_CLOSURE0, t_LFUN_SHIFT};
enum { t_LDESC_INIT, t_LDESC_THETA, t_LDESC_PRODUCT };

/* Elliptic curves */
/* common to Q and Rg */
enum { R_PERIODS = 1, R_ETA, R_ROOTS, R_AB };

enum { Qp_ROOT = 1, Qp_TATE };
enum { Q_GROUPGEN = 5, Q_GLOBALRED, Q_ROOTNO, Q_MINIMALMODEL };
enum { NF_MINIMALMODEL = 1, NF_GLOBALRED, NF_MINIMALPRIMES, NF_ROOTNO, NF_NF };

/* common to Fp and Fq */
enum { FF_CARD = 1, FF_GROUP, FF_GROUPGEN, FF_O };

/* for Buchall_param */
enum { fupb_NONE = 0, fupb_RELAT, fupb_LARGE, fupb_PRECI };

/* Polycyclic presentation for the classgroup of discriminant D */
typedef struct {
  int64_t D; /* Negative discriminant */
  int64_t h; /* Size of classgroup */
  int64_t enum_cnt; /* Either h or h/2 (if L0 is set) */
  /* If nonzero, L0=L[0] and n[0]=2 and classpoly is a perfect square
   * (and we enumerate each double root just once), default is 0 */
  int64_t L0;
  /* Product of primes L that are prohibited as norms of generators or
   * auxilliary prime forms (by default, primes that make enumeration hard) */
  int64_t Lfilter;
  /* Norms of implicit generators (primeforms a=(L*x^2+b*x*y+c*y^2)
   * with norm L and b >=0) */
  int64_t *L;
  int64_t *m; /* products of relative orders: m[i] is the order of <g_1,...,g_i> */
  int64_t *n; /* Relative orders */
  int64_t *o; /* Absolute orders */
  /* Power relations (a[i]^n[i] = a[0]^e[0]*...*a[i-1]^e[i-1], where e
   * is an exponent vector of length i stored at offset binom(i,2) of r) */
  int64_t *r;
  int64_t *orient_p; /* Optional list of norms of orienting primes p ... */
  int64_t *orient_q; /* or product of primes p*q (q=1 when only p is needed) */
  int64_t *orient_reps; /* Representation of orienting norm p*q in terms of Ls */
  int64_t inv; /* Attached invariant */
  int64_t k; /* Number of generators */
  GEN _data; /* Storage space for the above arrays */
} classgp_pcp_struct;
typedef classgp_pcp_struct classgp_pcp_t[1];

/* Represents the data in the equation(s)
 *   4p = t^2 - v^2 D = t^2 - v^2 u^2 D_K = w^2 D_K.
 * t is the absolute trace, so always > 0.
 * T is a twisting parameter, which satisfies (T|p) == -1. */
typedef struct {
  int64_t D, t, u, v;
  ulong p, pi, s2, T;
} norm_eqn_struct;
typedef norm_eqn_struct norm_eqn_t[1];

#define zv_to_longptr(v) (&((v)[1]))
#define zv_to_ulongptr(v) ((ulong *)&((v)[1]))

/* Modular invariants */
#define INV_J       0
#define INV_F       1
#define INV_F2      2
#define INV_F3      3
#define INV_F4      4
#define INV_G2      5
#define INV_W2W3    6
#define INV_F8      8
#define INV_W3W3    9
#define INV_W2W5    10
#define INV_W2W7    14
#define INV_W3W5    15
#define INV_W3W7    21
#define INV_W2W3E2  23
#define INV_W2W5E2  24
#define INV_W2W13   26
#define INV_W2W7E2  27
#define INV_W3W3E2  28
#define INV_W5W7    35
#define INV_W3W13   39

/* Get coefficient of x^d in f, assuming f is nonzero. */
INLINE ulong Flx_coeff(GEN f, int64_t d) { return f[d + 2]; }
/* Return the root of f, assuming deg(f) = 1. */
INLINE ulong Flx_deg1_root(GEN f, ulong p) {
  if (degpol(f) != 1) pari_err_BUG("Flx_deg1_root");
  return Fl_div(Fl_neg(Flx_coeff(f, 0), p), Flx_coeff(f, 1), p);
}

/* Allocation / gerepile */
int64_t   getdebugvar(void);
void   setdebugvar(int64_t n);
void   debug_stack(void);
void   fill_stack(void);
void   minim_alloc(int64_t n, double ***q, GEN *x, double **y,  double **z, double **v);
int    pop_entree_block(entree *ep, int64_t loc);
int    pop_val_if_newer(entree *ep, int64_t loc);

/* general printing */
PARILIB_API void print_errcontext(PariOUT *out, const char *msg, const char *s, const char *entry);    /* called from gp.c */
PARILIB_API void print_prefixed_text(PariOUT *out, const char *s, const char *prefix, const char *str);
INLINE void print_text(const char *s) { print_prefixed_text(pariOut, s,NULL,NULL); }    /* called from gp.c */
INLINE void out_print_text(PariOUT *out, const char *s) { print_prefixed_text(out, s,NULL,NULL); }    /* called from whatnow.c */
INLINE int64_t is_keyword_char(char c) { return (isalnum((int)c) || c=='_'); }

/* Interfaces (GP, etc.) */
hashtable *hash_from_link(GEN e, GEN names, int use_stack);
void gen_relink(GEN x, hashtable *table);
entree* do_alias(entree *ep);
char* get_sep(const char *t);
int64_t get_int(const char *s, int64_t dflt);
ulong get_uint(const char *s);

void pari_sigint(const char *s);
void* get_stack(double fraction, int64_t min);

void  free_graph(void);
void  initout(int initerr);
void  resetout(int initerr);
PARILIB_API void  init_linewrap(int64_t w);     /* called from gp.c */
void  print_functions_hash(const char *s);
PARILIB_API GEN   readbin(const char *name, FILE *f, int *vector);   /* called from gp.c */
int   term_height(void);
int   term_width(void);
/* gp_colors */
void decode_color(int64_t n, int64_t *c);


/* defaults */
extern int64_t precreal;

PARILIB_API void lim_lines_output(char *s, int64_t n, int64_t max);    /* called from gp.c */
PARILIB_API void gen_output(GEN x);          /* called from gp.c */
void fputGEN_pariout(GEN x, pariout_t *T, FILE *out);

void parsestate_reset(void);
void parsestate_save(struct pari_parsestate *state);
void parsestate_restore(struct pari_parsestate *state);

void compilestate_reset(void);
void compilestate_save(struct pari_compilestate *comp);
PARILIB_API void compilestate_restore(struct pari_compilestate *comp);     /* called from gp.c */

void filestate_save(struct pari_filestate *file);
void filestate_restore(struct pari_filestate *file);
PARILIB_API void tmp_restore(pariFILE *F);                 /* called from gp.c */

int64_t evalstate_get_trace(void);
void evalstate_set_trace(int64_t lvl);
PARILIB_API void evalstate_clone(void);                   /* called from gp.c */
void evalstate_reset(void);
void evalstate_restore(struct pari_evalstate *state);
GEN  evalstate_restore_err(struct pari_evalstate *state);
void evalstate_save(struct pari_evalstate *state);
void varstate_save(struct pari_varstate *s);
void varstate_restore(struct pari_varstate *s);

void mtstate_save(struct pari_mtstate *s);
void mtstate_reset(void);
void mtstate_restore(struct pari_mtstate *s);

void debug_context(void);

typedef struct {
  const char *s;
  size_t ls;
  char **dir;
} forpath_t;
void forpath_init(forpath_t *T, gp_path *path, const char *s);
char *forpath_next(forpath_t *T);

/* GP output && output format */
void gpwritebin(const char *s, GEN x);
extern char *current_logfile;

/* colors */
extern int64_t    gp_colors[];
PARILIB_API extern int     disable_color;      /* referenced from gp.c */

/* entrees */
#define EpVALENCE(ep) ((ep)->valence & 0xFF)
#define EpSTATIC(ep) ((ep)->valence & 0x100)
#define EpSETSTATIC(ep) ((ep)->valence |= 0x100)
enum { EpNEW = 100, EpALIAS, EpVAR, EpINSTALL };
#define initial_value(ep) ((ep)+1)

/* functions lists */
extern const int64_t functions_tblsz;  /* hashcodes table size */
extern entree **functions_hash;   /* functions hashtable */
extern entree **defaults_hash;    /* defaults hashtable */


PARILIB_API void init_filtre(filtre_t *F, Buffer *buf);      /* referenced from gp.c */


PARILIB_API extern int (*cb_pari_get_line_interactive)(const char*, const char*, filtre_t *F);
PARILIB_API extern char *(*cb_pari_fgets_interactive)(char *s, int n, FILE *f);


PARILIB_API char *pari_translate_string(const char *src, char *s, char *entry);

PARILIB_API gp_data *default_gp_data(void);         /* called from gp.c */

typedef char *(*fgets_t)(char *, int, void*);

typedef struct input_method {
/* optional */
  fgets_t myfgets;  /* like libc fgets() but last argument is (void*) */
/* mandatory */
  char * (*getline)(char**, int f, struct input_method*, filtre_t *F);
  int free; /* boolean: must we free the output of getline() ? */
/* optional */
  const char *prompt, *prompt_cont;
  void *file;  /* can be used as last argument for fgets() */
} input_method;

int input_loop(filtre_t *F, input_method *IM);
char *file_input(char **s0, int junk, input_method *IM, filtre_t *F);
char *file_getline(Buffer *b, char **s0, input_method *IM);

/* readline */
typedef struct {
  /* pointers to readline variables/functions */
  char **line_buffer;
  int *point;
  int *end;
  char **(*completion_matches)(const char *, char *(*)(const char*, int));
  char *(*filename_completion_function)(const char *, int);
  char *(*username_completion_function)(const char *, int);
  int (*insert)(int, int);
  int *completion_append_character;

  /* PARI-specific */
  int back;  /* rewind the cursor by this number of chars */
} pari_rl_interface;

/* Code which wants to use readline needs to do the following:

#include <readline.h>
#include <paripriv.h>
pari_rl_interface pari_rl;
pari_use_readline(pari_rl);

This will initialize the pari_rl structure. A pointer to this structure
must be given as first argument to all PARI readline functions. */

/* IMPLEMENTATION NOTE: this really must be a macro (not a function),
 * since we refer to readline symbols. */
#define pari_use_readline(pari_rl) do {\
    (pari_rl).line_buffer = &rl_line_buffer; \
    (pari_rl).point = &rl_point; \
    (pari_rl).end = &rl_end; \
    (pari_rl).completion_matches = &rl_completion_matches; \
    (pari_rl).filename_completion_function = &rl_filename_completion_function; \
    (pari_rl).username_completion_function = &rl_username_completion_function; \
    (pari_rl).insert = &rl_insert; \
    (pari_rl).completion_append_character = &rl_completion_append_character; \
    (pari_rl).back = 0; } while(0)

/* FIXME: EXPORT AND DOCUMENT THE FOLLOWING */

/* PROBABLY NOT IN THE RIGHT FILE, SORT BY THEME */

/* multiprecision */
GEN   adduispec_offset(ulong s, GEN x, int64_t offset, int64_t nx);
int   lgcdii(ulong* d, ulong* d1, ulong* u, ulong* u1, ulong* v, ulong* v1, ulong vmax);
ulong xxgcduu(ulong d, ulong d1, int f, ulong* u, ulong* u1, ulong* v, ulong* v1, int64_t* s);
ulong rgcduu(ulong d, ulong d1, ulong vmax, ulong* u, ulong* u1, ulong* v, ulong* v1, int64_t *s);
ulong xgcduu(ulong d, ulong d1, int f, ulong* v, ulong* v1, int64_t *s);

GEN   divgunu(GEN x, ulong i);
GEN   divrunu(GEN x, ulong i);
GEN   muliispec(GEN x, GEN y, int64_t nx, int64_t ny);
GEN   red_montgomery(GEN T, GEN N, ulong inv);
GEN   sqrispec(GEN x, int64_t nx);
ulong *convi(GEN x, int64_t *l);

int approx_0(GEN x, GEN y);

/* powers */
GEN    rpowuu(ulong a, ulong n, int64_t prec);

/* floats */
double dabs(double s, double t);
double darg(double s, double t);
void   dcxlog(double s, double t, double *a, double *b);
double dnorm(double s, double t);
double dbllog2(GEN z);

/* hnf */
GEN hnfadd(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,GEN extramat,GEN extraC);
GEN hnfadd_i(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,GEN extramat,GEN extraC);
GEN hnfspec_i(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,int64_t k0);
GEN hnfspec(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,int64_t k0);
GEN mathnfspec(GEN x, GEN *ptperm, GEN *ptdep, GEN *ptB, GEN *ptC);

GEN LLL_check_progress(GEN Bnorm, int64_t n0, GEN m, int final, int64_t *ti_LLL);

/* integer factorization / discrete log */
ulong is_kth_power(GEN x, ulong p, GEN *pt);
GEN   mpqs(GEN N);

/* Polynomials */
/* a) Arithmetic/conversions */
GEN  lift_if_rational(GEN x);
GEN  monomial(GEN a, int64_t degpol, int64_t v);
GEN  monomialcopy(GEN a, int64_t degpol, int64_t v);
GEN  ser2pol_i(GEN x, int64_t lx);
GEN  ser2rfrac_i(GEN x);
GEN  swap_vars(GEN b0, int64_t v);
GEN  RgX_recipspec_shallow(GEN x, int64_t l, int64_t n);

/* b) Modular */
GEN  bezout_lift_fact(GEN T, GEN Tmod, GEN p, int64_t e);
GEN  polsym_gen(GEN P, GEN y0, int64_t n, GEN T, GEN N);
GEN  ZXQ_charpoly_sqf(GEN A, GEN B, int64_t *lambda, int64_t v);
GEN  ZX_disc_all(GEN,ulong);
GEN  ZX_resultant_all(GEN A, GEN B, GEN dB, ulong bound);
GEN  ZX_ZXY_resultant_all(GEN A, GEN B, int64_t *lambda, GEN *LPRS);

GEN FlxqM_mul_Kronecker(GEN A, GEN B, GEN T, ulong p);
GEN FqM_mul_Kronecker(GEN x, GEN y, GEN T, GEN p);

/* c) factorization */
GEN chk_factors_get(GEN lt, GEN famod, GEN c, GEN T, GEN N);
int64_t cmbf_maxK(int64_t nb);
GEN ZX_DDF(GEN x);
GEN initgaloisborne(GEN T, GEN dn, int64_t prec, GEN *pL, GEN *pprep, GEN *pdis);
GEN quicktrace(GEN x, GEN sym);

/* pari_init / pari_close */
void pari_close_compiler(void);
void pari_close_evaluator(void);
void pari_close_files(void);
void pari_close_floats(void);
void pari_close_homedir(void);
void pari_close_parser(void);
void pari_close_paths(void);
void pari_close_primes(void);

void pari_init_compiler(void);
PARILIB_API void pari_init_defaults(void);        /* called from gp.c */
void pari_init_evaluator(void);
void pari_init_files(void);
void pari_init_floats(void);
void pari_init_homedir(void);
void pari_init_graphics(void);
void pari_init_parser(void);
void pari_init_rand(void);
PARILIB_API void pari_init_paths(void);         /* called from gp.c */
void pari_init_primetab(void);
void pari_init_seadata(void);
GEN pari_get_seadata(void);
void pari_set_primetab(GEN global_primetab);
void pari_set_seadata(GEN seadata);
void pari_set_varstate(int64_t *vp, struct pari_varstate *vs);
void pari_thread_close_files(void);

void export_add(const char *str, GEN val);
void export_del(const char *str);
GEN  export_get(const char *str);
void exportall(void);
void unexportall(void);

/* BY FILES */

/* parinf.h */

GEN fincke_pohst(GEN a,GEN BOUND,int64_t stockmax,int64_t PREC, FP_chk_fun *CHECK);
void init_zlog(zlog_S *S, GEN bid);
GEN  log_gen_arch(zlog_S *S, int64_t index);
GEN  log_gen_pr(zlog_S *S, int64_t index, GEN nf, int64_t e);
GEN  sprk_log_gen_pr(GEN nf, GEN sprk, int64_t e);
GEN  sprk_log_prk1(GEN nf, GEN a, GEN sprk);
GEN    poltobasis(GEN nf,GEN x);
GEN    coltoalg(GEN nf,GEN x);

GEN    rnfdisc_get_T(GEN nf, GEN P, GEN *lim);
GEN    make_integral(GEN nf, GEN L0, GEN f, GEN listpr);
GEN    rnfallbase(GEN nf, GEN pol, GEN lim, GEN rnfeq, GEN *pD, GEN *pfi,
                  GEN *pdKP);
GEN    subgroupcondlist(GEN cyc, GEN bound, GEN listKer);

/* Qfb.c */

GEN     redimagsl2(GEN q, GEN *U);
GEN     redrealsl2(GEN V, GEN d, GEN rd);
GEN     redrealsl2step(GEN A, GEN d, GEN rd);

/* alglin1.c */

typedef int64_t (*pivot_fun)(GEN,GEN,int64_t,GEN);
GEN ZM_pivots(GEN x0, int64_t *rr);
GEN RgM_pivots(GEN x0, GEN data, int64_t *rr, pivot_fun pivot);
void RgMs_structelim_col(GEN M, int64_t nbcol, int64_t nbrow, GEN A, GEN *p_col, GEN *p_lin);

/* arith1.c */

int     is_gener_Fp(GEN x, GEN p, GEN p_1, GEN L);
int     is_gener_Fl(ulong x, ulong p, ulong p_1, GEN L);

/* arith2.c */

int     divisors_init(GEN n, GEN *pP, GEN *pE);
int64_t    set_optimize(int64_t what, GEN g);

/* base1.c */

GEN zk_galoisapplymod(GEN nf, GEN z, GEN S, GEN p);

/* base2.c */

GEN     dim1proj(GEN prh);
GEN     gen_if_principal(GEN bnf, GEN x);

/* base3.c */

void    check_nfelt(GEN x, GEN *den);
GEN     zk_ei_mul(GEN nf, GEN x, int64_t i);
GEN     log_prk(GEN nf, GEN a, GEN sprk, GEN mod);
GEN     log_prk_units(GEN nf, GEN D, GEN sprk);
GEN     log_prk_units_init(GEN bnf);
GEN     veclog_prk(GEN nf, GEN v, GEN sprk);
GEN     log_prk_init(GEN nf, GEN pr, int64_t k, GEN mod);

/* base4.c */

GEN     factorbackprime(GEN nf, GEN L, GEN e);

/* bb_group.c */

GEN     producttree_scheme(int64_t n);

/* bern.c */
int64_t bernbitprec(int64_t N);

/* bibli2.c */

GEN sort_factor_pol(GEN y, int (*cmp)(GEN,GEN));

/* buch1.c */

int64_t   bnf_increase_LIMC(int64_t LIMC, int64_t LIMCMAX);

/* buch2.c */

typedef struct GRHprime_t { ulong p; double logp; GEN dec; } GRHprime_t;
typedef struct GRHcheck_t { double cD, cN; GRHprime_t *primes; int64_t clone, nprimes, maxprimes; ulong limp; forprime_t P; } GRHcheck_t;
void    free_GRHcheck(GRHcheck_t *S);
void    init_GRHcheck(GRHcheck_t *S, int64_t N, int64_t R1, double LOGD);
void    GRH_ensure(GRHcheck_t *S, int64_t nb);
ulong   GRH_last_prime(GRHcheck_t *S);
int     GRHok(GRHcheck_t *S, double L, double SA, double SB);
GEN     extract_full_lattice(GEN x);
GEN     init_red_mod_units(GEN bnf, int64_t prec);
GEN     isprincipalarch(GEN bnf, GEN col, GEN kNx, GEN e, GEN dx, int64_t *pe);
GEN     red_mod_units(GEN col, GEN z);

/* buch3.c */

GEN     minkowski_bound(GEN D, int64_t N, int64_t r2, int64_t prec);
int     subgroup_conductor_ok(GEN H, GEN L);
GEN     subgrouplist_cond_sub(GEN bnr, GEN C, GEN bound);

/* buch4.c */

GEN     nf_quadchar_modpr(GEN nf, GEN z, GEN modpr, GEN pstar);

/* crvwtors.c */

void random_curves_with_m_torsion(ulong *a4, ulong *a6, ulong *tx, ulong *ty, int64_t ncurves, int64_t m, ulong p);

/* dirichlet.c */
GEN direuler_factor(GEN s, int64_t n);

/* elliptic.c */

void ellprint(GEN e);
GEN  elltors_psylow(GEN e, ulong p);
GEN  ellintegralbmodel(GEN e, GEN *pv);
GEN  ellQ_genreduce(GEN E, GEN G, int64_t prec);

/* es.c */

PARILIB_API void    killallfiles(void);                                 /* called from gp.c */
PARILIB_API pariFILE* newfile(FILE *f, const char *name, int type);     /* called from gp.c */
PARILIB_API int     popinfile(void);                                    /* called from gp.c */
pariFILE* try_pipe(const char *cmd, int flag);

/* F2m.c */

GEN     F2m_gauss_pivot(GEN x, int64_t *rr);
GEN     F2m_gauss_sp(GEN a, GEN b);
GEN     F2m_invimage_i(GEN A, GEN B);

/* Fle.c */

void    FleV_add_pre_inplace(GEN P, GEN Q, GEN a4, ulong p, ulong pi);
void    FleV_dbl_pre_inplace(GEN P, GEN a4, ulong p, ulong pi);
void    FleV_mulu_pre_inplace(GEN P, ulong n, GEN a4, ulong p, ulong pi);
void    FleV_sub_pre_inplace(GEN P, GEN Q, GEN a4, ulong p, ulong pi);

/* Flv.c */

GEN     Flm_gauss_sp(GEN a, GEN b, ulong *detp, ulong p);
GEN     Flm_invimage_i(GEN A, GEN B, ulong p);
GEN     Flm_inv_sp(GEN a, ulong *detp, ulong p);
GEN     Flm_pivots(GEN x, ulong p, int64_t *rr, int64_t inplace);

/* Flxq_log.c */

GEN Flxq_log_index(GEN a0, GEN b0, GEN m, GEN T0, ulong p);
int Flxq_log_use_index(GEN m, GEN T0, ulong p);

/* FlxqE.c */

GEN     ZpXQ_norm_pcyc(GEN x, GEN T, GEN q, GEN p);
int64_t    zx_is_pcyc(GEN T);

/* FpV.c */

GEN FpMs_leftkernel_elt_col(GEN M, int64_t nbcol, int64_t nbrow, GEN p);
GEN FpX_to_mod_raw(GEN z, GEN p);

/* FpX.c */

GEN     ZlXQXn_expint(GEN h, int64_t e, GEN T, GEN p, ulong pp);

/* FpX_factor.c */

GEN     ddf_to_ddf2(GEN V);
int64_t    ddf_to_nbfact(GEN D);
GEN     vddf_to_simplefact(GEN V, int64_t d);

/* FpXQX_factor.c */

GEN     FpXQX_factor_Berlekamp(GEN x, GEN T, GEN p);

/* forprime.c*/

void    init_modular_big(forprime_t *S);
void    init_modular_small(forprime_t *S);

/* galconj.c */

GEN     galoiscosets(GEN O, GEN perm);
GEN     matrixnorm(GEN M, int64_t prec);

/* gen1.c */

GEN     gred_rfrac_simple(GEN n, GEN d);
GEN     sqr_ser_part(GEN x, int64_t l1, int64_t l2);

/* hash.c */

hashtable *hashstr_import_static(hashentry *e, ulong size);

/* hyperell.c */

GEN     ZlXQX_hyperellpadicfrobenius(GEN H, GEN T, ulong p, int64_t n);
GEN     hyperell_redsl2(GEN Q);

/* ifactor1.c */

ulong snextpr(ulong p, byteptr *d, int64_t *rcn, int64_t *q, int (*ispsp)(ulong));

/* intnum.c */

GEN     contfraceval_inv(GEN CF, GEN tinv, int64_t nlim);

/* mftrace.c */

void pari_close_mf(void);
int64_t polishomogeneous(GEN P);
GEN sertocol(GEN S);

/* prime.c */

int64_t    BPSW_psp_nosmalldiv(GEN N);
int     MR_Jaeschke(GEN n);
int64_t    isanypower_nosmalldiv(GEN N, GEN *px);
void    prime_table_next_p(ulong a, byteptr *pd, ulong *pp, ulong *pn);

/* perm.c */

int64_t    cosets_perm_search(GEN C, GEN p);
GEN     perm_generate(GEN S, GEN H, int64_t o);
int64_t    perm_relorder(GEN p, GEN S);
GEN     vecperm_extendschreier(GEN C, GEN v, int64_t n);

/* polclass.c */

GEN polclass0(int64_t D, int64_t inv, int64_t xvar, GEN *db);

/* polmodular.c */

GEN polmodular0_ZM(int64_t L, int64_t inv, GEN J, GEN Q, int compute_derivs, GEN *db);
GEN Flm_Fl_polmodular_evalx(GEN phi, int64_t L, ulong j, ulong p, ulong pi);
GEN polmodular_db_init(int64_t inv);
void polmodular_db_clear(GEN db);
void polmodular_db_add_level(GEN *db, int64_t L, int64_t inv);
void polmodular_db_add_levels(GEN *db, int64_t *levels, int64_t k, int64_t inv);
GEN polmodular_db_for_inv(GEN db, int64_t inv);
GEN polmodular_db_getp(GEN fdb, int64_t L, ulong p);

int64_t modinv_level(int64_t inv);
int64_t modinv_degree(int64_t *p1, int64_t *p2, int64_t inv);
int64_t modinv_ramified(int64_t D, int64_t inv);
int64_t modinv_j_from_2double_eta(GEN F, int64_t inv, ulong x0, ulong x1, ulong p, ulong pi);
GEN double_eta_raw(int64_t inv);
ulong modfn_root(ulong j, norm_eqn_t ne, int64_t inv);
int64_t modfn_unambiguous_root(ulong *r, int64_t inv, ulong j0, norm_eqn_t ne, GEN jdb);
GEN qfb_nform(int64_t D, int64_t n);

/* Fle.c */

ulong   Flj_order_ufact(GEN P, ulong n, GEN F, ulong a4, ulong p, ulong pi);

/* polarit3.c */

GEN     Flm_Frobenius_pow(GEN M, int64_t d, GEN T, ulong p);
GEN     FpM_Frobenius_pow(GEN M, int64_t d, GEN T, GEN p);
GEN     FpX_compositum(GEN A, GEN B, GEN p);
GEN     Flx_direct_compositum(GEN A, GEN B, ulong p);
GEN     FlxV_direct_compositum(GEN V, ulong p);
GEN     FlxqX_direct_compositum(GEN P, GEN Q, GEN T, ulong p);
GEN     FpX_direct_compositum(GEN A, GEN B, GEN p);
GEN     FpXV_direct_compositum(GEN V, GEN p);
GEN     nf_direct_compositum(GEN nf, GEN A, GEN B);
ulong   ZX_ZXY_ResBound(GEN A, GEN B, GEN dB);
GEN     ffinit_Artin_Schreier(ulong p, int64_t l);
GEN     ffinit_rand(GEN p, int64_t n);

/* readline.c */

char**  pari_completion(pari_rl_interface *pari_rl, char *text, int START, int END);
char**  pari_completion_matches(pari_rl_interface *pari_rl, const char *s, int64_t pos, int64_t *wordpos);

/* RgX.c */

GEN     RgX_homogenous_evalpow(GEN P, GEN A, GEN B);
GEN     QXQX_homogenous_evalpow(GEN P, GEN A, GEN B, GEN T);

/* subcyclo.c */

GEN     galoiscyclo(int64_t n, int64_t v);
GEN     znstar_bits(int64_t n, GEN H);
int64_t    znstar_conductor(GEN H);
int64_t    znstar_conductor_bits(GEN bits);
GEN     znstar_cosets(int64_t n, int64_t phi_n, GEN H);
GEN     znstar_elts(int64_t n, GEN H);
GEN     znstar_generate(int64_t n, GEN V);
GEN     znstar_hnf(GEN Z, GEN M);
GEN     znstar_hnf_elts(GEN Z, GEN H);
GEN     znstar_hnf_generators(GEN Z, GEN M);
GEN     znstar_reduce_modulus(GEN H, int64_t n);
GEN     znstar_small(GEN zn);

/* trans1.c */

struct abpq { GEN *a, *b, *p, *q; };
struct abpq_res { GEN P, Q, B, T; };
void    abpq_init(struct abpq *A, int64_t n);
void    abpq_sum(struct abpq_res *r, int64_t n1, int64_t n2, struct abpq *A);
GEN     logagmcx(GEN q, int64_t prec);
GEN     zellagmcx(GEN a0, GEN b0, GEN r, GEN t, int64_t prec);

/* trans2.c */

GEN     trans_fix_arg(int64_t *prec, GEN *s0, GEN *sig, GEN *tau, pari_sp *av, GEN *res);

/* trans3.c */

GEN     double_eta_quotient(GEN a, GEN w, GEN D, int64_t p, int64_t q, GEN pq, GEN sqrtD);
GEN     inv_szeta_euler(int64_t n, int64_t prec);

/* volcano.c */

int64_t j_level_in_volcano(GEN phi, ulong j, ulong p, ulong pi, int64_t L, int64_t depth);
ulong ascend_volcano(GEN phi, ulong j, ulong p, ulong pi, int64_t level, int64_t L, int64_t depth, int64_t steps);
ulong descend_volcano(GEN phi, ulong j, ulong p, ulong pi, int64_t level, int64_t L, int64_t depth, int64_t steps);
int64_t next_surface_nbr(ulong *nJ, GEN phi, int64_t L, int64_t h, ulong J, const ulong *pJ, ulong p, ulong pi);
GEN enum_roots(ulong j, norm_eqn_t ne, GEN fdb, classgp_pcp_t G);

#ifdef __cplusplus
ENDEXTERN
#endif
