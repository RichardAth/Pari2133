/* Copyright (C) 2000-2003  The PARI group.

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
/*        INITIALIZING THE SYSTEM, ERRORS, STACK MANAGEMENT        */
/*                                                                 */
/*******************************************************************/
/* _GNU_SOURCE is needed before first include to get RUSAGE_THREAD */
#undef _GNU_SOURCE /* avoid warning */
//#define _GNU_SOURCE
#include <string.h>

#if defined(_WIN32) || defined(__CYGWIN32__)
#  include <process.h>
#endif

#include "pari.h"
#include "int.h"
#define INLINE static

#if defined(_WIN32) || defined(__CYGWIN32__)
#include "mingw.h"
#endif

#if defined(STACK_CHECK) && !defined(__EMX__) && !defined(_WIN32)
#  include <sys/types.h>
#  include <sys/time.h>
#  include <sys/resource.h>
#endif
#if defined(HAS_WAITPID) && defined(HAS_SETSID)
#  include <sys/wait.h>
#endif
#ifdef HAS_MMAP
#  include <sys/mman.h>
#endif
#if defined(USE_GETTIMEOFDAY) || defined(USE_GETRUSAGE) || defined(USE_TIMES)
//#  include <sys/time.h>
#endif
#if defined(USE_GETRUSAGE)
#  include <sys/resource.h>
#endif
#if defined(USE_FTIME) || defined(USE_FTIMEFORWALLTIME)
#  include <sys/timeb.h>
#endif
#if defined(USE_CLOCK_GETTIME) || defined(USE_TIMES)
#  include <time.h>
#endif
#if defined(USE_TIMES)
#  include <time.h>
#endif
//#include "mpi.h"


#include "paripriv.h"
#include "anal.h"


const double LOG10_2 = 0.3010299956639812; /* log_10(2) */
const double LOG2_10 = 3.321928094887362;  /* log_2(10) */


GEN gnil, gen_0, gen_1, gen_m1, gen_2, gen_m2, ghalf, err_e_STACK;

static const ulong readonly_constants[] = {
  evaltyp(t_INT) | _evallg(2),  /* gen_0 */
  evallgefint(2),
  evaltyp(t_INT) | _evallg(2),  /* gnil */
  evallgefint(2),
  evaltyp(t_INT) | _evallg(3),  /* gen_1 */
  evalsigne(1) | evallgefint(3),
  1,
  evaltyp(t_INT) | _evallg(3),  /* gen_2 */
  evalsigne(1) | evallgefint(3),
  2,
  evaltyp(t_INT) | _evallg(3),  /* gen_m1 */
  evalsigne(-1) | evallgefint(3),
  1,
  evaltyp(t_INT) | _evallg(3),  /* gen_m2 */
  evalsigne(-1) | evallgefint(3),
  2,
  evaltyp(t_ERROR) | _evallg(2), /* err_e_STACK */
  e_STACK,
  evaltyp(t_FRAC) | _evallg(3), /* ghalf */
  (ulong)(readonly_constants+4),
  (ulong)(readonly_constants+7)
};
THREAD GEN zetazone, bernzone, primetab;
byteptr diffptr;
PARILIB_API FILE    *pari_outfile, *pari_errfile, *pari_logfile, *pari_infile;
char    *current_logfile, *current_psfile, *pari_datadir;
int64_t gp_colors[c_LAST];
PARILIB_API int     disable_color;
ulong   DEBUGFILES, DEBUGLEVEL, DEBUGMEM;
int64_t DEBUGVAR;
ulong   pari_mt_nbthreads;
int64_t precreal;
ulong   precdl;
PARILIB_API ulong pari_logstyle;
PARILIB_API gp_data *GP_DATA;

entree  **varentries;
THREAD int64_t *varpriority;

PARILIB_API THREAD pari_sp avma;


THREAD struct pari_mainstack *pari_mainstack;

static void ** MODULES;
static pari_stack s_MODULES;
const int64_t functions_tblsz = 135; /* size of functions_hash */
entree **functions_hash, **defaults_hash;

char *(*cb_pari_fgets_interactive)(char *s, int n, FILE *f);
int (*cb_pari_get_line_interactive)(const char*, const char*, filtre_t *F);
void (*cb_pari_quit)(int64_t);
void (*cb_pari_init_histfile)(void);
void (*cb_pari_ask_confirm)(const char *);
int  (*cb_pari_handle_exception)(int64_t);
int  (*cb_pari_err_handle)(GEN);
int  (*cb_pari_whatnow)(PariOUT *out, const char *, int);
PARILIB_API void (*cb_pari_sigint)(void);
void (*cb_pari_pre_recover)(int64_t);
PARILIB_API void (*cb_pari_err_recover)(int64_t);
int (*cb_pari_break_loop)(int);
int (*cb_pari_is_interactive)(void);
void (*cb_pari_start_output)();

const char * pari_library_path = NULL;

static THREAD GEN global_err_data;
THREAD jmp_buf *iferr_env;
const int64_t CATCH_ALL = -1;

static void pari_init_timer(void);

/*********************************************************************/
/*                                                                   */
/*                       BLOCKS & CLONES                             */
/*                                                                   */
/*********************************************************************/
/*#define DEBUG*/
static THREAD int64_t next_block;
static THREAD GEN cur_block; /* current block in block list */
static THREAD GEN root_block; /* current block in block list */
#ifdef DEBUG
static THREAD int64_t NUM;
#endif

static void
pari_init_blocks(void)
{
  next_block = 0; cur_block = NULL; root_block = NULL;
#ifdef DEBUG
  NUM = 0;
#endif
}

static void
pari_close_blocks(void)
{
  while (cur_block) killblock(cur_block);
}

static int64_t
blockheight(GEN bl) { return bl? bl_height(bl): 0; }

static int64_t
blockbalance(GEN bl)
{ return bl ? blockheight(bl_left(bl)) - blockheight(bl_right(bl)): 0; }

static void
fix_height(GEN bl)
{ bl_height(bl) = maxss(blockheight(bl_left(bl)), blockheight(bl_right(bl)))+1; }

static GEN
bl_rotright(GEN y)
{
  GEN x = bl_left(y), t = bl_right(x);
  bl_right(x) = y;
  bl_left(y)  = t;
  fix_height(y);
  fix_height(x);
  return x;
}

static GEN
bl_rotleft(GEN x)
{
  GEN y = bl_right(x), t = bl_left(y);
  bl_left(y)  = x;
  bl_right(x) = t;
  fix_height(x);
  fix_height(y);
  return y;
}

static GEN
blockinsert(GEN x, GEN bl, int64_t *d)
{
  int64_t b, c;
  if (!bl)
  {
    bl_left(x)=NULL; bl_right(x)=NULL;
    bl_height(x)=1; return x;
  }
  c = cmpuu((ulong)x, (ulong)bl);
  if (c < 0)
    bl_left(bl) = blockinsert(x, bl_left(bl), d);
  else if (c > 0)
    bl_right(bl) = blockinsert(x, bl_right(bl), d);
  else return bl; /* ??? Already exist in the tree ? */
  fix_height(bl);
  b = blockbalance(bl);
  if (b > 1)
  {
    if (*d > 0) bl_left(bl) = bl_rotleft(bl_left(bl));
    return bl_rotright(bl);
  }
  if (b < -1)
  {
    if (*d < 0) bl_right(bl) = bl_rotright(bl_right(bl));
    return bl_rotleft(bl);
  }
  *d = c; return bl;
}

static GEN
blockdelete(GEN x, GEN bl)
{
  int64_t b;
  if (!bl) return NULL; /* ??? Do not exist in the tree */
  if (x < bl)
    bl_left(bl) = blockdelete(x, bl_left(bl));
  else if (x > bl)
    bl_right(bl) = blockdelete(x, bl_right(bl));
  else
  {
    if (!bl_left(bl) && !bl_right(bl)) return NULL;
    else if (!bl_left(bl)) return bl_right(bl);
    else if (!bl_right(bl)) return bl_left(bl);
    else
    {
      GEN r = bl_right(bl);
      while (bl_left(r)) r = bl_left(r);
      bl_right(r) = blockdelete(r, bl_right(bl));
      bl_left(r) = bl_left(bl);
      bl = r;
    }
  }
  fix_height(bl);
  b = blockbalance(bl);
  if (b > 1)
  {
    if (blockbalance(bl_left(bl)) >= 0) return bl_rotright(bl);
    else
    { bl_left(bl) = bl_rotleft(bl_left(bl)); return bl_rotright(bl); }
  }
  if (b < -1)
  {
    if (blockbalance(bl_right(bl)) <= 0) return bl_rotleft(bl);
    else
    { bl_right(bl) = bl_rotright(bl_right(bl)); return bl_rotleft(bl); }
  }
  return bl;
}

static GEN
blocksearch(GEN x, GEN bl)
{
  if (isclone(x)) return x;
  if (isonstack(x) || is_universal_constant(x)) return NULL;
  while (bl)
  {
    if (x >= bl  && x < bl + bl_size(bl))
      return bl;
    bl = x < bl ? bl_left(bl): bl_right(bl);
  }
  return NULL; /* Unknown address */
}

void
clone_lock(GEN x)
{
  GEN y = blocksearch(x, root_block);
  if (y && isclone(y))
  {
    if (DEBUGMEM > 2)
      err_printf("locking block no %ld: %08lx from %08lx\n", bl_num(y), y, x);
    ++bl_refc(y);
  }
}

void
clone_unlock(GEN x)
{
  GEN y = blocksearch(x, root_block);
  if (y && isclone(y))
  {
    if (DEBUGMEM > 2)
      err_printf("unlocking block no %ld: %08lx from %08lx\n", bl_num(y), y, x);
    gunclone(y);
  }
}

void
clone_unlock_deep(GEN x)
{
  GEN y = blocksearch(x, root_block);
  if (y && isclone(y))
  {
    if (DEBUGMEM > 2)
      err_printf("unlocking deep block no %ld: %08lx from %08lx\n", bl_num(y), y, x);
    gunclone_deep(y);
  }
}

/* Return x, where:
 * x[-8]: AVL height
 * x[-7]: adress of left child or NULL
 * x[-6]: adress of right child or NULL
 * x[-5]: size
 * x[-4]: reference count
 * x[-3]: adress of next block
 * x[-2]: adress of preceding block.
 * x[-1]: number of allocated blocs.
 * x[0..n-1]: malloc-ed memory. */
GEN
newblock(size_t n)
{
  int64_t d = 0;
  int64_t *x = (int64_t *) pari_malloc((n + BL_HEAD)*sizeof(int64_t)) + BL_HEAD;

  bl_size(x) = n;
  bl_refc(x) = 1;
  bl_next(x) = NULL;
  bl_prev(x) = cur_block;
  bl_num(x)  = next_block++;
  if (cur_block) bl_next(cur_block) = x;
  root_block = blockinsert(x, root_block, &d);
#ifdef DEBUG
  err_printf("+ %ld\n", ++NUM);
#endif
  if (DEBUGMEM > 2)
    err_printf("new block, size %6lu (no %ld): %08lx\n", n, next_block-1, x);
  return cur_block = x;
}

GEN
gcloneref(GEN x)
{
  if (isclone(x)) { ++bl_refc(x); return x; }
  else return gclone(x);
}

void
gclone_refc(GEN x) { ++bl_refc(x); }

void
gunclone(GEN x)
{
  if (--bl_refc(x) > 0) return;
  BLOCK_SIGINT_START;
  root_block = blockdelete(x, root_block);
  if (bl_next(x)) bl_prev(bl_next(x)) = bl_prev(x);
  else
  {
    cur_block = bl_prev(x);
    next_block = bl_num(x);
  }
  if (bl_prev(x)) bl_next(bl_prev(x)) = bl_next(x);
  if (DEBUGMEM > 2)
    err_printf("killing block (no %ld): %08lx\n", bl_num(x), x);
  free((void*)bl_base(x)); /* pari_free not needed: we already block */
  BLOCK_SIGINT_END;
#ifdef DEBUG
  err_printf("- %ld\n", NUM--);
#endif
}

/* Recursively look for clones in the container and kill them. Then kill
 * container if clone. SIGINT could be blocked until it returns */
void
gunclone_deep(GEN x)
{
  int64_t i, lx;
  GEN v;
  if (isclone(x) && bl_refc(x) > 1) { --bl_refc(x); return; }
  BLOCK_SIGINT_START;
  switch(typ(x))
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x);
      for (i=1;i<lx;i++) gunclone_deep(gel(x,i));
      break;
    case t_LIST:
      v = list_data(x); lx = v? lg(v): 1;
      for (i=1;i<lx;i++) gunclone_deep(gel(v,i));
      if (v) killblock(v);
      break;
  }
  if (isclone(x)) gunclone(x);
  BLOCK_SIGINT_END;
}

int
pop_entree_block(entree *ep, int64_t loc)
{
  GEN x = (GEN)ep->value;
  if (bl_num(x) < loc) return 0; /* older */
  if (DEBUGMEM>2)
    err_printf("popping %s (block no %ld)\n", ep->name, bl_num(x));
  gunclone_deep(x); return 1;
}

/***************************************************************************
 **                                                                       **
 **                           Export                                      **
 **                                                                       **
 ***************************************************************************/

static hashtable *export_hash;
static void
pari_init_export(void)
{
  export_hash = hash_create_str(1,0);
}
static void
pari_close_export(void)
{
  hash_destroy(export_hash);
}

/* Exported values are blocks, but do not have the clone bit set so that they
 * are not affected by clone_lock and ensure_nb, etc. */

void
export_add(const char *str, GEN val)
{
  hashentry *h;
  val = gclone(val); unsetisclone(val);
  h = hash_search(export_hash, (void*) str);
  if (h)
  {
    GEN v = (GEN)h->val;
    h->val = val;
    setisclone(v); gunclone(v);
  }
  else
    hash_insert(export_hash,(void*)str, (void*) val);
}

void
export_del(const char *str)
{
  hashentry *h = hash_remove(export_hash,(void*)str);
  if (h)
  {
    GEN v = (GEN)h->val;
    setisclone(v); gunclone(v);
    pari_free(h);
  }
}

GEN
export_get(const char *str)
{
  return hash_haskey_GEN(export_hash,(void*)str);
}

void
unexportall(void)
{
  pari_sp av = avma;
  GEN keys = hash_keys(export_hash);
  int64_t i, l = lg(keys);
  for (i = 1; i < l; i++) mt_export_del((const char *)keys[i]);
  set_avma(av);
}

void
exportall(void)
{
  int64_t i;
  for (i = 0; i < functions_tblsz; i++)
  {
    entree *ep;
    for (ep = functions_hash[i]; ep; ep = ep->next)
      if (EpVALENCE(ep)==EpVAR) mt_export_add(ep->name, (GEN)ep->value);
  }
}

/*********************************************************************/
/*                                                                   */
/*                       C STACK SIZE CONTROL                        */
/*                                                                   */
/*********************************************************************/
/* Avoid core dump on deep recursion. Adapted Perl code by Dominic Dunlop */
THREAD void *PARI_stack_limit = NULL;

#ifdef STACK_CHECK

#  ifdef __EMX__                                /* Emulate */
void
pari_stackcheck_init(void *pari_stack_base)
{
  if (!pari_stack_base) { PARI_stack_limit = NULL; return; }
  PARI_stack_limit = get_stack(1./16, 32*1024);
}
#  elif _WIN32
void
pari_stackcheck_init(void *pari_stack_base)
{
  ulong size = 1ULL << 21;
  if (!pari_stack_base) { PARI_stack_limit = NULL; return; }
  if (size > (ulong)pari_stack_base)
    PARI_stack_limit = (void*)(((ulong)pari_stack_base) / 16);
  else
    PARI_stack_limit = (void*)((ulong)pari_stack_base - (size/16)*15);
}
#  else /* !__EMX__ && !_WIN32 */
/* Set PARI_stack_limit to (a little above) the lowest safe address that can be
 * used on the stack. Leave PARI_stack_limit at its initial value (NULL) to
 * show no check should be made [init failed]. Assume stack grows downward. */
void
pari_stackcheck_init(void *pari_stack_base)
{
  struct rlimit rip;
  ulong size;
  if (!pari_stack_base) { PARI_stack_limit = NULL; return; }
  if (getrlimit(RLIMIT_STACK, &rip)) return;
  size = rip.rlim_cur;
  if (size == (ulong)RLIM_INFINITY || size > (ulong)pari_stack_base)
    PARI_stack_limit = (void*)(((ulong)pari_stack_base) / 16);
  else
    PARI_stack_limit = (void*)((ulong)pari_stack_base - (size/16)*15);
}
#  endif /* !__EMX__ */

#else
void
pari_stackcheck_init(void *pari_stack_base)
{
  (void) pari_stack_base; PARI_stack_limit = NULL;
}
#endif /* STACK_CHECK */

/*******************************************************************/
/*                         HEAP TRAVERSAL                          */
/*******************************************************************/
struct getheap_t { int64_t n, l; };
/* x is a block, not necessarily a clone [x[0] may not be set] */
static void
f_getheap(GEN x, void *D)
{
  struct getheap_t *T = (struct getheap_t*)D;
  T->n++;
  T->l += bl_size(x) + BL_HEAD;
}
GEN
getheap(void)
{
  struct getheap_t T = { 0, 0 };
  traverseheap(&f_getheap, &T); return mkvec2s(T.n, T.l);
}

static void
traverseheap_r(GEN bl, void(*f)(GEN, void *), void *data)
{
  if (!bl) return;
  traverseheap_r(bl_left(bl), f, data);
  traverseheap_r(bl_right(bl), f, data);
  f(bl, data);
}

void
traverseheap( void(*f)(GEN, void *), void *data)
{
  traverseheap_r(root_block,f, data);
}

/*********************************************************************/
/*                          DAEMON / FORK                            */
/*********************************************************************/
#if defined(HAS_WAITPID) && defined(HAS_SETSID)
/* Properly fork a process, detaching from main process group without creating
 * zombies on exit. Parent returns 1, son returns 0 */
int
pari_daemon(void)
{
  pid_t pid = fork();
  switch(pid) {
      case -1: return 1; /* father, fork failed */
      case 0:
        (void)setsid(); /* son becomes process group leader */
        if (fork()) _exit(0); /* now son exits, also when fork fails */
        break; /* grandson: its father is the son, which exited,
                * hence father becomes 'init', that'll take care of it */
      default: /* father, fork succeeded */
        (void)waitpid(pid,NULL,0); /* wait for son to exit, immediate */
        return 1;
  }
  /* grandson. The silly '!' avoids a gcc-8 warning (unused value) */
  (void)!freopen("/dev/null","r",stdin);
  return 0;
}
#else
int
pari_daemon(void)
{
  pari_err_IMPL("pari_daemon without waitpid & setsid");
  return 0;
}
#endif

/*********************************************************************/
/*                                                                   */
/*                       SYSTEM INITIALIZATION                       */
/*                                                                   */
/*********************************************************************/
static int try_to_recover = 0;
PARILIB_API THREAD VOLATILE int PARI_SIGINT_block = 0, PARI_SIGINT_pending = 0;

/*********************************************************************/
/*                         SIGNAL HANDLERS                           */
/*********************************************************************/
static void
dflt_sigint_fun(void) { pari_err(e_MISC, "user interrupt"); }

#if defined(_WIN32) || defined(__CYGWIN32__)
PARILIB_API int win32ctrlc = 0, win32alrm = 0;
void
dowin32ctrlc(void)
{
  win32ctrlc = 0;
  cb_pari_sigint();
}
#endif

static void
pari_handle_SIGINT(void)
{
#ifdef _WIN32
  if (++win32ctrlc >= 5) _exit(3);
#else
  cb_pari_sigint();
#endif
}

typedef void (*pari_sighandler_t)(int);

pari_sighandler_t
os_signal(int sig, pari_sighandler_t f)
{
#ifdef HAS_SIGACTION
  struct sigaction sa, oldsa;

  sa.sa_handler = f;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_NODEFER;

  if (sigaction(sig, &sa, &oldsa)) return NULL;
  return oldsa.sa_handler;
#else
  return signal(sig,f);
#endif
}

void
pari_sighandler(int sig)
{
  const char *msg;
#ifndef HAS_SIGACTION
  /*SYSV reset the signal handler in the handler*/
  (void)os_signal(sig,pari_sighandler);
#endif
  switch(sig)
  {
#ifdef SIGBREAK
    case SIGBREAK:
      if (PARI_SIGINT_block==1)
      {
        PARI_SIGINT_pending=SIGBREAK;
        mt_sigint();
      }
      else pari_handle_SIGINT();
      return;
#endif

#ifdef SIGINT
    case SIGINT:
      if (PARI_SIGINT_block==1)
      {
        PARI_SIGINT_pending=SIGINT;
        mt_sigint();
      }
      else pari_handle_SIGINT();
      return;
#endif

#ifdef SIGSEGV
    case SIGSEGV:
      msg="PARI/GP (Segmentation Fault)"; break;
#endif
#ifdef SIGBUS
    case SIGBUS:
      msg="PARI/GP (Bus Error)"; break;
#endif
#ifdef SIGFPE
    case SIGFPE:
      msg="PARI/GP (Floating Point Exception)"; break;
#endif

#ifdef SIGPIPE
    case SIGPIPE:
    {
      pariFILE *f = GP_DATA->pp->file;
      if (f && pari_outfile == f->file)
      {
        GP_DATA->pp->file = NULL; /* to avoid oo recursion on error */
        pari_outfile = stdout; pari_fclose(f);
        pari_err(e_MISC, "Broken Pipe, resetting file stack...");
      }
      return; /* LCOV_EXCL_LINE */
    }
#endif

    default: msg="signal handling"; break;
  }
  pari_err_BUG(msg);
}

void
pari_sig_init(void (*f)(int))
{
#ifdef SIGBUS
  (void)os_signal(SIGBUS,f);
#endif
#ifdef SIGFPE
  (void)os_signal(SIGFPE,f);
#endif
#ifdef SIGINT
  (void)os_signal(SIGINT,f);
#endif
#ifdef SIGBREAK
  (void)os_signal(SIGBREAK,f);
#endif
#ifdef SIGPIPE
  (void)os_signal(SIGPIPE,f);
#endif
#ifdef SIGSEGV
  (void)os_signal(SIGSEGV,f);
#endif
}

/*********************************************************************/
/*                      STACK AND UNIVERSAL CONSTANTS                */
/*********************************************************************/
static void
init_universal_constants(void)
{
  gen_0  = (GEN)readonly_constants;
  gnil   = (GEN)readonly_constants+2;
  gen_1  = (GEN)readonly_constants+4;
  gen_2  = (GEN)readonly_constants+7;
  gen_m1 = (GEN)readonly_constants+10;
  gen_m2 = (GEN)readonly_constants+13;
  err_e_STACK = (GEN)readonly_constants+16;
  ghalf  = (GEN)readonly_constants+18;
}

static void
pari_init_errcatch(void)
{
  iferr_env = NULL;
  global_err_data = NULL;
}

/*********************************************************************/
/*                           INIT DEFAULTS                           */
/*********************************************************************/
void
pari_init_defaults(void)
{
  int64_t i;
  initout(1);

#ifdef LONG_IS_64BIT
  precreal = 128;
#else
  precreal = 96;
#endif

  precdl = 16;
  DEBUGFILES = DEBUGLEVEL = 0;
  DEBUGMEM = 1;
  disable_color = 1;
  pari_logstyle = logstyle_none;

  current_psfile = pari_strdup("pari.ps");
  current_logfile= pari_strdup("pari.log");
  pari_logfile = NULL;

  pari_datadir = os_getenv("GP_DATA_DIR");
  if (!pari_datadir)
  {
#if defined(_WIN32) || defined(__CYGWIN32__)
    if (paricfg_datadir[0]=='@' && paricfg_datadir[1]==0)
      pari_datadir = win32_datadir();
    else
#endif
      pari_datadir = pari_strdup(paricfg_datadir);
  }
  else pari_datadir= pari_strdup(pari_datadir);
  for (i=0; i<c_LAST; i++) gp_colors[i] = c_NONE;
}

/*********************************************************************/
/*                   FUNCTION HASHTABLES, MODULES                    */
/*********************************************************************/
extern entree functions_basic[], functions_default[];
static void
pari_init_functions(void)
{
  pari_stack_init(&s_MODULES, sizeof(*MODULES),(void**)&MODULES);
  pari_stack_pushp(&s_MODULES,functions_basic);
  functions_hash = (entree**) pari_calloc(sizeof(entree*)*functions_tblsz);
  pari_fill_hashtable(functions_hash, functions_basic);
  defaults_hash = (entree**) pari_calloc(sizeof(entree*)*functions_tblsz);
  pari_add_defaults_module(functions_default);
}

void
pari_add_module(entree *ep)
{
  pari_fill_hashtable(functions_hash, ep);
  pari_stack_pushp(&s_MODULES, ep);
}

void
pari_add_defaults_module(entree *ep)
{ pari_fill_hashtable(defaults_hash, ep); }

/*********************************************************************/
/*                       PARI MAIN STACK                             */
/*********************************************************************/

#ifdef HAS_MMAP
#define PARI_STACK_ALIGN (sysconf(_SC_PAGE_SIZE))
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif
#ifndef MAP_NORESERVE
#define MAP_NORESERVE 0
#endif
static void *
pari_mainstack_malloc(size_t size)
{
  void *b;
  /* Check that the system allows reserving "size" bytes. This is just
   * a check, we immediately free the memory. */
  BLOCK_SIGINT_START;
  b = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
  BLOCK_SIGINT_END;
  if (b == MAP_FAILED) return NULL;
  BLOCK_SIGINT_START;
  munmap(b, size);

  /* Map again, this time with MAP_NORESERVE. On some operating systems
   * like Cygwin, this is needed because remapping with PROT_NONE and
   * MAP_NORESERVE does not work as expected. */
  b = mmap(NULL, size, PROT_READ|PROT_WRITE,
                       MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0);
  BLOCK_SIGINT_END;
  if (b == MAP_FAILED) return NULL;
  return b;
}

static void
pari_mainstack_mfree(void *s, size_t size)
{
  BLOCK_SIGINT_START;
  munmap(s, size);
  BLOCK_SIGINT_END;
}

/* Completely discard the memory mapped between the addresses "from"
 * and "to" (which must be page-aligned).
 *
 * We use mmap() with PROT_NONE, which means that the underlying memory
 * is freed and that the kernel should not commit memory for it. We
 * still keep the mapping such that we can change the flags to
 * PROT_READ|PROT_WRITE later.
 *
 * NOTE: remapping with MAP_FIXED and PROT_NONE is not the same as
 * calling mprotect(..., PROT_NONE) because the latter will keep the
 * memory committed (this is in particular relevant on Linux with
 * vm.overcommit = 2). This remains true even when calling
 * madvise(..., MADV_DONTNEED). */
static void
pari_mainstack_mreset(pari_sp from, pari_sp to)
{
  size_t s = to - from;
  void *addr, *res;
  if (!s) return;

  addr = (void*)from;
  BLOCK_SIGINT_START;
  res = mmap(addr, s, PROT_NONE,
             MAP_FIXED|MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0);
  BLOCK_SIGINT_END;
  if (res != addr) pari_err(e_MEM);
}

/* Commit (make available) the virtual memory mapped between the
 * addresses "from" and "to" (which must be page-aligned).
 * Return 0 if successful, -1 if failed. */
static int
pari_mainstack_mextend(pari_sp from, pari_sp to)
{
  size_t s = to - from;
  int ret;
  BLOCK_SIGINT_START;
  ret = mprotect((void*)from, s, PROT_READ|PROT_WRITE);
  BLOCK_SIGINT_END;
  return ret;
}

/* Set actual stack size to the given size. This sets st->size and
 * st->bot. If not enough system memory is available, this can fail.
 * Return 1 if successful, 0 if failed (in that case, st->size is not
 * changed) */
static int
pari_mainstack_setsize(struct pari_mainstack *st, size_t size)
{
  pari_sp newbot = st->top - size;
  /* Align newbot to pagesize */
  pari_sp alignbot = newbot & ~(pari_sp)(PARI_STACK_ALIGN - 1);
  if (pari_mainstack_mextend(alignbot, st->top))
  {
    /* Making the memory available did not work: limit vsize to the
     * current actual stack size. */
    st->vsize = st->size;
    pari_warn(warnstack, st->vsize);
    return 0;
  }
  pari_mainstack_mreset(st->vbot, alignbot);
  st->bot = newbot;
  st->size = size;
  return 1;
}

#else
#define PARI_STACK_ALIGN (0x40UL)
static void *
pari_mainstack_malloc(size_t s)
{
  char * tmp;
  BLOCK_SIGINT_START;
  tmp = malloc(s); /* NOT pari_malloc, e_MEM would be deadly */
  BLOCK_SIGINT_END;
  return tmp;
}

static void
pari_mainstack_mfree(void *s, size_t size) { (void) size; pari_free(s); }

static int
pari_mainstack_setsize(struct pari_mainstack *st, size_t size)
{
  st->bot = st->top - size;
  st->size = size;
  return 1;
}

#endif

static const size_t MIN_STACK = 500032UL;
static size_t
fix_size(size_t a)
{
  size_t ps = PARI_STACK_ALIGN;
  size_t b = a & ~(ps - 1); /* Align */
  if (b < a && b < ~(ps - 1)) b += ps;
  if (b < MIN_STACK) b = MIN_STACK;
  return b;
}

static void
pari_mainstack_alloc(int numerr, struct pari_mainstack *st, size_t rsize, size_t vsize)
{
  size_t sizemax = vsize ? vsize: rsize, s = fix_size(sizemax);
  for (;;)
  {
    st->vbot = (pari_sp)pari_mainstack_malloc(s);
    if (st->vbot) break;
    if (s == MIN_STACK) pari_err(e_MEM); /* no way out. Die */
    s = fix_size(s >> 1);
    pari_warn(numerr, s);
  }
  st->vsize = vsize ? s: 0;
  st->rsize = minuu(rsize, s);
  st->top = st->vbot+s;
  if (!pari_mainstack_setsize(st, st->rsize))
  {
    /* This should never happen since we only decrease the allocated space */
    pari_err(e_MEM);
  }
  st->memused = 0;
}

static void
pari_mainstack_free(struct pari_mainstack *st)
{
  pari_mainstack_mfree((void*)st->vbot, st->vsize ? st->vsize : fix_size(st->rsize));
  st->top = st->bot = st->vbot = 0;
  st->size = st->vsize = 0;
}

static void
pari_mainstack_resize(struct pari_mainstack *st, size_t rsize, size_t vsize)
{
  BLOCK_SIGINT_START;
  pari_mainstack_free(st);
  pari_mainstack_alloc(warnstack, st, rsize, vsize);
  BLOCK_SIGINT_END;
}

static void
pari_mainstack_use(struct pari_mainstack *st)
{
  pari_mainstack = st;
  avma = st->top; /* don't use set_avma */
}

static void
paristack_alloc(size_t rsize, size_t vsize)
{
  pari_mainstack_alloc(warnstack, pari_mainstack, rsize, vsize);
  pari_mainstack_use(pari_mainstack);
}

void
paristack_setsize(size_t rsize, size_t vsize)
{
  pari_mainstack_resize(pari_mainstack, rsize, vsize);
  pari_mainstack_use(pari_mainstack);
}

void
parivstack_resize(ulong newsize)
{
  size_t s;
  if (newsize && newsize < pari_mainstack->rsize)
    pari_err_DIM("stack sizes [parisizemax < parisize]");
  if (newsize == pari_mainstack->vsize) return;
  evalstate_reset();
  paristack_setsize(pari_mainstack->rsize, newsize);
  s = pari_mainstack->vsize ? pari_mainstack->vsize : pari_mainstack->rsize;
  if (DEBUGMEM)
    pari_warn(warner,"new maximum stack size = %lu (%.3f Mbytes)",
              s, s/1048576.);
  pari_init_errcatch();
  cb_pari_err_recover(-1);
}

void
paristack_newrsize(ulong newsize)
{
  size_t s, vsize = pari_mainstack->vsize;
  if (!newsize) newsize = pari_mainstack->rsize << 1;
  if (newsize != pari_mainstack->rsize)
    pari_mainstack_resize(pari_mainstack, newsize, vsize);
  evalstate_reset();
  s = pari_mainstack->rsize;
  if (DEBUGMEM)
    pari_warn(warner,"new stack size = %lu (%.3f Mbytes)", s, s/1048576.);
  pari_init_errcatch();
  cb_pari_err_recover(-1);
}

void
paristack_resize(ulong newsize)
{
  int64_t size = pari_mainstack->size;
  if (!newsize)
    newsize = 2 * size;
  newsize = minuu(newsize, pari_mainstack->vsize);
  if (newsize <= pari_mainstack->size) return;
  if (pari_mainstack_setsize(pari_mainstack, newsize))
  {
    if (DEBUGMEM)
      pari_warn(warner, "increasing stack size to %lu", pari_mainstack->size);
  }
  else
  {
    pari_mainstack_setsize(pari_mainstack, size);
    pari_err(e_STACK);
  }
}

void
parivstack_reset(void)
{
  pari_mainstack_setsize(pari_mainstack, pari_mainstack->rsize);
  if (avma < pari_mainstack->bot)
    pari_err_BUG("parivstack_reset [avma < bot]");
}

/* Enlarge the stack if needed such that the unused portion of the stack
 * (between bot and avma) is large enough to contain x longs. */
void
new_chunk_resize(size_t x)
{
  if (pari_mainstack->vsize==0
    || x > (avma-pari_mainstack->vbot) / sizeof(int64_t)) pari_err(e_STACK);
  while (x > (avma-pari_mainstack->bot) / sizeof(int64_t))
    paristack_resize(0);
}

/*********************************************************************/
/*                       PARI THREAD                                 */
/*********************************************************************/

/* Initial PARI thread structure t with a stack of size s and
 * argument arg */

static void
pari_thread_set_global(struct pari_global_state *gs)
{
  push_localbitprec(gs->bitprec);
  pari_set_primetab(gs->primetab);
  pari_set_seadata(gs->seadata);
  pari_set_varstate(gs->varpriority, &gs->varstate);
}

static void
pari_thread_get_global(struct pari_global_state *gs)
{
  gs->bitprec = get_localbitprec();
  gs->primetab = primetab;
  gs->seadata = pari_get_seadata();
  varstate_save(&gs->varstate);
  gs->varpriority = varpriority;
}

void
pari_thread_alloc(struct pari_thread *t, size_t s, GEN arg)
{
  pari_mainstack_alloc(warnstackthread, &t->st,s,0);
  pari_thread_get_global(&t->gs);
  t->data = arg;
}

/* Initial PARI thread structure t with a stack of size s and virtual size v
 * and argument arg */

void
pari_thread_valloc(struct pari_thread *t, size_t s, size_t v, GEN arg)
{
  pari_mainstack_alloc(warnstackthread, &t->st,s,v);
  pari_thread_get_global(&t->gs);
  t->data = arg;
}

void
pari_thread_free(struct pari_thread *t)
{
  pari_mainstack_free(&t->st);
}

void
pari_thread_init(void)
{
  int64_t var;
  pari_stackcheck_init((void*)&var);
  pari_init_blocks();
  pari_init_errcatch();
  pari_init_rand();
  pari_init_floats();
  pari_init_parser();
  pari_init_compiler();
  pari_init_evaluator();
  pari_init_files();
}

void
pari_thread_close(void)
{
  pari_thread_close_files();
  pari_close_evaluator();
  pari_close_compiler();
  pari_close_parser();
  pari_close_floats();
  pari_close_blocks();
}

GEN
pari_thread_start(struct pari_thread *t)
{
  pari_mainstack_use(&t->st);
  pari_thread_init();
  pari_thread_set_global(&t->gs);
  return t->data;
}

/*********************************************************************/
/*                       LIBPARI INIT / CLOSE                        */
/*********************************************************************/

static void
pari_exit(void)
{
  err_printf("  ***   Error in the PARI system. End of program.\n");
  exit(1);
}

static void
dflt_err_recover(int64_t errnum) { (void) errnum; pari_exit(); }

static void
dflt_pari_quit(int64_t err) { (void)err; /*do nothing*/; }

static int pari_err_display(GEN err);

/* initialize PARI data. Initialize [new|old]fun to NULL for default set. */
void
pari_init_opts(size_t parisize, ulong maxprime, ulong init_opts)
{
  ulong u;

  pari_mt_nbthreads = 0;
  cb_pari_quit = dflt_pari_quit;
  cb_pari_init_histfile = NULL;
  cb_pari_get_line_interactive = NULL;
  cb_pari_fgets_interactive = NULL;
  cb_pari_whatnow = NULL;
  cb_pari_handle_exception = NULL;
  cb_pari_err_handle = pari_err_display;
  cb_pari_pre_recover = NULL;
  cb_pari_break_loop = NULL;
  cb_pari_is_interactive = NULL;
  cb_pari_start_output = NULL;
  cb_pari_sigint = dflt_sigint_fun;
  if (init_opts&INIT_JMPm) cb_pari_err_recover = dflt_err_recover;

  pari_stackcheck_init(&u);
  pari_init_homedir();
  if (init_opts&INIT_DFTm) {
    pari_init_defaults();
    GP_DATA = default_gp_data();
    pari_init_paths();
  }

  pari_mainstack = (struct pari_mainstack *) malloc(sizeof(*pari_mainstack));
  paristack_alloc(parisize, 0);
  init_universal_constants();
  diffptr = NULL;
  if (!(init_opts&INIT_noPRIMEm))
  {
    GP_DATA->primelimit = maxprime;
    pari_init_primes(GP_DATA->primelimit);
  }
  if (!(init_opts&INIT_noINTGMPm)) pari_kernel_init();
  pari_init_graphics();
  pari_thread_init();
  pari_set_primetab(NULL);
  pari_set_seadata(NULL);
  pari_init_functions();
  pari_init_export();
  pari_var_init();
  pari_init_timer();
  pari_init_buffers();
  (void)getabstime();
  try_to_recover = 1;
  if (!(init_opts&INIT_noIMTm)) pari_mt_init();
  if ((init_opts&INIT_SIGm)) pari_sig_init(pari_sighandler);
}

void
pari_init(size_t parisize, ulong maxprime)
{ pari_init_opts(parisize, maxprime, INIT_JMPm | INIT_SIGm | INIT_DFTm); }

void
pari_close_opts(ulong init_opts)
{
  int64_t i;

  BLOCK_SIGINT_START;
  if ((init_opts&INIT_SIGm)) pari_sig_init(SIG_DFL);
  if (!(init_opts&INIT_noIMTm)) pari_mt_close();

  pari_var_close(); /* must come before destruction of functions_hash */
  for (i = 0; i < functions_tblsz; i++)
  {
    entree *ep = functions_hash[i];
    while (ep) {
      entree *EP = ep->next;
      if (!EpSTATIC(ep)) { freeep(ep); free(ep); }
      ep = EP;
    }
  }
  pari_close_mf();
  pari_thread_close();
  pari_close_export();
  pari_close_files();
  pari_close_homedir();
  if (!(init_opts&INIT_noINTGMPm)) pari_kernel_close();

  free((void*)functions_hash);
  free((void*)defaults_hash);
  if (diffptr) pari_close_primes();
  free(current_logfile);
  free(current_psfile);
  pari_mainstack_free(pari_mainstack);
  free((void*)pari_mainstack);
  pari_stack_delete(&s_MODULES);
  if (pari_datadir) free(pari_datadir);
  if (init_opts&INIT_DFTm)
  { /* delete GP_DATA */
    pari_close_paths();
    if (GP_DATA->hist->v) free((void*)GP_DATA->hist->v);
    if (GP_DATA->pp->cmd) free((void*)GP_DATA->pp->cmd);
    if (GP_DATA->help) free((void*)GP_DATA->help);
    if (GP_DATA->plothsizes) free((void*)GP_DATA->plothsizes);
    if (GP_DATA->colormap) pari_free(GP_DATA->colormap);
    if (GP_DATA->graphcolors) pari_free(GP_DATA->graphcolors);
    free((void*)GP_DATA->prompt);
    free((void*)GP_DATA->prompt_cont);
    free((void*)GP_DATA->histfile);
  }
  BLOCK_SIGINT_END;
}

void
pari_close(void)
{ pari_close_opts(INIT_JMPm | INIT_SIGm | INIT_DFTm); }

/*******************************************************************/
/*                                                                 */
/*                         ERROR RECOVERY                          */
/*                                                                 */
/*******************************************************************/
void
gp_context_save(struct gp_context* rec)
{
  rec->prettyp = GP_DATA->fmt->prettyp;
  rec->listloc = next_block;
  rec->iferr_env = iferr_env;
  rec->err_data  = global_err_data;
  varstate_save(&rec->var);
  evalstate_save(&rec->eval);
  parsestate_save(&rec->parse);
  filestate_save(&rec->file);
}

void
gp_context_restore(struct gp_context* rec)
{
  int64_t i;

  if (!try_to_recover) return;
  /* disable gp_context_restore() and SIGINT */
  try_to_recover = 0;
  BLOCK_SIGINT_START
  if (DEBUGMEM>2) err_printf("entering recover(), loc = %ld\n", rec->listloc);
  evalstate_restore(&rec->eval);
  parsestate_restore(&rec->parse);
  filestate_restore(&rec->file);
  global_err_data = rec->err_data;
  iferr_env = rec->iferr_env;
  GP_DATA->fmt->prettyp = rec->prettyp;

  for (i = 0; i < functions_tblsz; i++)
  {
    entree *ep = functions_hash[i];
    while (ep)
    {
      entree *EP = ep->next;
      switch(EpVALENCE(ep))
      {
        case EpVAR:
          while (pop_val_if_newer(ep,rec->listloc)) /* empty */;
          break;
        case EpNEW: break;
      }
      ep = EP;
    }
  }
  varstate_restore(&rec->var);
  if (DEBUGMEM>2) err_printf("leaving recover()\n");
  BLOCK_SIGINT_END
  try_to_recover = 1;
}

static void
err_recover(int64_t numerr)
{
  if (cb_pari_pre_recover)
    cb_pari_pre_recover(numerr);
  evalstate_reset();
  killallfiles();
  pari_init_errcatch();
  cb_pari_err_recover(numerr);
}

static void
err_init(void)
{
  /* make sure pari_err msg starts at the beginning of line */
  if (!pari_last_was_newline()) pari_putc('\n');
  pariOut->flush();
  pariErr->flush();
  out_term_color(pariErr, c_ERR);
}

static void
err_init_msg(int user)
{
  const char *gp_function_name;
  out_puts(pariErr, "  *** ");
  if (!user && (gp_function_name = closure_func_err()))
    out_printf(pariErr, "%s: ", gp_function_name);
  else
    out_puts(pariErr, "  ");
}

void
pari_warn(int numerr, ...)
{
  char *ch1;
  va_list ap;

  va_start(ap,numerr);

  err_init();
  err_init_msg(numerr==warnuser || numerr==warnstack);
  switch (numerr)
  {
    case warnuser:
      out_puts(pariErr, "user warning: ");
      out_print0(pariErr, NULL, va_arg(ap, GEN), f_RAW);
      break;

    case warnmem:
      out_puts(pariErr, "collecting garbage in "); ch1=va_arg(ap, char*);
      out_vprintf(pariErr, ch1,ap); out_putc(pariErr, '.');
      break;

    case warner:
      out_puts(pariErr, "Warning: "); ch1=va_arg(ap, char*);
      out_vprintf(pariErr, ch1,ap); out_putc(pariErr, '.');
      break;

    case warnprec:
      out_vprintf(pariErr, "Warning: increasing prec in %s; new prec = %ld",
                      ap);
      break;

    case warnfile:
      out_puts(pariErr, "Warning: failed to "),
      ch1 = va_arg(ap, char*);
      out_printf(pariErr, "%s: %s", ch1, va_arg(ap, char*));
      break;

    case warnstack:
    case warnstackthread:
    {
      ulong  s = va_arg(ap, ulong);
      char buf[128];
      const char * stk = numerr == warnstackthread
                         || mt_is_thread() ? "thread": "PARI";
      sprintf_s(buf, sizeof(buf), "Warning: not enough memory, new %s stack %llu", stk, s);
      out_puts(pariErr,buf);
      break;
    }
  }
  va_end(ap);
  out_term_color(pariErr, c_NONE);
  out_putc(pariErr, '\n');
  pariErr->flush();
}
void
pari_sigint(const char *time_s)
{
  int recover=0;
  BLOCK_SIGALRM_START
  err_init();
  closure_err(0);
  err_init_msg(0);
  out_puts(pariErr, "user interrupt after ");
  out_puts(pariErr, time_s);
  out_term_color(pariErr, c_NONE);
  pariErr->flush();
  if (cb_pari_handle_exception)
    recover = cb_pari_handle_exception(-1);
  if (!recover && !block)
    PARI_SIGINT_pending = 0;
  BLOCK_SIGINT_END
  if (!recover) err_recover(e_MISC);
}

#define retmkerr2(x,y)\
  do { GEN _v = cgetg(3, t_ERROR);\
       _v[1] = (x);\
       gel(_v,2) = (y); return _v; } while(0)
#define retmkerr3(x,y,z)\
  do { GEN _v = cgetg(4, t_ERROR);\
       _v[1] = (x);\
       gel(_v,2) = (y);\
       gel(_v,3) = (z); return _v; } while(0)
#define retmkerr4(x,y,z,t)\
  do { GEN _v = cgetg(5, t_ERROR);\
       _v[1] = (x);\
       gel(_v,2) = (y);\
       gel(_v,3) = (z);\
       gel(_v,4) = (t); return _v; } while(0)
#define retmkerr5(x,y,z,t,u)\
  do { GEN _v = cgetg(6, t_ERROR);\
       _v[1] = (x);\
       gel(_v,2) = (y);\
       gel(_v,3) = (z);\
       gel(_v,4) = (t);\
       gel(_v,5) = (u); return _v; } while(0)
#define retmkerr6(x,y,z,t,u,v)\
  do { GEN _v = cgetg(7, t_ERROR);\
       _v[1] = (x);\
       gel(_v,2) = (y);\
       gel(_v,3) = (z);\
       gel(_v,4) = (t);\
       gel(_v,5) = (u);\
       gel(_v,6) = (v); return _v; } while(0)

static GEN
pari_err2GEN(int64_t numerr, va_list ap)
{
  switch ((enum err_list) numerr)
  {
  case e_SYNTAX:
    {
      const char *msg = va_arg(ap, char*);
      const char *s = va_arg(ap,char *);
      const char *entry = va_arg(ap,char *);
      retmkerr3(numerr,strtoGENstr(msg), mkvecsmall2((int64_t)s,(int64_t)entry));
    }
  case e_MISC: case e_ALARM:
    {
      const char *ch1 = va_arg(ap, char*);
      retmkerr2(numerr, gvsprintf(ch1,ap));
    }
  case e_NOTFUNC:
  case e_USER:
    retmkerr2(numerr,va_arg(ap, GEN));
  case e_FILE:
    {
      const char *f = va_arg(ap, const char*);
      retmkerr3(numerr, strtoGENstr(f), strtoGENstr(va_arg(ap, char*)));
    }
  case e_FILEDESC:
    {
      const char *f = va_arg(ap, const char*);
      retmkerr3(numerr, strtoGENstr(f), stoi(va_arg(ap, int64_t)));
    }
  case e_OVERFLOW:
  case e_IMPL:
  case e_DIM:
  case e_CONSTPOL:
  case e_ROOTS0:
  case e_FLAG:
  case e_PREC:
  case e_BUG:
  case e_ARCH:
  case e_PACKAGE:
    retmkerr2(numerr, strtoGENstr(va_arg(ap, char*)));
  case e_MODULUS:
  case e_VAR:
    {
      const char *f = va_arg(ap, const char*);
      GEN x = va_arg(ap, GEN);
      GEN y = va_arg(ap, GEN);
      retmkerr4(numerr, strtoGENstr(f), x,y);
    }
  case e_INV:
  case e_IRREDPOL:
  case e_PRIME:
  case e_SQRTN:
  case e_TYPE:
    {
      const char *f = va_arg(ap, const char*);
      GEN x = va_arg(ap, GEN);
      retmkerr3(numerr, strtoGENstr(f), x);
    }
  case e_COPRIME: case e_OP: case e_TYPE2:
    {
      const char *f = va_arg(ap, const char*);
      GEN x = va_arg(ap, GEN);
      GEN y = va_arg(ap, GEN);
      retmkerr4(numerr,strtoGENstr(f),x,y);
    }
  case e_COMPONENT:
    {
      const char *f= va_arg(ap, const char *);
      const char *op = va_arg(ap, const char *);
      GEN l = va_arg(ap, GEN);
      GEN x = va_arg(ap, GEN);
      retmkerr5(numerr,strtoGENstr(f),strtoGENstr(op),l,x);
    }
  case e_DOMAIN:
    {
      const char *f = va_arg(ap, const char*);
      const char *v = va_arg(ap, const char *);
      const char *op = va_arg(ap, const char *);
      GEN l = va_arg(ap, GEN);
      GEN x = va_arg(ap, GEN);
      retmkerr6(numerr,strtoGENstr(f),strtoGENstr(v),strtoGENstr(op),l,x);
    }
  case e_PRIORITY:
    {
      const char *f = va_arg(ap, const char*);
      GEN x = va_arg(ap, GEN);
      const char *op = va_arg(ap, const char *);
      int64_t v = va_arg(ap, int64_t);
      retmkerr5(numerr,strtoGENstr(f),x,strtoGENstr(op),stoi(v));
    }
  case e_MAXPRIME:
    retmkerr2(numerr, utoi(va_arg(ap, ulong)));
  case e_STACK:
    return err_e_STACK;
  case e_STACKTHREAD:
    retmkerr3(numerr, utoi(va_arg(ap, ulong)), utoi(va_arg(ap, ulong)));
  default:
    return mkerr(numerr);
  }
}

static char *
type_dim(GEN x)
{
  char *v = stack_malloc(64);
  switch(typ(x))
  {
    case t_MAT:
    {
      int64_t l = lg(x), r = (l == 1)? 1: lgcols(x);
      sprintf_s(v, 64, "t_MAT (%lldx%lld)", r-1,l-1);
      break;
    }
    case t_COL:
      sprintf_s(v, 64, "t_COL (%lld elts)", lg(x)-1);
      break;
    case t_VEC:
      sprintf_s(v, 64, "t_VEC (%lld elts)", lg(x)-1);
      break;
    default:
      v = (char*)type_name(typ(x));
  }
  return v;
}

static char *
gdisplay(GEN x)
{
  char *s = GENtostr_raw(x);
  if (strlen(s) < 1600) return s;
  if (! GP_DATA->breakloop) return (char*)"(...)";
  return stack_sprintf("\n  ***  (...) Huge %s omitted; you can access it via dbg_err()", type_name(typ(x)));
}

char *
pari_err2str(GEN e)
{
  int64_t numerr = err_get_num(e);
  switch ((enum err_list) numerr)
  {
  case e_ALARM:
    return pari_sprintf("alarm interrupt after %Ps.",gel(e,2));
  case e_MISC:
    return pari_sprintf("%Ps.",gel(e,2));

  case e_ARCH:
    return pari_sprintf("sorry, '%Ps' not available on this system.",gel(e,2));
  case e_BUG:
    return pari_sprintf("bug in %Ps, please report.",gel(e,2));
  case e_CONSTPOL:
    return pari_sprintf("constant polynomial in %Ps.", gel(e,2));
  case e_COPRIME:
    return pari_sprintf("elements not coprime in %Ps:\n    %s\n    %s",
                        gel(e,2), gdisplay(gel(e,3)), gdisplay(gel(e,4)));
  case e_DIM:
    return pari_sprintf("inconsistent dimensions in %Ps.", gel(e,2));
  case e_FILE:
    return pari_sprintf("error opening %Ps: `%Ps'.", gel(e,2), gel(e,3));
  case e_FILEDESC:
    return pari_sprintf("invalid file descriptor in %Ps [%Ps]", gel(e,2), gel(e,3));
  case e_FLAG:
    return pari_sprintf("invalid flag in %Ps.", gel(e,2));
  case e_IMPL:
    return pari_sprintf("sorry, %Ps is not yet implemented.", gel(e,2));
  case e_PACKAGE:
    return pari_sprintf("package %Ps is required, please install it.", gel(e,2));
  case e_INV:
    return pari_sprintf("impossible inverse in %Ps: %s.", gel(e,2),
                        gdisplay(gel(e,3)));
  case e_IRREDPOL:
    return pari_sprintf("not an irreducible polynomial in %Ps: %s.",
                        gel(e,2), gdisplay(gel(e,3)));
  case e_MAXPRIME:
    {
      const char * msg = "not enough precomputed primes";
      ulong c = itou(gel(e,2));
      if (c) return pari_sprintf("%s, need primelimit ~ %lu.",msg, c);
      else   return pari_strdup(msg);
    }
  case e_MEM:
    return pari_strdup("not enough memory");
  case e_MODULUS:
    {
      GEN x = gel(e,3), y = gel(e,4);
      return pari_sprintf("inconsistent moduli in %Ps: %s != %s",
                          gel(e,2), gdisplay(x), gdisplay(y));
    }
  case e_NONE: return NULL;
  case e_NOTFUNC:
    return pari_strdup("not a function in function call");
  case e_OP: case e_TYPE2:
    {
      pari_sp av = avma;
      char *v;
      const char *f, *op = GSTR(gel(e,2));
      const char *what = numerr == e_OP? "inconsistent": "forbidden";
      GEN x = gel(e,3);
      GEN y = gel(e,4);
      switch(*op)
      {
      case '+': f = "addition"; break;
      case '*': f = "multiplication"; break;
      case '/': case '%': case '\\': f = "division"; break;
      case '=': op = "-->"; f = "assignment"; break;
      default:  f = op; op = ","; break;
      }
      v = pari_sprintf("%s %s %s %s %s.", what,f,type_dim(x),op,type_dim(y));
      set_avma(av); return v;
    }
  case e_COMPONENT:
    {
      const char *f= GSTR(gel(e,2));
      const char *op= GSTR(gel(e,3));
      GEN l = gel(e,4);
      if (!*f)
        return pari_sprintf("nonexistent component: index %s %Ps",op,l);
      return pari_sprintf("nonexistent component in %s: index %s %Ps",f,op,l);
    }
  case e_DOMAIN:
    {
      const char *f = GSTR(gel(e,2));
      const char *v = GSTR(gel(e,3));
      const char *op= GSTR(gel(e,4));
      GEN l = gel(e,5);
      if (!*op)
        return pari_sprintf("domain error in %s: %s out of range",f,v);
      return pari_sprintf("domain error in %s: %s %s %Ps",f,v,op,l);
    }
  case e_PRIORITY:
    {
      const char *f = GSTR(gel(e,2));
      int64_t vx = gvar(gel(e,3));
      const char *op= GSTR(gel(e,4));
      int64_t v = itos(gel(e,5));
      return pari_sprintf("incorrect priority in %s: variable %Ps %s %Ps",f,
             pol_x(vx), op, pol_x(v));
    }
  case e_OVERFLOW:
    return pari_sprintf("overflow in %Ps.", gel(e,2));
  case e_PREC:
    return pari_sprintf("precision too low in %Ps.", gel(e,2));
  case e_PRIME:
    return pari_sprintf("not a prime number in %Ps: %s.",
                        gel(e,2), gdisplay(gel(e,3)));
  case e_ROOTS0:
    return pari_sprintf("zero polynomial in %Ps.", gel(e,2));
  case e_SQRTN:
    return pari_sprintf("not an n-th power residue in %Ps: %s.",
                        gel(e,2), gdisplay(gel(e,3)));
  case e_STACK:
  case e_STACKTHREAD:
    {
      const char *stack = numerr == e_STACK? "PARI": "thread";
      const char *var = numerr == e_STACK? "parisizemax": "threadsizemax";
      size_t rsize = numerr == e_STACKTHREAD && GP_DATA->threadsize ?
                                GP_DATA->threadsize: pari_mainstack->rsize;
      size_t vsize = numerr == e_STACK? pari_mainstack->vsize:
                                        GP_DATA->threadsizemax;
      char *buf = (char *) pari_malloc(512*sizeof(char));
      if (vsize)
      {
        sprintf_s(buf, 512*sizeof(char), "the %s stack overflows !\n"
            "  current stack size: %llu (%.3f Mbytes)\n"
            "  [hint] you can increase '%s' using default()\n",
            stack, (ulong)vsize, (double)vsize/1048576., var);
      }
      else
      {
        sprintf_s(buf, 512 * sizeof(char), "the %s stack overflows !\n"
            "  current stack size: %llu (%.3f Mbytes)\n"
            "  [hint] set '%s' to a nonzero value in your GPRC\n",
            stack, (ulong)rsize, (double)rsize/1048576., var);
      }
      return buf;
    }
  case e_SYNTAX:
    return pari_strdup(GSTR(gel(e,2)));
  case e_TYPE:
    return pari_sprintf("incorrect type in %Ps (%s).",
                        gel(e,2), type_name(typ(gel(e,3))));
  case e_USER:
    return pari_sprint0("user error: ", gel(e,2), f_RAW);
  case e_VAR:
    {
      GEN x = gel(e,3), y = gel(e,4);
      return pari_sprintf("inconsistent variables in %Ps, %Ps != %Ps.",
                          gel(e,2), pol_x(varn(x)), pol_x(varn(y)));
    }
  }
  return NULL; /*LCOV_EXCL_LINE*/
}

static int
pari_err_display(GEN err)
{
  int64_t numerr=err_get_num(err);
  err_init();
  if (numerr==e_SYNTAX)
  {
    const char *msg = GSTR(gel(err,2));
    const char *s     = (const char *) gmael(err,3,1);
    const char *entry = (const char *) gmael(err,3,2);
    print_errcontext(pariErr, msg, s, entry);
  }
  else
  {
    char *s;
    closure_err(0);
    err_init_msg(numerr==e_USER);
    s = pari_err2str(err); pariErr->puts(s); pari_free(s);
    if (numerr==e_NOTFUNC)
    {
      GEN fun = gel(err,2);
      if (gequalX(fun))
      {
        entree *ep = varentries[varn(fun)];
        const char *t = ep->name;
        if (cb_pari_whatnow) cb_pari_whatnow(pariErr,t,1);
      }
    }
  }
  out_term_color(pariErr, c_NONE);
  pariErr->flush(); return 0;
}

void
pari_err(int numerr, ...)
{
  va_list ap;
  GEN E;

  va_start(ap,numerr);

  if (numerr)
    E = pari_err2GEN(numerr,ap);
  else
  {
    E = va_arg(ap,GEN);
    numerr = err_get_num(E);
  }
  global_err_data = E;
  if (*iferr_env) longjmp(*iferr_env, numerr);
  mt_err_recover(numerr);
  va_end(ap);
  if (cb_pari_err_handle &&
      cb_pari_err_handle(E)) return;
  if (cb_pari_handle_exception &&
      cb_pari_handle_exception(numerr)) return;
  err_recover(numerr);
}

GEN
pari_err_last(void) { return global_err_data; }

const char *
numerr_name(int64_t numerr)
{
  switch ((enum err_list) numerr)
  {
  case e_ALARM:    return "e_ALARM";
  case e_ARCH:     return "e_ARCH";
  case e_BUG:      return "e_BUG";
  case e_COMPONENT: return "e_COMPONENT";
  case e_CONSTPOL: return "e_CONSTPOL";
  case e_COPRIME:  return "e_COPRIME";
  case e_DIM:      return "e_DIM";
  case e_DOMAIN:   return "e_DOMAIN";
  case e_FILE:     return "e_FILE";
  case e_FILEDESC: return "e_FILEDESC";
  case e_FLAG:     return "e_FLAG";
  case e_IMPL:     return "e_IMPL";
  case e_INV:      return "e_INV";
  case e_IRREDPOL: return "e_IRREDPOL";
  case e_MAXPRIME: return "e_MAXPRIME";
  case e_MEM:      return "e_MEM";
  case e_MISC:     return "e_MISC";
  case e_MODULUS:  return "e_MODULUS";
  case e_NONE:     return "e_NONE";
  case e_NOTFUNC:  return "e_NOTFUNC";
  case e_OP:       return "e_OP";
  case e_OVERFLOW: return "e_OVERFLOW";
  case e_PACKAGE:  return "e_PACKAGE";
  case e_PREC:     return "e_PREC";
  case e_PRIME:    return "e_PRIME";
  case e_PRIORITY: return "e_PRIORITY";
  case e_ROOTS0:   return "e_ROOTS0";
  case e_SQRTN:    return "e_SQRTN";
  case e_STACK:    return "e_STACK";
  case e_SYNTAX:   return "e_SYNTAX";
  case e_STACKTHREAD:   return "e_STACKTHREAD";
  case e_TYPE2:    return "e_TYPE2";
  case e_TYPE:     return "e_TYPE";
  case e_USER:     return "e_USER";
  case e_VAR:      return "e_VAR";
  }
  return "invalid error number";
}

int64_t
name_numerr(const char *s)
{
  if (!strcmp(s,"e_ALARM"))    return e_ALARM;
  if (!strcmp(s,"e_ARCH"))     return e_ARCH;
  if (!strcmp(s,"e_BUG"))      return e_BUG;
  if (!strcmp(s,"e_COMPONENT")) return e_COMPONENT;
  if (!strcmp(s,"e_CONSTPOL")) return e_CONSTPOL;
  if (!strcmp(s,"e_COPRIME"))  return e_COPRIME;
  if (!strcmp(s,"e_DIM"))      return e_DIM;
  if (!strcmp(s,"e_DOMAIN"))   return e_DOMAIN;
  if (!strcmp(s,"e_FILE"))     return e_FILE;
  if (!strcmp(s,"e_FILEDESC")) return e_FILEDESC;
  if (!strcmp(s,"e_FLAG"))     return e_FLAG;
  if (!strcmp(s,"e_IMPL"))     return e_IMPL;
  if (!strcmp(s,"e_INV"))      return e_INV;
  if (!strcmp(s,"e_IRREDPOL")) return e_IRREDPOL;
  if (!strcmp(s,"e_MAXPRIME")) return e_MAXPRIME;
  if (!strcmp(s,"e_MEM"))      return e_MEM;
  if (!strcmp(s,"e_MISC"))     return e_MISC;
  if (!strcmp(s,"e_MODULUS"))  return e_MODULUS;
  if (!strcmp(s,"e_NONE"))     return e_NONE;
  if (!strcmp(s,"e_NOTFUNC"))  return e_NOTFUNC;
  if (!strcmp(s,"e_OP"))       return e_OP;
  if (!strcmp(s,"e_OVERFLOW")) return e_OVERFLOW;
  if (!strcmp(s,"e_PACKAGE"))  return e_PACKAGE;
  if (!strcmp(s,"e_PREC"))     return e_PREC;
  if (!strcmp(s,"e_PRIME"))    return e_PRIME;
  if (!strcmp(s,"e_PRIORITY")) return e_PRIORITY;
  if (!strcmp(s,"e_ROOTS0"))   return e_ROOTS0;
  if (!strcmp(s,"e_SQRTN"))    return e_SQRTN;
  if (!strcmp(s,"e_STACK"))    return e_STACK;
  if (!strcmp(s,"e_SYNTAX"))   return e_SYNTAX;
  if (!strcmp(s,"e_TYPE"))     return e_TYPE;
  if (!strcmp(s,"e_TYPE2"))    return e_TYPE2;
  if (!strcmp(s,"e_USER"))     return e_USER;
  if (!strcmp(s,"e_VAR"))      return e_VAR;
  pari_err(e_MISC,"unknown error name");
  return -1; /* LCOV_EXCL_LINE */
}

GEN
errname(GEN err)
{
  if (typ(err)!=t_ERROR) pari_err_TYPE("errname",err);
  return strtoGENstr(numerr_name(err_get_num(err)));
}

/* Try f (trapping error e), recover using r (break_loop, if NULL) */
GEN
trap0(const char *e, GEN r, GEN f)
{
  int64_t numerr = CATCH_ALL;
  GEN x;
  if (!e || !*e) numerr = CATCH_ALL;
  else numerr = name_numerr(e);
  if (!f) {
    pari_warn(warner,"default handlers are no longer supported --> ignored");
    return gnil;
  }
  x = closure_trapgen(f, numerr);
  if (x == (GEN)1L) x = r? closure_evalgen(r): gnil;
  return x;
}

/*******************************************************************/
/*                                                                */
/*                       CLONING & COPY                            */
/*                  Replicate an existing GEN                      */
/*                                                                 */
/*******************************************************************/
/* lontyp[tx] = 0 (non recursive type) or number of codewords for type tx */
const  int64_t lontyp[] = { 0,0,0,1,1,2,1,2,1,1, 2,2,0,1,1,1,1,1,1,1, 2,0,0,2,2,1 };

static GEN
list_internal_copy(GEN z, int64_t nmax)
{
  int64_t i, l;
  GEN a;
  if (!z) return NULL;
  l = lg(z);
  a = newblock(nmax+1);
  for (i = 1; i < l; i++) gel(a,i) = gel(z,i)? gclone(gel(z,i)): gen_0;
  a[0] = z[0]; setisclone(a); return a;
}

static void
listassign(GEN x, GEN y)
{
  int64_t nmax = list_nmax(x);
  GEN L = list_data(x);
  if (!nmax && L) nmax = lg(L) + 32; /* not malloc'ed yet */
  y[1] = evaltyp(list_typ(x))|evallg(nmax);
  list_data(y) = list_internal_copy(L, nmax);
}

/* transform a non-malloced list (e.g. from gtolist or gtomap) to a malloced
 * list suitable for listput */
GEN
listinit(GEN x)
{
  GEN y = cgetg(3, t_LIST);
  listassign(x, y); return y;
}

/* copy list on the PARI stack */
GEN
listcopy(GEN x)
{
  GEN y = mklist(), L = list_data(x);
  if (L) list_data(y) = gcopy(L);
  y[1] = evaltyp(list_typ(x));
  return y;
}

GEN
gcopy(GEN x)
{
  int64_t tx = typ(x), lx, i;
  GEN y;
  switch(tx)
  { /* non recursive types */
    case t_INT: return signe(x)? icopy(x): gen_0;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return leafcopy(x);
    /* one more special case */
    case t_LIST: return listcopy(x);
  }
  y = cgetg_copy(x, &lx);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy(gel(x,i));
  return y;
}

/* as gcopy, but truncate to the first lx components if recursive type
 * [ leaves use their own lg ]. No checks. */
GEN
gcopy_lg(GEN x, int64_t lx)
{
  int64_t tx = typ(x), i;
  GEN y;
  switch(tx)
  { /* non recursive types */
    case t_INT: return signe(x)? icopy(x): gen_0;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return leafcopy(x);
    /* one more special case */
    case t_LIST: return listcopy(x);
  }
  y = cgetg(lx, tx);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy(gel(x,i));
  return y;
}

/* cf cgetg_copy: "allocate" (by updating first codeword only) for subsequent
 * copy of x, as if avma = *AVMA */
INLINE GEN
cgetg_copy_avma(GEN x, int64_t *plx, pari_sp *AVMA) {
  GEN z;
  *plx = lg(x);
  z = ((GEN)*AVMA) - *plx;
  z[0] = x[0] & (TYPBITS|LGBITS);
  *AVMA = (pari_sp)z; return z;
}
INLINE GEN
cgetlist_avma(pari_sp *AVMA)
{
  GEN y = ((GEN)*AVMA) - 3;
  y[0] = _evallg(3) | evaltyp(t_LIST);
  *AVMA = (pari_sp)y; return y;
}

/* copy x as if avma = *AVMA, update *AVMA */
GEN
gcopy_avma(GEN x, pari_sp *AVMA)
{
  int64_t i, lx, tx = typ(x);
  GEN y;

  switch(typ(x))
  { /* non recursive types */
    case t_INT:
      if (lgefint(x) == 2) return gen_0;
      *AVMA = (pari_sp)icopy_avma(x, *AVMA);
      return (GEN)*AVMA;
    case t_REAL: case t_STR: case t_VECSMALL:
      *AVMA = (pari_sp)leafcopy_avma(x, *AVMA);
      return (GEN)*AVMA;

    /* one more special case */
    case t_LIST:
      y = cgetlist_avma(AVMA);
      listassign(x, y); return y;

  }
  y = cgetg_copy_avma(x, &lx, AVMA);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy_avma(gel(x,i), AVMA);
  return y;
}

/* [copy_bin/bin_copy:] same as gcopy_avma but use NULL to code an exact 0, and
 * make shallow copies of finalized t_LISTs */
static GEN
gcopy_av0(GEN x, pari_sp *AVMA)
{
  int64_t i, lx, tx = typ(x);
  GEN y;

  switch(tx)
  { /* non recursive types */
    case t_INT:
      if (!signe(x)) return NULL; /* special marker */
      *AVMA = (pari_sp)icopy_avma(x, *AVMA);
      return (GEN)*AVMA;
    case t_LIST:
      if (list_data(x) && !list_nmax(x)) break; /* not finalized, need copy */
      /* else finalized: shallow copy */
    case t_REAL: case t_STR: case t_VECSMALL:
      *AVMA = (pari_sp)leafcopy_avma(x, *AVMA);
      return (GEN)*AVMA;
  }
  y = cgetg_copy_avma(x, &lx, AVMA);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy_av0(gel(x,i), AVMA);
  return y;
}

INLINE GEN
icopy_avma_canon(GEN x, pari_sp AVMA)
{
  int64_t i, lx = lgefint(x);
  GEN y = ((GEN)AVMA) - lx;
  y[0] = evaltyp(t_INT)|evallg(lx); /* kills isclone */
  y[1] = x[1]; x = int_MSW(x);
  for (i=2; i<lx; i++, x = int_precW(x)) y[i] = *x;
  return y;
}

/* [copy_bin_canon:] same as gcopy_av0, but copy integers in
 * canonical (native kernel) form and make a full copy of t_LISTs */
static GEN
gcopy_av0_canon(GEN x, pari_sp *AVMA)
{
  int64_t i, lx, tx = typ(x);
  GEN y;

  switch(tx)
  { /* non recursive types */
    case t_INT:
      if (!signe(x)) return NULL; /* special marker */
      *AVMA = (pari_sp)icopy_avma_canon(x, *AVMA);
      return (GEN)*AVMA;
    case t_REAL: case t_STR: case t_VECSMALL:
      *AVMA = (pari_sp)leafcopy_avma(x, *AVMA);
      return (GEN)*AVMA;

    /* one more special case */
    case t_LIST:
    {
      int64_t t = list_typ(x);
      GEN y = cgetlist_avma(AVMA), z = list_data(x);
      if (z) {
        list_data(y) = gcopy_av0_canon(z, AVMA);
        y[1] = evaltyp(t)|evallg(lg(z)-1);
      } else {
        list_data(y) = NULL;
        y[1] = evaltyp(t);
      }
      return y;
    }
  }
  y = cgetg_copy_avma(x, &lx, AVMA);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy_av0_canon(gel(x,i), AVMA);
  return y;
}

/* [copy_bin/bin_copy:] size (number of words) required for
 * gcopy_av0_canon(x) */
static int64_t
taille0_canon(GEN x)
{
  int64_t i,n,lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT: return signe(x)? lgefint(x): 0;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return lg(x);

    /* one more special case */
    case t_LIST:
    {
      GEN L = list_data(x);
      return L? 3 + taille0_canon(L): 3;
    }
  }
  n = lx = lg(x);
  for (i=lontyp[tx]; i<lx; i++) n += taille0_canon(gel(x,i));
  return n;
}

/* [copy_bin/bin_copy:] size (number of words) required for gcopy_av0(x) */
static int64_t
taille0(GEN x)
{
  int64_t i,n,lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT:
      lx = lgefint(x);
      return lx == 2? 0: lx;
    case t_LIST:
    {
      GEN L = list_data(x);
      if (L && !list_nmax(x)) break; /* not finalized, deep copy */
    }
    /* else finalized: shallow */
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      return lg(x);
  }
  n = lx = lg(x);
  for (i=lontyp[tx]; i<lx; i++) n += taille0(gel(x,i));
  return n;
}

static int64_t
gsizeclone_i(GEN x)
{
  int64_t i,n,lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT: lx = lgefint(x); return lx == 2? 0: lx;;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return lg(x);

    case t_LIST: return 3;
    default:
      n = lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) n += gsizeclone_i(gel(x,i));
      return n;
  }
}

/* #words needed to clone x; t_LIST is a special case since list_data() is
 * malloc'ed later, in list_internal_copy() */
static int64_t
gsizeclone(GEN x) { return (typ(x) == t_INT)? lgefint(x): gsizeclone_i(x); }

int64_t
gsizeword(GEN x)
{
  int64_t i, n, lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT:
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return lg(x);

    case t_LIST:
      x = list_data(x);
      return x? 3 + gsizeword(x): 3;

    default:
      n = lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) n += gsizeword(gel(x,i));
      return n;
  }
}
int64_t
gsizebyte(GEN x) { return gsizeword(x) * sizeof(int64_t); }

/* return a clone of x structured as a gcopy */
GENbin*
copy_bin(GEN x)
{
  int64_t t = taille0(x);
  GENbin *p = (GENbin*)pari_malloc(sizeof(GENbin) + t*sizeof(int64_t));
  pari_sp AVMA = (pari_sp)(GENbinbase(p) + t);
  p->rebase = &shiftaddress;
  p->len = t;
  p->x   = gcopy_av0(x, &AVMA);
  p->base= (GEN)AVMA; return p;
}

/* same, writing t_INT in canonical native form */
GENbin*
copy_bin_canon(GEN x)
{
  int64_t t = taille0_canon(x);
  GENbin *p = (GENbin*)pari_malloc(sizeof(GENbin) + t*sizeof(int64_t));
  pari_sp AVMA = (pari_sp)(GENbinbase(p) + t);
  p->rebase = &shiftaddress_canon;
  p->len = t;
  p->x   = gcopy_av0_canon(x, &AVMA);
  p->base= (GEN)AVMA; return p;
}

GEN
gclone(GEN x)
{
  int64_t i,lx,tx = typ(x), t = gsizeclone(x);
  GEN y = newblock(t);
  switch(tx)
  { /* non recursive types */
    case t_INT:
      lx = lgefint(x);
      y[0] = evaltyp(t_INT)|evallg(lx);
      for (i=1; i<lx; i++) y[i] = x[i];
      break;
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      lx = lg(x);
      for (i=0; i<lx; i++) y[i] = x[i];
      break;

    /* one more special case */
    case t_LIST:
      y[0] = evaltyp(t_LIST)|_evallg(3);
      listassign(x, y);
      break;
    default: {
      pari_sp AVMA = (pari_sp)(y + t);
      lx = lg(x);
      y[0] = x[0];
      if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
      for (; i<lx; i++) gel(y,i) = gcopy_avma(gel(x,i), &AVMA);
    }
  }
  setisclone(y); return y;
}

void
shiftaddress(GEN x, int64_t dec)
{
  int64_t i, lx, tx = typ(x);
  if (is_recursive_t(tx))
  {
    if (tx == t_LIST)
    {
      if (!list_data(x) || list_nmax(x)) return; /* empty or finalized */
      /* not finalized, update pointers  */
    }
    lx = lg(x);
    for (i=lontyp[tx]; i<lx; i++) {
      if (!x[i]) gel(x,i) = gen_0;
      else
      {
        x[i] += dec;
        shiftaddress(gel(x,i), dec);
      }
    }
  }
}

void
shiftaddress_canon(GEN x, int64_t dec)
{
  int64_t i, lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT: {
      GEN y;
      lx = lgefint(x); if (lx <= 3) return;
      y = x + 2;
      x = int_MSW(x);  if (x == y) return;
      while (x > y) { lswap(*x, *y); x = int_precW(x); y++; }
      break;
    }
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      break;

    /* one more special case */
    case t_LIST:
      if (!list_data(x)) break;
    default: /* Fall through */
      lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) {
        if (!x[i]) gel(x,i) = gen_0;
        else
        {
          x[i] += dec;
          shiftaddress_canon(gel(x,i), dec);
        }
      }
  }
}

/********************************************************************/
/**                                                                **/
/**                INSERT DYNAMIC OBJECT IN STRUCTURE              **/
/**                                                                **/
/********************************************************************/
GEN
obj_reinit(GEN S)
{
  GEN s, T = leafcopy(S);
  int64_t a = lg(T)-1;
  s = gel(T,a);
  gel(T,a) = zerovec(lg(s)-1);
  return T;
}

GEN
obj_init(int64_t d, int64_t n)
{
  GEN S = cgetg(d+2, t_VEC);
  gel(S, d+1) = zerovec(n);
  return S;
}
/* insert O in S [last position] at position K, return it */
GEN
obj_insert(GEN S, int64_t K, GEN O)
{ return obj_insert_shallow(S, K, gclone(O)); }
/* as obj_insert. WITHOUT cloning (for libpari, when creating a *new* obj
 * from an existing one) */
GEN
obj_insert_shallow(GEN S, int64_t K, GEN O)
{
  GEN o, v = gel(S, lg(S)-1);
  if (typ(v) != t_VEC) pari_err_TYPE("obj_insert", S);
  o = gel(v,K);
  gel(v,K) = O; /*SIGINT: before unclone(o)*/
  if (isclone(o)) gunclone(o); return gel(v,K);
}

/* Does S [last position] contain data at position K ? Return it, or NULL */
GEN
obj_check(GEN S, int64_t K)
{
  GEN O, v = gel(S, lg(S)-1);
  if (typ(v) != t_VEC || K >= lg(v)) pari_err_TYPE("obj_check", S);
  O = gel(v,K); return isintzero(O)? NULL: O;
}

GEN
obj_checkbuild(GEN S, int64_t tag, GEN (*build)(GEN))
{
  GEN O = obj_check(S, tag);
  if (!O)
  { pari_sp av = avma; O = obj_insert(S, tag, build(S)); set_avma(av); }
  return O;
}

GEN
obj_checkbuild_prec(GEN S, int64_t tag, GEN (*build)(GEN,int64_t),
  int64_t (*pr)(GEN), int64_t prec)
{
  pari_sp av = avma;
  GEN w = obj_check(S, tag);
  if (!w || pr(w) < prec) w = obj_insert(S, tag, build(S, prec));
  set_avma(av); return gcopy(w);
}
GEN
obj_checkbuild_realprec(GEN S, int64_t tag, GEN (*build)(GEN,int64_t), int64_t prec)
{ return obj_checkbuild_prec(S,tag,build,gprecision,prec); }
GEN
obj_checkbuild_padicprec(GEN S, int64_t tag, GEN (*build)(GEN,int64_t), int64_t prec)
{ return obj_checkbuild_prec(S,tag,build,padicprec_relative,prec); }

/* Reset S [last position], freeing all clones */
void
obj_free(GEN S)
{
  GEN v = gel(S, lg(S)-1);
  int64_t i;
  if (typ(v) != t_VEC) pari_err_TYPE("obj_free", S);
  for (i = 1; i < lg(v); i++)
  {
    GEN o = gel(v,i);
    gel(v,i) = gen_0;
    gunclone_deep(o);
  }
}

/*******************************************************************/
/*                                                                 */
/*                         STACK MANAGEMENT                        */
/*                                                                 */
/*******************************************************************/
INLINE void
dec_gerepile(pari_sp *x, pari_sp av0, pari_sp av, pari_sp tetpil, size_t dec)
{
  if (*x < av && *x >= av0)
  { /* update address if in stack */
    if (*x < tetpil) *x += dec;
    else pari_err_BUG("gerepile, significant pointers lost");
  }
}

void
gerepileallsp(pari_sp av, pari_sp tetpil, int n, ...)
{
  const pari_sp av0 = avma;
  const size_t dec = av-tetpil;
  int i;
  va_list a; va_start(a, n);
  (void)gerepile(av,tetpil,NULL);
  for (i=0; i<n; i++) dec_gerepile((pari_sp*)va_arg(a,GEN*), av0,av,tetpil,dec);
  va_end(a);
}

/* Takes an array of pointers to GENs, of length n.
 * Cleans up the stack between av and tetpil, updating those GENs. */
void
gerepilemanysp(pari_sp av, pari_sp tetpil, GEN* gptr[], int n)
{
  const pari_sp av0 = avma;
  const size_t dec = av-tetpil;
  int i;
  (void)gerepile(av,tetpil,NULL);
  for (i=0; i<n; i++) dec_gerepile((pari_sp*)gptr[i], av0, av, tetpil, dec);
}

/* Takes an array of GENs (cast to longs), of length n.
 * Cleans up the stack between av and tetpil, updating those GENs. */
void
gerepilecoeffssp(pari_sp av, pari_sp tetpil, int64_t *g, int n)
{
  const pari_sp av0 = avma;
  const size_t dec = av-tetpil;
  int i;
  (void)gerepile(av,tetpil,NULL);
  for (i=0; i<n; i++,g++) dec_gerepile((pari_sp*)g, av0, av, tetpil, dec);
}

static int
dochk_gerepileupto(GEN av, GEN x)
{
  int64_t i,lx,tx;
  if (!isonstack(x)) return 1;
  if (x > av)
  {
    pari_warn(warner,"bad object %Ps",x);
    return 0;
  }
  tx = typ(x);
  if (! is_recursive_t(tx)) return 1;

  lx = lg(x);
  for (i=lontyp[tx]; i<lx; i++)
    if (!dochk_gerepileupto(av, gel(x,i)))
    {
      pari_warn(warner,"bad component %ld in object %Ps",i,x);
      return 0;
    }
  return 1;
}
/* check that x and all its components are out of stack, or have been
 * created after av */
int
chk_gerepileupto(GEN x) { return dochk_gerepileupto(x, x); }

/* print stack between avma & av */
void
dbg_gerepile(pari_sp av)
{
  GEN x = (GEN)avma;
  while (x < (GEN)av)
  {
    const int64_t tx = typ(x), lx = lg(x);
    GEN *a;

    pari_printf(" [%ld] %Ps:", x - (GEN)avma, x);
    if (! is_recursive_t(tx)) { pari_putc('\n'); x += lx; continue; }
    a = (GEN*)x + lontyp[tx]; x += lx;
    for (  ; a < (GEN*)x; a++)
    {
      if (*a == gen_0)
        pari_puts("  gen_0");
      else if (*a == gen_1)
        pari_puts("  gen_1");
      else if (*a == gen_m1)
        pari_puts("  gen_m1");
      else if (*a == gen_2)
        pari_puts("  gen_2");
      else if (*a == gen_m2)
        pari_puts("  gen_m2");
      else if (*a == ghalf)
        pari_puts("  ghalf");
      else if (isclone(*a))
        pari_printf("  %Ps (clone)", *a);
      else
        pari_printf("  %Ps [%ld]", *a, *a - (GEN)avma);
      if (a+1 < (GEN*)x) pari_putc(',');
    }
    pari_printf("\n");
  }
}
void
dbg_gerepileupto(GEN q)
{
  err_printf("%Ps:\n", q);
  dbg_gerepile((pari_sp) (q+lg(q)));
}

GEN
gerepile(pari_sp av, pari_sp tetpil, GEN q)
{
  const size_t dec = av - tetpil;
  const pari_sp av0 = avma;
  GEN x, a;

  if (dec == 0) return q;
  if ((int64_t)dec < 0) pari_err(e_MISC,"lbot>ltop in gerepile");

  /* dec_gerepile(&q, av0, av, tetpil, dec), saving 1 comparison */
  if (q >= (GEN)av0 && q < (GEN)tetpil)
    q = (GEN) (((pari_sp)q) + dec);

  for (x = (GEN)av, a = (GEN)tetpil; a > (GEN)av0; ) *--x = *--a;
  set_avma((pari_sp)x);
  while (x < (GEN)av)
  {
    const int64_t tx = typ(x), lx = lg(x);

    if (! is_recursive_t(tx)) { x += lx; continue; }
    a = x + lontyp[tx]; x += lx;
    for (  ; a < x; a++) dec_gerepile((pari_sp*)a, av0, av, tetpil, dec);
  }
  return q;
}

void
fill_stack(void)
{
  GEN x = ((GEN)pari_mainstack->bot);
  while (x < (GEN)avma) *x++ = 0xfefefefeUL;
}

void
debug_stack(void)
{
  pari_sp top = pari_mainstack->top, bot = pari_mainstack->bot;
  GEN z;
  err_printf("bot=0x%lx\ttop=0x%lx\tavma=0x%lx\n", bot, top, avma);
  for (z = ((GEN)top)-1; z >= (GEN)avma; z--)
    err_printf("%p:\t0x%lx\t%lu\n",z,*z,*z);
}

void
setdebugvar(int64_t n) { DEBUGVAR=n; }

int64_t
getdebugvar(void) { return DEBUGVAR; }

int64_t
getstack(void) { return pari_mainstack->top-avma; }

/*******************************************************************/
/*                                                                 */
/*                               timer_delay                             */
/*                                                                 */
/*******************************************************************/

#if defined(USE_CLOCK_GETTIME)
#if defined(_POSIX_THREAD_CPUTIME)
static THREAD clockid_t time_type = CLOCK_THREAD_CPUTIME_ID;
#else
static const THREAD clockid_t time_type = CLOCK_PROCESS_CPUTIME_ID;
#endif
static void
pari_init_timer(void)
{
#if defined(_POSIX_THREAD_CPUTIME)
  time_type = CLOCK_PROCESS_CPUTIME_ID;
#endif
}

void
timer_start(pari_timer *T)
{
  struct timespec t;
  clock_gettime(time_type,&t);
  T->us = t.tv_nsec / 1000;
  T->s  = t.tv_sec;
}
#elif defined(USE_GETRUSAGE)
#ifdef RUSAGE_THREAD
static THREAD int rusage_type = RUSAGE_THREAD;
#else
static const THREAD int rusage_type = RUSAGE_SELF;
#endif /*RUSAGE_THREAD*/
static void
pari_init_timer(void)
{
#ifdef RUSAGE_THREAD
  rusage_type = RUSAGE_SELF;
#endif
}

void
timer_start(pari_timer *T)
{
  struct rusage r;
  getrusage(rusage_type,&r);
  T->us = r.ru_utime.tv_usec;
  T->s  = r.ru_utime.tv_sec;
}
#elif defined(USE_FTIME)

static void
pari_init_timer(void) { }

void
timer_start(pari_timer *T)
{
  struct timeb t;
  ftime(&t);
  T->us = ((int64_t)t.millitm) * 1000;
  T->s  = t.time;
}

#else

static void
_get_time(pari_timer *T, int64_t Ticks, int64_t TickPerSecond)
{
  T->us = (int64_t) ((Ticks % TickPerSecond) * (1000000. / TickPerSecond));
  T->s  = Ticks / TickPerSecond;
}

# ifdef USE_TIMES
static void
pari_init_timer(void) { }

void
timer_start(pari_timer *T)
{
# ifdef _SC_CLK_TCK
  int64_t tck = sysconf(_SC_CLK_TCK);
# else
  int64_t tck = CLK_TCK;
# endif
  struct tms t; times(&t);
  _get_time(T, t.tms_utime, tck);
}
# elif defined(_WIN32)
static void
pari_init_timer(void) { }

void
timer_start(pari_timer *T)
{ _get_time(T, win32_timer(), 1000); }
# else
#  include <time.h>
#  ifndef CLOCKS_PER_SEC
#   define CLOCKS_PER_SEC 1000000 /* may be false on YOUR system */
#  endif
static void
pari_init_timer(void) { }

void
timer_start(pari_timer *T)
{ _get_time(T, clock(), CLOCKS_PER_SEC); }
# endif
#endif

/* round microseconds to milliseconds */
static int64_t
rndus(int64_t x) { return (x + 500) / 1000; }
static int64_t
timer_aux(pari_timer *T, pari_timer *U, void (*settime)(pari_timer *))
{
  int64_t s = T->s, us = T->us;
  settime(U); return 1000 * (U->s - s) + rndus(U->us - us);
}

/* return delay, set timer checkpoint */
int64_t
timer_delay(pari_timer *T) { return timer_aux(T, T, &timer_start); }
/* return delay, don't set checkpoint */
int64_t
timer_get(pari_timer *T) {pari_timer t; return timer_aux(T, &t, &timer_start);}

static void
timer_vprintf(pari_timer *T, const char *format, va_list args)
{
  out_puts(pariErr, "Time ");
  out_vprintf(pariErr, format,args);
  out_printf(pariErr, ": %ld\n", timer_delay(T));
  pariErr->flush();
}
void
timer_printf(pari_timer *T, const char *format, ...)
{
  va_list args; va_start(args, format);
  timer_vprintf(T, format, args);
  va_end(args);
}

int64_t
timer(void)  { static THREAD pari_timer T; return timer_delay(&T);}
int64_t
gettime(void)  { static THREAD pari_timer T; return timer_delay(&T);}

static THREAD pari_timer timer2_T, abstimer_T;
int64_t
timer2(void) {  return timer_delay(&timer2_T);}
void
msgtimer(const char *format, ...)
{
  va_list args; va_start(args, format);
  timer_vprintf(&timer2_T, format, args);
  va_end(args);
}
int64_t
getabstime(void)  { return timer_get(&abstimer_T);}

void
walltimer_start(pari_timer *ti)
{
#if defined(USE_CLOCK_GETTIME)
  struct timespec t;
  if (!clock_gettime(CLOCK_REALTIME,&t))
  { ti->s = t.tv_sec; ti->us = rndus(t.tv_nsec); return; }
#elif defined(USE_GETTIMEOFDAY)
  struct timeval tv;
  if (!gettimeofday(&tv, NULL))
  {  ti->s = tv.tv_sec; ti->us = tv.tv_usec; return; }
#elif defined(USE_FTIMEFORWALLTIME)
  struct timeb tp;
  if (!ftime(&tp))
  { ti->s = tp.time; ti->us = tp.millitm*1000; return; }
#endif
  timer_start(ti);
}
/* return delay, set timer checkpoint */
int64_t
walltimer_delay(pari_timer *T) { return timer_aux(T, T, &walltimer_start); }
/* return delay, don't set checkpoint */
int64_t
walltimer_get(pari_timer *T)
{
  pari_timer t;
  return timer_aux(T, &t, &walltimer_start);
}

static GEN
timetoi(ulong s, ulong m)
{
  pari_sp av = avma;
  return gerepileuptoint(av, addiu(muluu(s, 1000), m));
}
GEN
getwalltime(void)
{
  pari_timer ti;
  walltimer_start(&ti);
  return timetoi(ti.s, rndus(ti.us));
}

/*******************************************************************/
/*                                                                 */
/*                   FUNCTIONS KNOWN TO THE ANALYZER               */
/*                                                                 */
/*******************************************************************/
GEN
pari_version(void)
{
  const ulong mask = (1ULL<<PARI_VERSION_SHIFT) - 1;
  ulong major, minor, patch, n = paricfg_version_code;
  patch = n & mask; n >>= PARI_VERSION_SHIFT;
  minor = n & mask; n >>= PARI_VERSION_SHIFT;
  major = n;
  if (*paricfg_vcsversion) {
    const char *ver = paricfg_vcsversion;
    const char *s = strchr(ver, '-');
    char t[8];
    const int64_t len = s-ver;
    GEN v;
    if (!s || len > 6) pari_err_BUG("pari_version()"); /* paranoia */
    memcpy(t, ver, len); t[len] = 0;
    v = cgetg(6, t_VEC);
    gel(v,1) = utoi(major);
    gel(v,2) = utoi(minor);
    gel(v,3) = utoi(patch);
    gel(v,4) = stoi( atoi(t) );
    gel(v,5) = strtoGENstr(s+1);
    return v;
  } else {
    GEN v = cgetg(4, t_VEC);
    gel(v,1) = utoi(major);
    gel(v,2) = utoi(minor);
    gel(v,3) = utoi(patch);
    return v;
  }
}

/* List of GP functions: generated from the description system.
 * Format (struct entree) :
 *   char *name   : name (under GP).
 *   ulong valence: (EpNEW, EpALIAS,EpVAR, EpINSTALL)|EpSTATIC
 *   void *value  : For PREDEFINED FUNCTIONS: C function to call.
 *                  For USER FUNCTIONS: pointer to defining data (block) =
 *                   entree*: NULL, list of entree (arguments), NULL
 *                   char*  : function text
 *   int64_t menu    : which help section do we belong to
 *                   1: Standard monadic or dyadic OPERATORS
 *                   2: CONVERSIONS and similar elementary functions
 *                   3: functions related to COMBINATORICS
 *                   4: TRANSCENDENTAL functions, etc.
 *   char *code   : GP prototype, aka Parser Code (see libpari's manual)
 *                  if NULL, use valence instead.
 *   char *help   : short help text (init to NULL).
 *   void *pvalue : push_val history.
 *   int64_t arity   : maximum number of arguments.
 *   entree *next : next entree (init to NULL, used in hashing code). */
#include "init.h"
#include "default.h"
