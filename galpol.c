/* Copyright (C) 2000-2018  The PARI group.

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

GEN
galoisnbpol(int64_t a)
{
  GEN n;
  pariFILE *F;
  char *s = stack_malloc(strlen(pari_datadir) + 11 + 20 + 1);
  sprintf_s (s, strlen(pari_datadir) + 11 + 20, 
      "%s/galpol/%lld/nb", pari_datadir, a);
  F = pari_fopengz(s);
  if (!F) pari_err_FILE("galpol file",s);
  n = gp_read_stream(F->file);
  if (!n || typ(n)!=t_INT) pari_err_FILE("galpol file [incompatible]",s);
  pari_fclose(F); return n;
}

GEN
galoisgetpol(int64_t a, int64_t b, int64_t sig)
{
  pariFILE *F;
  GEN V;
  const char *si;
  char *s;
  if (a<=0) pari_err_DOMAIN("galoisgetpol", "degree", "<=", gen_0, stoi(a));
  if (b<0) pari_err_DOMAIN("galoisgetpol", "index", "<", gen_0, stoi(b));
  if (!b) return galoisnbpol(a);
  switch(sig)
  {
    case 1: si="real"; break;
    case 2: if (a%2==0) { si="complex"; break; }
      pari_err_DOMAIN("galoisgetpol", "s", ">", gen_1, stoi(sig));
    default:
      pari_err_FLAG("galoisgetpol");
      return NULL;/*LCOV_EXCL_LINE*/
  }
  /* left on stack */
  s = stack_sprintf("%s/galpol/%lld/%lld/%s", pari_datadir, a,b,si);
  F = pari_fopengz(s);
  if (!F)
  {
    int64_t n = itos(galoisnbpol(a));
    if (b > n)
      pari_err_DOMAIN("galoisgetpol", "group index", ">", stoi(n), stoi(b));
    else pari_err_FILE("galpol file", s);
  }
  V = gp_read_stream(F->file);
  if (!V || typ(V)!=t_VEC) pari_err_FILE("galpol file", F->name);
  pari_fclose(F); return V;
}

GEN
galoisgetgroup(int64_t a, int64_t b)
{
  pariFILE *F;
  GEN V;
  char *s;
  if (a<=0) pari_err_DOMAIN("galoisgetgroup", "degree", "<=", gen_0, stoi(a));
  if (b<0) pari_err_DOMAIN("galoisgetgroup", "index", "<", gen_0, stoi(b));
  if (!b) return galoisnbpol(a);
  /* left on stack */
  s = stack_sprintf("%s/galpol/%lld/%lld/group", pari_datadir, a,b);
  F = pari_fopengz(s);
  if (!F)
  {
    int64_t n = itos(galoisnbpol(a));
    if (b > n)
      pari_err_DOMAIN("galoisgetgroup", "group index", ">", stoi(n), stoi(b));
    else pari_err_FILE("galpol file", s);
  }
  V = gp_read_stream(F->file);
  if (!V || typ(V)!=t_VEC) pari_err_FILE("galpol file", F->name);
  pari_fclose(F); return V;
}

GEN
galoisgetname(int64_t a, int64_t b)
{
  pariFILE *F;
  GEN V;
  char *s;
  if (a<=0) pari_err_DOMAIN("galoisgetname", "degree", "<=", gen_0, stoi(a));
  if (b<0) pari_err_DOMAIN("galoisgetname", "index", "<", gen_0, stoi(b));
  /* left on stack */
  s = stack_sprintf("%s/galpol/%lld/%lld/name", pari_datadir, a,b);
  F = pari_fopengz(s);
  if (!F)
  {
    int64_t n = itos(galoisnbpol(a));
    if (b > n)
      pari_err_DOMAIN("galoisgetname", "group index", ">", stoi(n), stoi(b));
    else pari_err_FILE("galpol file", s);
  }
  V = gp_read_stream(F->file);
  if (!V || typ(V)!=t_STR) pari_err_FILE("galpol file", F->name);
  pari_fclose(F); return V;
}
