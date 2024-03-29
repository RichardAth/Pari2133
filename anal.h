#pragma once
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

/*************************************************************************/
/*                                                                       */
/*                 Declarations specific to the analyzer                 */
/*                                                                       */
/*************************************************************************/
BEGINEXTERN
/* functions */
void   changevalue(entree *ep, GEN val);
void    freeep(entree *ep);
void   pari_fill_hashtable(entree **table, entree *ep);

void compile_err(const char *msg, const char *str);
void compile_varerr(const char *str);

#ifdef STACK_CHECK
extern THREAD void *PARI_stack_limit;
#endif

extern entree  **varentries;

struct node_loc
{
  const char *start,*end;
};

union token_value { int64_t val; };

int pari_lex(union token_value *yylval, struct node_loc *yylloc, char **lex);
int pari_parse(char **lex);
entree* fetch_entry_raw(const char *s, int64_t len);
entree* fetch_entry(const char *s);
void optimizenode(int64_t n);
void push_frame(GEN C, int64_t lpc, int64_t flag);
GEN  gp_closure(int64_t n);
int64_t eval_mnemonic(GEN str, const char *tmplate);

ENDEXTERN
