/* Copyright (C) 2006-2008  The PARI group.

This file is part of the PARI package.

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
BEGINEXTERN
#include "parse.h"
ENDEXTERN
#include "anal.h"
#include "tree.h"

static THREAD int pari_once;
static THREAD int64_t pari_discarded;
static THREAD const char *pari_lex_start;
static THREAD GEN pari_lasterror;

static void pari_error(PARI_LTYPE *yylloc, char **lex, const char *s)
{
  (void) yylloc; (void) lex;
  if (pari_lasterror) cgiv(pari_lasterror);
  pari_lasterror=strtoGENstr(s);
}

static THREAD pari_stack s_node;
THREAD node *pari_tree;

void
pari_init_parser(void)
{
  int64_t i;
  const char *opname[]={"_||_", "_&&_", "_===_", "_==_", "_!=_", "_>=_", "_>_", "_<=_", "_<_", "_-_","_+_","_<<_", "_>>_", "_%_", "_\\/_", "_\\_", "_/_", "_*_","_^_","__","_--","_++","_-=_", "_+=_", "_<<=_", "_>>=_", "_%=_", "_\\/=_", "_\\=_", "_/=_", "_*=_","+_","-_","!_","_!","_'_","_~","[_.._]","[_|_<-_,_]","[_|_<-_,_;_]","%","%#","#_",""};

  pari_stack_init(&s_node,sizeof(*pari_tree),(void **)&pari_tree);
  pari_stack_alloc(&s_node,OPnboperator);
  parsestate_reset();
  for (i=0;i<OPnboperator;i++)
  {
    pari_tree[i].f    = Fconst;
    pari_tree[i].x    = CSTentry;
    pari_tree[i].y    = -1;
    pari_tree[i].str  = opname[i];
    pari_tree[i].len  = strlen(opname[i]);
    pari_tree[i].flags= 0;
  }
}
void
pari_close_parser(void) { pari_stack_delete(&s_node); }

void
compile_err(const char *msg, const char *str)
{
  pari_err(e_SYNTAX, msg, str, pari_lex_start);
}

void
compile_varerr(const char *str)
{
  pari_err(e_SYNTAX, "variable name expected", str, pari_lex_start);
}

void
parsestate_reset(void)
{
  s_node.n = OPnboperator;
  pari_lex_start = NULL;
  pari_once=1;
  pari_discarded=0;
  pari_lasterror=NULL;
}
void
parsestate_save(struct pari_parsestate *state)
{
  state->node = s_node.n;
  state->lex_start = pari_lex_start;
  state->once = pari_once;
  state->discarded = pari_discarded;
  state->lasterror = pari_lasterror;
}
void
parsestate_restore(struct pari_parsestate *state)
{
  s_node.n = state->node;
  pari_lex_start = state->lex_start;
  pari_once = state->once;
  pari_discarded = state->discarded;
  pari_lasterror = state->lasterror;
}

GEN
pari_compile_str(const char *lex)
{
  pari_sp ltop=avma;
  GEN code;
  struct pari_parsestate state;
  parsestate_save(&state);
  pari_lex_start = lex;
  pari_once=1;
  pari_discarded=0;
  pari_lasterror=NULL;
  if (pari_parse((char**)&lex) || pari_discarded)
  {
    if (pari_lasterror)
      compile_err(GSTR(pari_lasterror),lex-1);
    else /* should not happen */
      compile_err("syntax error",lex-1);
  }
  set_avma(ltop);
  optimizenode(s_node.n-1);
  code=gp_closure(s_node.n-1);
  parsestate_restore(&state);
  return code;
}

static int64_t
newnode(Ffunc f, int64_t x, int64_t y, struct node_loc *loc)
{
  int64_t n=pari_stack_new(&s_node);
  pari_tree[n].f=f;
  pari_tree[n].x=x;
  pari_tree[n].y=y;
  pari_tree[n].str=loc->start;
  pari_tree[n].len=loc->end-loc->start;
  pari_tree[n].flags=0;
  return n;
}

static int64_t
newconst(int64_t x, struct node_loc *loc)
{
  return newnode(Fconst,x,-1,loc);
}

static int64_t
newopcall(OPerator op, int64_t x, int64_t y, struct node_loc *loc)
{
  if (y==-1)
    return newnode(Ffunction,op,x,loc);
  else
    return newnode(Ffunction,op,newnode(Flistarg,x,y,loc),loc);
}

static int64_t
newopcall3(OPerator op, int64_t x, int64_t y, int64_t z, struct node_loc *loc)
{
  return newopcall(op,newnode(Flistarg,x,y,loc),z,loc);
}

static int64_t
countarg(int64_t n)
{
  int64_t i;
  for(i=1; pari_tree[n].f==Flistarg; i++)
    n = pari_tree[n].x;
  return i;
}

static int64_t
addcurrexpr(int64_t n, int64_t currexpr, struct node_loc *loc)
{
  int64_t y, m = n;
  while (pari_tree[m].x==OPcomprc)
  {
    y = pari_tree[m].y; if (countarg(y)==4) y = pari_tree[y].x;
    m = pari_tree[y].y;
  }
  y = pari_tree[m].y; if (countarg(y)==4) y = pari_tree[y].x;
  pari_tree[y].y = currexpr;
  pari_tree[n].str=loc->start;
  pari_tree[n].len=loc->end-loc->start;
  return n;
}

static int64_t
newintnode(struct node_loc *loc)
{
  if (loc->end-loc->start<=(int64_t)(1+LOG10_2*BITS_IN_LONG))
  {
    pari_sp ltop=avma;
    GEN g=strtoi(loc->start);
    int64_t s = itos_or_0(g), sg = signe(g);
    set_avma(ltop);
    if (sg==0 || s) return newnode(Fsmall,s,-1,loc);
  }
  return newconst(CSTint,loc);
}

static int64_t
newfunc(CSTtype t, struct node_loc *func, int64_t args, int64_t code,
                   struct node_loc *loc)
{
  int64_t name=newnode(Fentry,newconst(t,func),-1,func);
  return newnode(Fassign,name,newnode(Flambda,args,code,loc),loc);
}


