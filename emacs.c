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
/*                           EMACS FRONTEND                        */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "int.h"

#include "gp.h"

void
init_emacs(void)
{
  GP_DATA->breakloop = 0;
  disable_color = 1;
}
