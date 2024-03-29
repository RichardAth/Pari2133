/* Copyright (C) 2009  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

PARILIB_API const char* win32_basedir(void);
PARILIB_API char* win32_datadir(void);
PARILIB_API void win32_ansi_fputs(const char* s, void* f);
PARILIB_API int win32_terminal_width(void);
PARILIB_API int win32_terminal_height(void);
PARILIB_API void win32_set_codepage(void);
PARILIB_API char *win32_set_pdf_viewer(void);
PARILIB_API void win32_alarm(unsigned long long s);
PARILIB_API long win32_timer(void);
PARILIB_API long win32_nbthreads(void);
