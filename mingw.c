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

/* Originally written by Vasili Burdo */

#include "pwinver.h"
#include <windows.h>
#include <stdio.h>

/* definition below is required when parilib is built as a .dll file
otherwise comment it out and replace with"#DEFINE PARILIB_API " */

#ifdef _WIN32
#ifdef PARILIBDLL_EXPORTS
#define PARILIB_API __declspec(dllexport)
#else
#ifdef _LIB
#define PARILIB_API
#else
#define PARILIB_API __declspec(dllimport)
#endif
#endif
#else
#define PARILIB_API
#endif

#include "mingw.h"

/* for Microsoft Windows Visual Studio, we need to supress error 4996
to allow strcpy, sprintf, strcpy, getenv etc  which Microsoft consider to be unsafe
(no checks for buffer overflow) */
#ifdef _WIN32
#pragma warning (disable: 4996)
#endif

static const char * pariwin32_basedir = NULL;

const char*
win32_basedir(void)
{
  if (!pariwin32_basedir)
  {
    char basedir[1024];
    char* slash;
    GetModuleFileNameA(0, basedir, sizeof(basedir) );
    slash = strrchr(basedir, '\\');
    if (slash) slash[1] = 0;
    pariwin32_basedir = _strdup(basedir);
  }
  return pariwin32_basedir;
}

char*
win32_datadir(void)
{
  char datadir[1029];
  const char * basedir = win32_basedir();
  sprintf_s(datadir, sizeof(datadir), "%sdata", basedir);
  return _strdup(datadir);
}

static WORD
win32_console_color(unsigned long c)
{
  int shift, intense = 0;
  if( c >= 30 && c <= 37 ) { 
      shift = 0; c -= 30; } 
  else if( c >= 40 && c <= 47 ) { 
      shift = 4; c -= 40; } 
  else if( c >= 90 && c <= 97 ) { 
      shift = 0; intense = 8; c -= 90; } 
  else if(c >= 100 && c <= 107) { 
      shift = 4; intense = 8; c -= 100; } 
  else
    return 0;

  WORD w = 0;
  switch(c) {
  case 0: w = 0; break; /* black      */
  case 1: w = 4; break; /* red        */
  case 2: w = 2; break; /* green      */
  case 3: w = 6; break; /* yellow RG  */
  case 4: w = 1; break; /* blue       */
  case 5: w = 5; break; /* magenta RB */
  case 6: w = 3; break; /* cyan GB    */
  case 7: w = 7; break; /* white RGB  */
  }
  return (w|intense) << shift;
}

void
win32_ansi_fputs(const char* s, void* f)
{
  WORD color;
  unsigned long c[3];
  long nbarg;
  if( !(f == stdout || f == stderr) ) {
    fputs(s,f);
    return;
  }

  while(1) {
    char *p;
    p = strstr(s, "\x1b[");
    if( p > s )
      fwrite(s,p-s,1,f);

    if( p )
      p += 2;
    else {
      fputs(s,f);
      return;
    }
    nbarg = 0;
    c[nbarg++] = strtoul(p,&p,10);
    if( *p == ';' ) c[nbarg++] = strtoul(p+1,&p,10);
    if( *p == ';' ) c[nbarg++] = strtoul(p+1,&p,10);
    if( *p++ == 'm' ) {
      switch(nbarg)
      {
      case 1:
        color = 7;
        break;
      case 2:
        color = win32_console_color(c[1]);
        if (c[0]&4) color |= 0x8000;
        break;
      case 3:
        color = win32_console_color(c[1]) | win32_console_color(c[2]);
        if (c[0]&4) color |= 0x8000;
      }
      fflush(f);
      SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),color);
    }
    s = p;
  }
}

int
win32_terminal_width(void)
{
  CONSOLE_SCREEN_BUFFER_INFO sbi;
  if (!GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &sbi))
    return 0;
  return sbi.srWindow.Right - sbi.srWindow.Left + 1;
}

int
win32_terminal_height(void)
{
  CONSOLE_SCREEN_BUFFER_INFO sbi;
  if (!GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &sbi))
    return 0;
  return sbi.srWindow.Bottom - sbi.srWindow.Top + 1;
}

void
win32_set_codepage(void)
{
  SetConsoleCP( GetACP() );
  SetConsoleOutputCP( GetACP() );
}

char *
win32_set_pdf_viewer(void)
{
  char *s = getenv("GP_PDF_VIEWER");
  char* buf = s;
  if (!s)
  {
    HKEY handle;
    const char *key = "AcroExch.Document\\shell\\open\\command";
    const long SZ = 512;
    char str[512];
    int status;
    DWORD L = SZ;
    DWORD type;

    LSTATUS rv = RegOpenKeyExA(HKEY_CLASSES_ROOT, key, 0, KEY_READ, &handle);
    if (rv == ERROR_SUCCESS) {
        status = RegQueryValueExA(handle, NULL, 0, &type, (LPBYTE)str, &L);
        RegCloseKey(handle);
        if (status != ERROR_SUCCESS) 
            return NULL;
        /* str now contains the text from the registry key */
        buf = malloc(strlen(str) + 16); /*must not be freed*/
        if (str[0] = '"') {
             char *cp = strchr(str + 1, '"');
            if (cp != NULL)
                cp[1] = '\0';   /* replace char after 2nd " with null terminator */
        }
        sprintf_s(buf, strlen(str) + 16, "GP_PDF_VIEWER=%s", str);
        _putenv(buf);
    }
  }
  return buf;
}

PARILIB_API extern int win32ctrlc, win32alrm;
static HANDLE hTimerQueue = NULL;

static void CALLBACK
win32_cb_alarm(void *lpParam, BOOLEAN TimerOrWaitFired)
{
  (void) lpParam; (void) TimerOrWaitFired;
  win32ctrlc++;
  win32alrm = 1;
}

void
win32_alarm(unsigned long long s)
{
  if (hTimerQueue)
  {
    HANDLE oldhTimerQueue = hTimerQueue;
    hTimerQueue = NULL;
    DeleteTimerQueue(oldhTimerQueue);
  }
  if (s)
  {
    void *arg = NULL;
    HANDLE hTimer = NULL;
    hTimerQueue = CreateTimerQueue();
    CreateTimerQueueTimer( &hTimer, hTimerQueue,
        (WAITORTIMERCALLBACK)win32_cb_alarm, &arg , s*1000, 0, 0);
  }
}

#define WIN32_FILETIME_PER_MILLISECOND  10000

long
win32_timer(void)
{
  FILETIME lpCreation, lpExit, lpKernel, lpUser;
  LARGE_INTEGER time;
  GetProcessTimes(
    GetCurrentProcess(),
    &lpCreation, &lpExit, &lpKernel, &lpUser
    );
  time.HighPart = lpUser.dwHighDateTime;
  time.LowPart = lpUser.dwLowDateTime;
  time.QuadPart /= WIN32_FILETIME_PER_MILLISECOND;
  return time.LowPart;
}

long
win32_nbthreads(void)
{
  SYSTEM_INFO system_info;
  GetSystemInfo(&system_info);
  return system_info.dwNumberOfProcessors;
}
