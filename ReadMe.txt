This project takes the source code of Pari/GP and converts it so that it can
be compiled and built using Microsft Visual Studio.

The main benefit of this is that it is then possible to link to ParilLib and 
call Pari functions in the normal way from programs compiled using Visual Studio.

It is also possible to link dynamically to the dll that can be downloaded from 
the Pari/GP website, and then create a function pointer for every Pari function 
you want to call. It is not possible to link statically to the downloaded dll
because it was compiled and built using gcc and MSYS2, so is not compatible
with the Microsoft linker.

To download Pari/GP go to http://pari.math.u-bordeaux.fr/download.html

There are a number of reference manuals in PDF format that can be downloaded.
The most importand (and by far the largest!) are users.pdf and libpari.pdf

users.pdf is mainly targeted at users of GP, whereas parilib.pdf explains
how to use the pari library from a C program.

Many changes to the source code were needed to compile with Visual Studio:

In Visual Studio C a 'long' is 32 bits and in GCC a long is 64 bits. The trick
in the original code of using a macro to define long as long long is both
confusing and caused errors when using Microsoft Windows headers.

Several things were changed:

'long' was replaced everywhere by 'int64_t'
Numeric constants ending in L were changed to end with LL e.g.
~0UL becomes ~0ULL, etc

In format strings for printf and related functions %d or %ld generally
becomes %lld, %lu becone %llu and so on.

The level1.h header was replaced by a level1.c module. This also gets round
the problem that '__extension__' does not exist in visual studio C.

Variable length arrays in a couple of places had to be replaced with 
an appropriate pointer and storage was allocated using malloc.
elements affected: ratpoints.c

In many places strcat was replaced with strcat_s and strcpy with strcpy_s.

In other places "#pragma warning (disable: 4996)" was used to allow strcat,
strcpy etc to be used.

There is no 'configure' process here. 
MPIR is required.

In a few places intrinsics were used to replace C code to improve performance.
the intrinsics used are:
_BitScanReverse64       (function bfffo)
_udiv128                (function divll)
_umul128                (function addmul & mulll)
_addcarry_u64           (function addmul)

Notes on building:

To build GP separately from PariLib compile and linkonly the following elements 
in GP:

emacs.c
gp.c
texmacs.c
whatnow.c

Compile and link all the other source elements separately into a .obj or .dll 
file, then include this .obj or .dll file when linking GP.

This .obj or .dll can also be linked into any program that uses parilib functions.

When using GP it is important that the gprc.txt file is set up. By default GP 
expects this file to be in the same folder as the gp.exe file.

currently some of the extended help functions in GP do not work. 
I cannot get the perl script to work. As a workaround the ?? command to open the 
users.pdf file, and the ??  tutorial / refcard / libpari 
(tutorial/reference card/libpari manual) should work provided that Acrobat Reader 
and the pdf files are installed in the expected places. The paths can be 
specified in the gprc.txt file.
e.g 
docpath = "C:\Program Files(x86)\Pari64-2-13-2\doc\"
acrobat = "C:\Program Files\Adobe\Acrobat DC\Acrobat\Acrobat.exe"

Note that in the gprc.txt file spaces in a path name require that the path name be
enclosed with quotes " "


I am not aware of any other problems.

Pari GP comand line:
Usage: <path>\PariGp.exe [options] [GP files]
Available Command line options for GP:
  [-f,--fast]           Fast start: do not read .gprc
  [-q,--quiet]          Quiet mode: do not print banner and history numbers
  [-s stacksize]        Start with the PARI stack of given size (in bytes)
  [--default key=val]   Execute default(key,val) on startup
  [--emacs]             Run as if in Emacs shell
  [--help]              Print this message
  [--test]              Test mode. No history, wrap long lines (bench only)
  [--texmacs]           Run as if using TeXmacs frontend
  [--version]           Output version info and exit
  [--version-short]     Output version number and exit
  
  Defaults set in the command line by --default over-ride settings in the gprc file.

  Quiet mode can only be set from the command line, not by default() or \d.
