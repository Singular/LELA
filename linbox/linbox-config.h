#ifndef _LINBOX_LINBOX_CONFIG_H
#define _LINBOX_LINBOX_CONFIG_H 1
 
/* linbox/linbox-config.h. Generated automatically at end of configure. */
/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.in by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef __LINBOX_AC_APPLE_UNIVERSAL_BUILD */

/* Enable Autoimplementation of dgetri routine with dtrti and dtrsm */
/* #undef __LINBOX_AUTOIMPLEMENT_DGETRI */

/* Define if BLAS routines are available */
#ifndef __LINBOX_BLAS_AVAILABLE 
#define __LINBOX_BLAS_AVAILABLE  /**/ 
#endif

/* Define if GMP is version 3.xxx */
/* #undef __LINBOX_GMP_VERSION_3 */

/* Define that architecture uses big endian storage */
/* #undef __LINBOX_HAVE_BIG_ENDIAN */

/* Define if BLAS is installed */
#ifndef __LINBOX_HAVE_BLAS 
#define __LINBOX_HAVE_BLAS  1 
#endif

/* Define if C interface to BLAS is available */
/* #undef __LINBOX_HAVE_CBLAS */

/* Define if dgetrf is available */
/* #undef __LINBOX_HAVE_DGETRF */

/* Define if dgetri is available */
/* #undef __LINBOX_HAVE_DGETRI */

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef __LINBOX_HAVE_DLFCN_H 
#define __LINBOX_HAVE_DLFCN_H  1 
#endif

/* Define if dtrtri is available */
/* #undef __LINBOX_HAVE_DTRTRI */

/* Define if GIVARO is installed */
#ifndef __LINBOX_HAVE_GIVARO 
#define __LINBOX_HAVE_GIVARO  1 
#endif

/* Define if GMP is installed */
#ifndef __LINBOX_HAVE_GMP 
#define __LINBOX_HAVE_GMP  1 
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef __LINBOX_HAVE_INTTYPES_H 
#define __LINBOX_HAVE_INTTYPES_H  1 
#endif

/* Define if LIDIA is installed */
/* #undef __LINBOX_HAVE_LIDIA */

/* Define that architecture uses little endian storage */
#ifndef __LINBOX_HAVE_LITTLE_ENDIAN 
#define __LINBOX_HAVE_LITTLE_ENDIAN  1 
#endif

/* Define if MAPLE is installed */
/* #undef __LINBOX_HAVE_MAPLE */

/* Define to 1 if you have the <memory.h> header file. */
#ifndef __LINBOX_HAVE_MEMORY_H 
#define __LINBOX_HAVE_MEMORY_H  1 
#endif

/* Define if NTL is installed */
/* #undef __LINBOX_HAVE_NTL */

/* Define if SACLIB is installed */
/* #undef __LINBOX_HAVE_SACLIB */

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef __LINBOX_HAVE_STDINT_H 
#define __LINBOX_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef __LINBOX_HAVE_STDLIB_H 
#define __LINBOX_HAVE_STDLIB_H  1 
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef __LINBOX_HAVE_STRINGS_H 
#define __LINBOX_HAVE_STRINGS_H  1 
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef __LINBOX_HAVE_STRING_H 
#define __LINBOX_HAVE_STRING_H  1 
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef __LINBOX_HAVE_SYS_STAT_H 
#define __LINBOX_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef __LINBOX_HAVE_SYS_TYPES_H 
#define __LINBOX_HAVE_SYS_TYPES_H  1 
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef __LINBOX_HAVE_UNISTD_H 
#define __LINBOX_HAVE_UNISTD_H  1 
#endif

/* Canonical 16-bit data type */
#ifndef __LINBOX_INT16 
#define __LINBOX_INT16  short 
#endif

/* Canonical 32-bit data type */
#ifndef __LINBOX_INT32 
#define __LINBOX_INT32  int 
#endif

/* Canonical 64-bit data type */
#ifndef __LINBOX_INT64 
#define __LINBOX_INT64  long 
#endif

/* Canonical 8-bit data type */
#ifndef __LINBOX_INT8 
#define __LINBOX_INT8  char 
#endif

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#ifndef __LINBOX_LT_OBJDIR 
#define __LINBOX_LT_OBJDIR  ".libs/" 
#endif

/* define is the version of Maple have access function to gmp data */
/* #undef __LINBOX_MAPLE_GMP_ACCESS */

/* Name of package */
#ifndef __LINBOX_PACKAGE 
#define __LINBOX_PACKAGE  "linbox" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef __LINBOX_PACKAGE_BUGREPORT 
#define __LINBOX_PACKAGE_BUGREPORT  "linbox-use@googlegroups.com" 
#endif

/* Define to the full name of this package. */
#ifndef __LINBOX_PACKAGE_NAME 
#define __LINBOX_PACKAGE_NAME  "linbox" 
#endif

/* Define to the full name and version of this package. */
#ifndef __LINBOX_PACKAGE_STRING 
#define __LINBOX_PACKAGE_STRING  "linbox 1.1.7" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef __LINBOX_PACKAGE_TARNAME 
#define __LINBOX_PACKAGE_TARNAME  "linbox" 
#endif

/* Define to the home page for this package. */
#ifndef __LINBOX_PACKAGE_URL 
#define __LINBOX_PACKAGE_URL  "" 
#endif

/* Define to the version of this package. */
#ifndef __LINBOX_PACKAGE_VERSION 
#define __LINBOX_PACKAGE_VERSION  "1.1.7" 
#endif

/* The size of `char', as computed by sizeof. */
#ifndef __LINBOX_SIZEOF_CHAR 
#define __LINBOX_SIZEOF_CHAR  1 
#endif

/* The size of `int', as computed by sizeof. */
#ifndef __LINBOX_SIZEOF_INT 
#define __LINBOX_SIZEOF_INT  4 
#endif

/* The size of `long', as computed by sizeof. */
#ifndef __LINBOX_SIZEOF_LONG 
#define __LINBOX_SIZEOF_LONG  8 
#endif

/* The size of `long long', as computed by sizeof. */
#ifndef __LINBOX_SIZEOF_LONG_LONG 
#define __LINBOX_SIZEOF_LONG_LONG  8 
#endif

/* The size of `short', as computed by sizeof. */
#ifndef __LINBOX_SIZEOF_SHORT 
#define __LINBOX_SIZEOF_SHORT  2 
#endif

/* The size of `__int64', as computed by sizeof. */
#ifndef __LINBOX_SIZEOF___INT64 
#define __LINBOX_SIZEOF___INT64  0 
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef __LINBOX_STDC_HEADERS 
#define __LINBOX_STDC_HEADERS  1 
#endif

/* Define if optimized threshold for Strassen-Winograd matrix multiplication
   is available */
/* #undef __LINBOX_STRASSEN_OPTIMIZATION */

/* Version number of package */
#ifndef __LINBOX_VERSION 
#define __LINBOX_VERSION  "1.1.7" 
#endif

/* optimized threshold for switching to strassen matrix multiplication */
/* #undef __LINBOX_WINOTHRESHOLD */

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define if Expat is installed */
/* #undef __LINBOX_XMLENABLED */
 
/* once: _LINBOX_LINBOX_CONFIG_H */
#endif
