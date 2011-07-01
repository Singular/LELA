INCLUDES=-I$(top_srcdir) -I$(top_builddir)

pkgincludesubdir=$(pkgincludedir)/algorithms

pkgincludesub_HEADERS =		\
	strassen-winograd.h	\
	strassen-winograd.tcc	\
	elimination.h		\
	elimination.tcc		\
	gauss-jordan.h 		\
	gauss-jordan.tcc	\
	faugere-lachartre.h	\
	faugere-lachartre.tcc

AM_CPPFLAGS= $(CBLAS_FLAG) $(GMP_CFLAGS)

