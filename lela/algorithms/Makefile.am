# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_builddir)

pkgincludesubdir=$(pkgincludedir)/algorithms

pkgincludesub_HEADERS =		\
	strassen-winograd.h	\
	strassen-winograd.tcc	\
	pivot-strategy.h	\
	pivot-strategy.tcc	\
	elimination.h		\
	elimination.tcc		\
	gauss-jordan.h 		\
	gauss-jordan.tcc	\
	faugere-lachartre.h	\
	faugere-lachartre.tcc

AM_CPPFLAGS= $(CBLAS_FLAG) $(GMP_CFLAGS)

