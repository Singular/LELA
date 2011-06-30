# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir) -I$(top_builddir) $(GMP_CFLAGS)
pkgincludesubdir=$(pkgincludedir)/randiter

AM_CXXFLAGS = -O2 -Wall -DDISABLE_COMMENTATOR $(GMP_CFLAGS) $(ATLAS_CFLAGS)

BASIC_HDRS =				\
	interface.h			\
	gmp-rational.h			\
	modular.h			\
	unparametric.h			\
	gf2.h				\
	nonzero.h			\
	mersenne-twister.h

pkgincludesub_HEADERS =			\
	$(BASIC_HDRS)

noinst_LTLIBRARIES = libranditer.la

libranditer_la_SOURCES =	\
	mersenne-twister.C

