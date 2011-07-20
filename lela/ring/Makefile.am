# Copyright (c) 2010 the LELA group
# This file is part of LELA
# see COPYING for licence
INCLUDES=-I$(top_srcdir) -I$(top_builddir) $(GMP_CFLAGS)
AM_CXXFLAGS = -O2 -Wall $(GMP_CFLAGS) $(ATLAS_CFLAGS)

pkgincludesubdir=$(pkgincludedir)/ring

BASIC_HDRS =			\
	interface.h		\
	type-wrapper.h		\
	rationals.h		\
	modular.h		\
	gf2.h

pkgincludesub_HEADERS =			\
	$(BASIC_HDRS)

noinst_LTLIBRARIES = libring.la

libring_la_SOURCES =	\
	rationals.C
