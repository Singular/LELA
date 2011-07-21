# Copyright (c) 2010 the LELA group
# This file is part of LELA
# see COPYING for licence

INCLUDES=-I$(top_srcdir) -I$(top_builddir) $(GMP_CFLAGS)
AM_CXXFLAGS = -O2 -Wall

pkgincludesubdir=$(pkgincludedir)/randiter

BASIC_HDRS =				\
	interface.h			\
	integers.h			\
	rationals.h			\
	modular.h			\
	type-wrapper.h			\
	gf2.h				\
	nonzero.h			\
	mersenne-twister.h

pkgincludesub_HEADERS =			\
	$(BASIC_HDRS)

noinst_LTLIBRARIES = libranditer.la

libranditer_la_SOURCES =	\
	mersenne-twister.C

