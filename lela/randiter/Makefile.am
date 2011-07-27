# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

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

