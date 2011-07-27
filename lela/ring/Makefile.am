# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_builddir) $(GMP_CFLAGS)
AM_CXXFLAGS = -O2 -Wall

pkgincludesubdir=$(pkgincludedir)/ring

BASIC_HDRS =			\
	interface.h		\
	type-wrapper.h		\
	integers.h		\
	rationals.h		\
	modular.h		\
	gf2.h

pkgincludesub_HEADERS =			\
	$(BASIC_HDRS)

noinst_LTLIBRARIES = libring.la

libring_la_SOURCES =	\
	rationals.C
