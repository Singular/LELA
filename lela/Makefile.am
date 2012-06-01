# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_builddir)

SUBDIRS=util element ring blas vector matrix randiter algorithms solutions

BASIC_HDRS =		\
	integer.h	\
	cblas.h		\
	lela-config.h

pkginclude_HEADERS =			\
	$(BASIC_HDRS)

lib_LTLIBRARIES=liblela.la

liblela_la_SOURCES = dummy.C

liblela_la_LIBADD = \
	util/libutil.la		\
	ring/libring.la		\
	randiter/libranditer.la
