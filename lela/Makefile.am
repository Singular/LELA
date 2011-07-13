# Copyright (c) 2010 the LELA group
# This file is part of LELA
# see COPYING for licence
INCLUDES=-I$(top_srcdir)

SUBDIRS=util element ring blas vector matrix randiter algorithms solutions

BASIC_HDRS =				\
	integer.h	\
	config-blas.h	lela-config.h

pkginclude_HEADERS =			\
	$(BASIC_HDRS)

lib_LTLIBRARIES=liblela.la

liblela_la_SOURCES = dummy.C

liblela_la_LIBADD = \
	util/libutil.la		\
	randiter/libranditer.la
