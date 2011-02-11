# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir)

SUBDIRS=util element field randiter vector matrix algorithms fflas ffpack

BASIC_HDRS =				\
	integer.h	\
	config-blas.h	linbox-config.h

pkginclude_HEADERS =			\
	$(BASIC_HDRS)

lib_LTLIBRARIES=liblinbox.la

liblinbox_la_SOURCES = dummy.C

liblinbox_la_LIBADD = \
	util/libutil.la		\
	randiter/libranditer.la \
	algorithms/libalgorithms.la
