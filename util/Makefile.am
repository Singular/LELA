# Copyright 2011 Bradford Hovinen
# This file is part of LELA
# See COPYING for license

INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -Wall -O2 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la

noinst_PROGRAMS = \
	row-echelon-form \
	equal		\
	diff
#	convert

# convert_SOURCES = \
# 	convert.C	\
# 	support.C

equal_SOURCES = \
	equal.C	\
	support.C

diff_SOURCES = \
	diff.C	\
	support.C

row_echelon_form_SOURCES = \
	row-echelon-form.C	\
	support.C

noinst_HEADERS = \
	support.h
