# Copyright 2011 Bradford Hovinen
# This file is part of LinBox
# See COPYING for license

INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -Wall -O2 

AM_CPPFLAGS= $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
LDADD = $(BLAS_LIBS) $(NTL_LIBS) $(GIVARO_LIBS) $(GMP_LIBS)  $(LIDIA_LIBS) $(M4RI_LIBS) $(top_builddir)/linbox/liblinbox.la

noinst_PROGRAMS = \
	convert \
	row-echelon-form \
	equal

convert_SOURCES = \
	convert.C	\
	support.C

equal_SOURCES = \
	equal.C	\
	support.C

row_echelon_form_SOURCES = \
	row-echelon-form.C	\
	support.C

noinst_HEADERS = \
	support.h
