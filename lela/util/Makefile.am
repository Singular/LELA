# Copyright (c) 2010 the LELA group
# This file is part of LELA
# see COPYING for licence

INCLUDES=-I$(top_srcdir) -I$(top_builddir) $(GMP_CFLAGS)
AM_CXXFLAGS=-Wall -O2

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la

libutil_la_SOURCES=	\
	timer.C		\
	error.C		\
	commentator.C	\
	debug.C		\
	splicer.C

pkgincludesub_HEADERS=\
	debug.h		\
	error.h		\
	commentator.h 	\
	timer.h		\
	splicer.h	\
	splicer.tcc	\
	double-word.h	\
	property.h
