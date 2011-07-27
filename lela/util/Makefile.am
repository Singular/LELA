# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

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
