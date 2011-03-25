# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir) -I$(top_builddir)

AM_CXXFLAGS=-Wall -O2

SUBDIRS=formats

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la

libutil_la_SOURCES=	\
	timer.C		\
	error.C		\
	commentator.C	\
	debug.C

pkgincludesub_HEADERS=\
	debug.h				\
	error.h				\
	commentator.h 			\
	timer.h				\
	field-axpy.h			\
	matrix-stream.h			\
	matrix-stream.inl
