# Copyright 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir) -I$(top_builddir)

pkgincludesubdir=$(pkgincludedir)/vector

EXTRA_DIST=			\
	vector.doxy

pkgincludesub_HEADERS=		\
	sparse.h		\
	vector-traits.h		\
	subvector.h		\
	subiterator.h		\
	reverse.h		\
	random.h                \
	stream.h		\
	bit-iterator.h		\
	bit-vector.h		\
	bit-vector.inl		\
	bit-subvector.h		\
	bit-subvector.inl	\
	bit-subvector-word-aligned.h	\
	shifted-vector.h	\
	sparse-subvector.h	\
	sparse-subvector-hybrid.h	\
	sparse-subvector-hybrid.inl	\
	vector-domain.h		\
	vector-domain.inl	\
	vector-domain-gf2.h	\
	vector-domain-gf2.inl

noinst_LTLIBRARIES = libvector.la

libvector_la_SOURCES =	\
	vector-domain-gf2.C
