# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/matrix

BASIC_HDRS =			\
	interface.h		\
	traits.h		\
	io.h			\
	io.tcc			\
	io-png.tcc		\
	raw-iterator.h		\
	sparse.h		\
	sparse.tcc		\
	transpose.h		\
	transpose-submatrix.h	\
	dense.h  		\
	dense-zero-one.h	\
	sparse-zero-one.h	\
	sparse-zero-one.tcc	\
	m4ri-matrix.h		\
	submatrix.h

pkgincludesub_HEADERS =		\
	$(BASIC_HDRS)		\
	$(USE_NTL_HDRS)		\
	$(USE_LIDIA_HDRS)

EXTRA_DIST = 			\
	$(NTL_HDRS)		\
	$(LIDIA_HDRS)		\
	matrix.doxy

