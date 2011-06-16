# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/matrix

BASIC_HDRS =			\
	archetype.h		\
	traits.h		\
	io.h			\
	io.tcc			\
	io-png.tcc		\
	raw-iterator.h		\
	sparse.h		\
	sparse.tcc		\
	transpose.h		\
	dense.h  		\
	dense-zero-one.h	\
	sparse-zero-one.h	\
	sparse-zero-one.tcc	\
	m4ri-matrix.h		\
	submatrix.h

NTL_HDRS =
LIDIA_HDRS = 

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

pkgincludesub_HEADERS =		\
	$(BASIC_HDRS)		\
	$(USE_NTL_HDRS)		\
	$(USE_LIDIA_HDRS)

EXTRA_DIST = 			\
	$(NTL_HDRS)		\
	$(LIDIA_HDRS)		\
	matrix.doxy

