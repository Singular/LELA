# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/matrix

BASIC_HDRS =			\
	archetype.h		\
	matrix-category.h	\
	matrix-traits.h		\
	io.h			\
	io.tcc			\
	io-png.tcc		\
	raw-iterator.h		\
	sparse.h		\
	sparse.tcc		\
	transpose-matrix.h	\
	dense.h  		\
	dense-rows-matrix.h	\
	dense-zero-one.h	\
	sparse-zero-one.h	\
	sparse-zero-one.tcc	\
	m4ri-matrix.h		\
	submatrix.h		\
	matrix-domain.h		\
	matrix-domain.tcc	\
	matrix-domain-gf2.h	\
	matrix-domain-gf2.tcc	\
	matrix-domain-m4ri.h	\
	blas-matrix.h		\
	factorized-matrix.h	\
	factorized-matrix.inl	

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

