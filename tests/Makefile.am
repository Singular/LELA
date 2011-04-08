# Copyright 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

INCLUDES=-I$(top_srcdir) -I$(top_builddir)
#AM_LDFLAGS = -llapack -lcblas -latlas 
#AM_CXXFLAGS = -g -Wall -O2 
AM_CXXFLAGS = -g -Wall -DDEBUG -O0
#AM_CXXFLAGS = -g -Wall -D__LINBOX_HAVE_DGETRF -D__LINBOX_HAVE_DGETRI 
#AM_CFLAGS = -g -Wall

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = data

AM_CPPFLAGS= $(CBLAS_FLAG) $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
LDADD = $(BLAS_LIBS) $(NTL_LIBS) $(GIVARO_LIBS) $(GMP_LIBS)  $(LIDIA_LIBS) $(M4RI_LIBS) $(top_builddir)/linbox/liblinbox.la

EXTRA_DIST=\
test-common.C

# Put an entry for your test-mycomponent.C here (or in one of the other test groups).  
# Don't forget to add a test_mycomponent_SOURCES entry below.
#
# Once present, tests should remain present, being moved from one group to another.
# Note that among the available groups are NONCOMPILING_TESTS and FAILING_TESTS, which
# will be ignored by "make check" and can be temporary storage for problematic tests.
BASIC_TESTS =  \
        test-commentator        \
        test-ffpack             \
        test-gf2                \
        test-gmp-rational       \
	test-matrix		\
        test-modular            \
        test-modular-byte       \
        test-modular-double     \
        test-modular-short      \
        test-modular-int32      \
        test-modular-balanced-int32 \
        test-matrix-domain      \
	test-matrix-domain-gf2	\
        test-randiter-nonzero   \
        test-subiterator        \
        test-subvector          \
	test-bit-subvector	\
	test-hybrid-vector	\
        test-vector-domain      \
        test-unparametric-field \
	test-gauss-jordan

# a benchmarker, not to be included in check.
BENCHMARKS =                    \
        benchmark-fields        \
	benchmark-matrix-domain	\
        benchmark-blas-domain

NON_COMPILING_TESTS = \
        test-blas-domain        

EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =                               \
        $(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_commentator_SOURCES =                \
        test-commentator.C                \
        test-common.C

test_modular_SOURCES =                        \
        test-modular.C                        \
        test-generic.h                        \
        test-common.C

test_gf2_SOURCES =                        \
        test-gf2.C                        \
        test-generic.h                    \
        test-common.C

test_gmp_rational_SOURCES =               \
        test-gmp-rational.C               \
        test-generic.h                    \
        test-common.C

test_unparametric_field_SOURCES =        \
        test-unparametric-field.C        \
        test-generic.h                   \
        test-common.C

test_randiter_nonzero_SOURCES =          \
        test-randiter-nonzero.C          \
        test-common.C

test_subiterator_SOURCES =                 \
        test-subiterator.C                \
        test-common.C

test_subvector_SOURCES =                 \
        test-subvector.C                \
        test-common.C

test_vector_domain_SOURCES =                 \
        test-vector-domain.C                \
        test-vector-domain.h                \
	test-common.C

test_matrix_domain_SOURCES =                 \
        test-matrix-domain.C                \
        test-matrix-domain.h                \
        test-common.C

test_matrix_domain_gf2_SOURCES =                 \
        test-matrix-domain-gf2.C                \
        test-matrix-domain.h                \
        test-common.C

test_modular_int32_SOURCES =                \
        test-modular-int32.C                 \
        test-common.C        

test_modular_short_SOURCES =                \
        test-modular-short.C                \
        test-generic.h                        \
        test-common.C

test_modular_double_SOURCES =                \
        test-modular-double.C                \
        test-generic.h                        \
        test-common.C

test_modular_byte_SOURCES =                \
        test-modular-byte.C                \
        test-generic.h                        \
        test-common.C

test_modular_balanced_int32_SOURCES =             \
        test-modular-balanced-int32.C             \
        test-generic.h                  \
        test-common.C

test_matrix_SOURCES = \
	test-matrix.C \
	test-common.C \
	test-matrix.h

test_bit_subvector_SOURCES = \
	test-bit-subvector.C \
	test-common.C

test_hybrid_vector_SOURCES = \
	test-hybrid-vector.C \
	test-common.C

test_blas_domain_SOURCES=                \
        test-common.C   \
        test-blas-domain.C

test_ffpack_SOURCES=                \
        test-common.C   \
        test-ffpack.C

test_gauss_jordan_SOURCES = \
	test-common.C		\
	test-gauss-jordan.C

benchmark_matrix_domain_SOURCES = \
	benchmark-matrix-domain.C \
	test-common.C		\
	test-matrix-domain.h>

benchmark_matrix_domain_CXXFLAGS = ${BENCHMARK_CXXFLAGS}