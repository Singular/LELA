# Copyright 2010 the LELA group
# This file is part of LELA
# see COPYING for licence

INCLUDES=-I$(top_srcdir) -I$(top_builddir)
#AM_LDFLAGS = -llapack -lcblas -latlas 
#AM_CXXFLAGS = -g -Wall -O2 
AM_CXXFLAGS = -g -Wall -DDEBUG -O0
#AM_CXXFLAGS = -g -Wall -D__LELA_HAVE_DGETRF -D__LELA_HAVE_DGETRI 
#AM_CFLAGS = -g -Wall

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = data

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la

# Put an entry for your test-mycomponent.C here (or in one of the other test groups).  
# Don't forget to add a test_mycomponent_SOURCES entry below.
#
# Once present, tests should remain present, being moved from one group to another.
# Note that among the available groups are NONCOMPILING_TESTS and FAILING_TESTS, which
# will be ignored by "make check" and can be temporary storage for problematic tests.
BASIC_TESTS =  \
        test-commentator        \
        test-gf2                \
        test-gmp-rational       \
	test-matrix		\
        test-modular            \
        test-modular-double     \
        test-blas-generic-module      \
	test-blas-generic-module-gf2	\
        test-randiter-nonzero   \
        test-subiterator        \
        test-subvector          \
	test-bit-subvector	\
	test-hybrid-vector	\
        test-unparametric-ring	\
	test-strassen-winograd	\
	test-elimination	\
	test-gauss-jordan	\
	test-faugere-lachartre

# a benchmarker, not to be included in check.
BENCHMARKS =            \
        benchmark-rings \
	benchmark-blas

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
        test-common.C

test_gf2_SOURCES =                        \
        test-gf2.C                        \
        test-common.C

test_gmp_rational_SOURCES =               \
        test-gmp-rational.C               \
        test-common.C

test_unparametric_ring_SOURCES =        \
        test-unparametric-ring.C        \
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

test_blas_generic_module_SOURCES =   \
        test-blas-generic-module.C   \
        test-blas-level1.h           \
        test-blas-level2.h           \
        test-blas-level3.h           \
        test-common.C

test_blas_generic_module_gf2_SOURCES = \
        test-blas-generic-module-gf2.C \
        test-blas-level1.h             \
        test-blas-level2.h             \
        test-blas-level3.h             \
        test-common.C

test_modular_double_SOURCES =                \
        test-modular-double.C                \
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

test_strassen_winograd_SOURCES = \
	test-common.C		\
	test-strassen-winograd.C

test_elimination_SOURCES = \
	test-common.C		\
	test-elimination.C

test_gauss_jordan_SOURCES = \
	test-common.C		\
	test-gauss-jordan.C

test_faugere_lachartre_SOURCES = \
	test-common.C		\
	test-faugere-lachartre.C

benchmark_blas_SOURCES =    \
	benchmark-blas.C    \
	test-common.C	    \
	test-blas-level3.h

benchmark_matrix_domain_CXXFLAGS = ${BENCHMARK_CXXFLAGS}