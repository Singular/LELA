# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -g -Wall -DDEBUG -O0

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
        test-integers		\
        test-rationals		\
        test-modular            \
        test-modular-double     \
        test-gf2                \
        test-type-wrapper-ring	\
        test-randiter-nonzero   \
	test-vector		\
        test-subiterator        \
        test-subvector          \
	test-bit-subvector	\
	test-hybrid-vector	\
	test-matrix		\
        test-blas-generic-module      \
	test-blas-generic-module-gf2	\
        test-blas-zp-module     \
	test-blas-cblas-module	\
	test-strassen-winograd	\
	test-elimination	\
	test-gauss-jordan	\
	test-splicer		\
	test-faugere-lachartre

# a benchmarker, not to be included in check.
BENCHMARKS =            \
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

test_integers_SOURCES =               \
        test-integers.C               \
        test-common.C

test_rationals_SOURCES =               \
        test-rationals.C               \
        test-common.C

test_type_wrapper_ring_SOURCES =        \
        test-type-wrapper-ring.C        \
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

test_blas_zp_module_SOURCES =   \
        test-blas-zp-module.C   \
        test-blas-level1.h           \
        test-blas-level2.h           \
        test-blas-level3.h           \
        test-common.C

test_blas_cblas_module_SOURCES =   \
        test-blas-cblas-module.C   \
        test-blas-level1.h           \
        test-blas-level2.h           \
        test-blas-level3.h           \
        test-common.C

test_modular_double_SOURCES =                \
        test-modular-double.C                \
        test-common.C

test_vector_SOURCES = \
	test-vector.C \
	test-common.C \
	test-vector.h

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

test_splicer_SOURCES = \
	test-splicer.C	\
	test-common.C

test_faugere_lachartre_SOURCES = \
	test-common.C		\
	test-faugere-lachartre.C

benchmark_blas_CXXFLAGS = -O2

benchmark_blas_SOURCES =    \
	benchmark-blas.C    \
	test-common.C	    \
	test-blas-level3.h

benchmark_matrix_domain_CXXFLAGS = ${BENCHMARK_CXXFLAGS}

noinst_HEADERS =	\
	test-common.h	\
	test-ring.h	\
	test-vector.h	\
	test-matrix.h	\
	test-blas-level1.h	\
	test-blas-level2.h	\
	test-blas-level3.h
