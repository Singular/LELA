/* linbox/blas/level3-m4ri.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_M4RI_H
#define __BLAS_LEVEL3_M4RI_H

#ifndef __LINBOX_HAVE_M4RI
#  error "This header file requires that LinBox be configured with libm4ri enabled. Please ensure that libm4ri is properly installed and re-run configure."
#endif

#include <m4ri/m4ri.h>

#include "linbox/blas/context.h"
#include "linbox/blas/level3-ll.h"
#include "linbox/matrix/traits.h"
#include "linbox/matrix/m4ri-matrix.h"
#include "linbox/matrix/submatrix.h"

#ifdef __BLAS_LEVEL3_H
#  warning "linbox/blas/level3.h has already been included by this point. The M4RI-specialisations may not work."
#endif // __BLAS_LEVEL3_H

namespace LinBox
{

namespace BLAS3
{

M4RIMatrixBase &copy_impl (const GF2 &F, M4RIModule &M, const M4RIMatrixBase &A, M4RIMatrixBase &B,
			   MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ mzd_copy (B._rep, A._rep); return B; }

DenseMatrix<bool> &copy_impl (const GF2 &F, M4RIModule &M, const DenseMatrix<bool> &A, DenseMatrix<bool> &B,
			   MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _copy (F, M, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }

DenseMatrix<bool> &copy_impl (const GF2 &F, M4RIModule &M, const M4RISubmatrix &A, DenseMatrix<bool> &B,
			   MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _copy (F, M, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }

M4RISubmatrix &copy_impl (const GF2 &F, M4RIModule &M, const M4RISubmatrix &A, M4RISubmatrix &B,
			  MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _copy (F, M, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }

M4RIMatrixBase &scal_impl (const GF2 &F, M4RIModule &M, bool a, M4RIMatrixBase &A, MatrixCategories::RowMatrixTag);

DenseMatrix<bool> &scal_impl (const GF2 &F, M4RIModule &M, bool a, DenseMatrix<bool> &A, MatrixCategories::RowMatrixTag)
	{ _scal (F, M, a, (M4RIMatrixBase &) A); return A; }

M4RIMatrixBase &axpy_impl (const GF2 &F, M4RIModule &M, bool a, const M4RIMatrixBase &A, M4RIMatrixBase &B,
			   MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

DenseMatrix<bool> &axpy_impl (const GF2 &F, M4RIModule &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _axpy (F, M, a, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B); return B; }

M4RIMatrixBase &gemm_impl (const GF2 &F, M4RIModule &M,
			   bool a, const M4RIMatrixBase &A, const M4RIMatrix &B, bool b, M4RIMatrixBase &C,
			   MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

DenseMatrix<bool> &gemm_impl (const GF2 &F, M4RIModule &M,
			      bool a, const DenseMatrix<bool> &A, const DenseMatrix<bool> &B, bool b, DenseMatrix<bool> &C,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _gemm (F, M, a, (const M4RIMatrixBase &) A, (const M4RIMatrix &) B, b, (M4RIMatrixBase &) C); return C; }

DenseMatrix<bool> &gemm_impl (const GF2 &F, M4RIModule &M,
			      bool a, const M4RISubmatrix &A, const DenseMatrix<bool> &B, bool b, DenseMatrix<bool> &C,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _gemm (F, M, a, (const M4RIMatrixBase &) A, (const M4RIMatrix &) B, b, (M4RIMatrixBase &) C); return C; }

M4RISubmatrix &gemm_impl (const GF2 &F, M4RIModule &M,
			  bool a, const M4RISubmatrix &A, const DenseMatrix<bool> &B, bool b, M4RISubmatrix &C,
			  MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _gemm (F, M, a, (const M4RIMatrixBase &) A, (const M4RIMatrix &) B, b, (M4RIMatrixBase &) C); return C; }

DenseMatrix<bool> &trmm_impl (const GF2 &F, M4RIModule &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B, TriangularMatrixType type, bool diagIsOne,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

M4RIMatrixBase &trsm_impl (const GF2 &F, M4RIModule &M, bool a, const M4RIMatrixBase &A, M4RIMatrixBase &B, TriangularMatrixType type, bool diagIsOne,
			   MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

DenseMatrix<bool> &trsm_impl (const GF2 &F, AllModules<GF2> &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B, TriangularMatrixType type, bool diagIsOne,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ _trsm (F, (M4RIModule &) M, a, (const M4RIMatrixBase &) A, (M4RIMatrixBase &) B, type, diagIsOne); return B; }

template <class Iterator>
M4RIMatrixBase &permute_rows_impl (const GF2 &F, M4RIModule &M, Iterator P_begin, Iterator P_end, M4RIMatrixBase &A, MatrixCategories::RowMatrixTag);

template <class Iterator>
M4RIMatrixBase &permute_cols_impl (const GF2 &F, M4RIModule &M, Iterator P_begin, Iterator P_end, M4RIMatrixBase &A, MatrixCategories::RowMatrixTag);

bool equal_impl (const GF2 &F, M4RIModule &M, const M4RIMatrix &A, const M4RIMatrix &B,
		 MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ return mzd_equal (A._rep, B._rep); }

bool equal_impl (const GF2 &F, M4RIModule &M, const DenseMatrix<bool> &A, const DenseMatrix<bool> &B,
		 MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
	{ return _equal (F, M, (const M4RIMatrix &) A, (const M4RIMatrix &) B); }

bool is_zero_impl (const GF2 &F, M4RIModule &M, const M4RIMatrixBase &A, MatrixCategories::RowMatrixTag)
	{ return mzd_is_zero (A._rep); }

bool is_zero_impl (const GF2 &F, M4RIModule &M, const DenseMatrix<bool> &A, MatrixCategories::RowMatrixTag)
	{ return _is_zero (F, M, (const M4RIMatrixBase &) A); }

bool is_zero_impl (const GF2 &F, M4RIModule &M, const M4RISubmatrix &A, MatrixCategories::RowMatrixTag)
	{ return _is_zero (F, M, (const M4RIMatrixBase &) A); }

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_M4RI_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
