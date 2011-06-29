/* linbox/blas/level2-cblas.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Cblas implementations of level 2 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_CBLAS_H
#define __BLAS_LEVEL2_CBLAS_H

#include <algorithm>

#ifndef __LINBOX_BLAS_AVAILABLE
#  error "A working installation of BLAS is required to use this header-file."
#endif // __LINBOX_BLAS_AVAILABLE

#include <cblas.h>

#include "linbox/blas/context.h"
#include "linbox/ring/unparametric.h"
#include "linbox/vector/traits.h"
#include "linbox/matrix/traits.h"
#include "linbox/matrix/dense.h"

namespace LinBox
{

namespace BLAS2
{

template <class Vector1, class Vector2>
Vector2 &gemv_impl (const UnparametricRing<float> &F, BLASModule &M,
		    float a, const DenseMatrix<float> &A, const Vector1 &x, float b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowColMatrixTag,
		    VectorRepresentationTypes::Dense,
		    VectorRepresentationTypes::Dense)
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());
	cblas_sgemv (CblasRowMajor, CblasNoTrans, A.rowdim (), A.coldim (), a, &A[0][0], A.disp (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
	return y;
}

template <class Vector>
Vector &trmv_impl (const UnparametricRing<float> &F, BLASModule &M, const DenseMatrix<float> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowColMatrixTag,
		   VectorRepresentationTypes::Dense)
{
	linbox_check (A.rowdim () == x.size ());
	linbox_check (A.coldim () == x.size ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_strmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
	return x;
}

template <class Vector>
Vector &trsv_impl (const UnparametricRing<float> &F, BLASModule &M, const DenseMatrix<float> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowColMatrixTag,
		   VectorRepresentationTypes::Dense)
{
	linbox_check (A.rowdim () == x.size ());
	linbox_check (A.coldim () == x.size ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_strsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
	return x;
}

template <class Vector1, class Vector2>
DenseMatrix<float> &ger_impl (const UnparametricRing<float> &F, BLASModule &M, float a, const Vector1 &x, const Vector2 &y, DenseMatrix<float> &A,
			      VectorRepresentationTypes::Dense,
			      VectorRepresentationTypes::Dense,
			      MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == x.size ());
	linbox_check (A.coldim () == y.size ());
	cblas_sger (CblasRowMajor, A.rowdim (), A.coldim (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0], &A[0][0], A.disp ());
	return A;
}

template <class Vector1, class Vector2>
Vector2 &gemv_impl (const UnparametricRing<double> &F, BLASModule &M,
		    double a, const DenseMatrix<double> &A, const Vector1 &x, double b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowColMatrixTag,
		    VectorRepresentationTypes::Dense,
		    VectorRepresentationTypes::Dense)
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());
	cblas_dgemv (CblasRowMajor, CblasNoTrans, A.rowdim (), A.coldim (), a, &A[0][0], A.disp (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
	return y;
}

template <class Vector>
Vector &trmv_impl (const UnparametricRing<double> &F, BLASModule &M, const DenseMatrix<double> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowColMatrixTag,
		   VectorRepresentationTypes::Dense)
{
	linbox_check (A.rowdim () == x.size ());
	linbox_check (A.coldim () == x.size ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_dtrmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
	return x;
}

template <class Vector>
Vector &trsv_impl (const UnparametricRing<double> &F, BLASModule &M, const DenseMatrix<double> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowColMatrixTag,
		   VectorRepresentationTypes::Dense)
{
	linbox_check (A.rowdim () == x.size ());
	linbox_check (A.coldim () == x.size ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_dtrsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
	return x;
}

template <class Vector1, class Vector2>
DenseMatrix<double> &ger_impl (const UnparametricRing<double> &F, BLASModule &M, double a, const Vector1 &x, const Vector2 &y, DenseMatrix<double> &A,
			       VectorRepresentationTypes::Dense,
			       VectorRepresentationTypes::Dense,
			       MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == x.size ());
	linbox_check (A.coldim () == y.size ());
	cblas_dger (CblasRowMajor, A.rowdim (), A.coldim (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0], &A[0][0], A.disp ());
	return A;
}

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_CBLAS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
