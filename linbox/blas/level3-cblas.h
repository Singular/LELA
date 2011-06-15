/* linbox/blas/level3-cblas.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Cblas implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_CBLAS_H
#define __BLAS_LEVEL3_CBLAS_H

#include <algorithm>

#ifndef __LINBOX_BLAS_AVAILABLE
#  error "A working installation of BLAS is required to use this header-file."
#endif // __LINBOX_BLAS_AVAILABLE

#include <cblas.h>

#include "linbox/blas/context.h"
#include "linbox/matrix/matrix-traits.h"
#include "linbox/matrix/dense.h"

namespace LinBox
{

namespace BLAS3
{

DenseMatrix<float> &gemm_impl (const Modular<float> &F, BLASModule &M,
			       float a, const DenseMatrix<float> &A, const DenseMatrix<float> &B, float b, DenseMatrix<float> &C,
			       MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());
	linbox_check (A.coldim () == B.rowdim ());
	cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
		     a, &A[0][0], A.disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
	return C;
}

DenseMatrix<float> &trmm_impl (const Modular<float> &F, BLASModule &M, float a, const DenseMatrix<float> &A, DenseMatrix<float> &B,
			       TriangularMatrixType type, bool diagIsOne,
			       MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_strmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
	return B;
}

DenseMatrix<float> &trsm_impl (const Modular<float> &F, BLASModule &M, float a, const DenseMatrix<float> &A, DenseMatrix<float> &B,
			       TriangularMatrixType type, bool diagIsOne,
			       MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_strsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
	return B;
}

DenseMatrix<double> &gemm_impl (const Modular<double> &F, BLASModule &M,
				double a, const DenseMatrix<double> &A, const DenseMatrix<double> &B, double b, DenseMatrix<double> &C,
				MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());
	linbox_check (A.coldim () == B.rowdim ());
	cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
		     a, &A[0][0], A.disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
	return C;
}

DenseMatrix<double> &trmm_impl (const Modular<double> &F, BLASModule &M, double a, const DenseMatrix<double> &A, DenseMatrix<double> &B,
				TriangularMatrixType type, bool diagIsOne,
				MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_dtrmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
	return B;
}

DenseMatrix<double> &trsm_impl (const Modular<double> &F, BLASModule &M, double a, const DenseMatrix<double> &A, DenseMatrix<double> &B,
				TriangularMatrixType type, bool diagIsOne,
				MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (type == UpperTriangular || type == LowerTriangular);
	cblas_dtrsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
		     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
	return B;
}

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_CBLAS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
