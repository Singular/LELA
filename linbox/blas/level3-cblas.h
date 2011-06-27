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
#include "linbox/ring/unparametric.h"
#include "linbox/matrix/traits.h"
#include "linbox/matrix/dense.h"
#include "linbox/blas/level3-ll.h"

namespace LinBox
{

namespace BLAS3
{

template <>
class _gemm<UnparametricRing<float>, BLASModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const UnparametricRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C)
		{ return _gemm<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, A, B, b, C); }

	static DenseMatrix<float> &op (const UnparametricRing<float> &F, BLASModule &M,
				       float a, const DenseMatrix<float> &A, const DenseMatrix<float> &B, float b, DenseMatrix<float> &C)
	{
		linbox_check (A.rowdim () == C.rowdim ());
		linbox_check (B.coldim () == C.coldim ());
		linbox_check (A.coldim () == B.rowdim ());
		cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &A[0][0], A.disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
		return C;
	}
};

template <>
class _trmm<UnparametricRing<float>, BLASModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const UnparametricRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trmm<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	static DenseMatrix<float> &op (const UnparametricRing<float> &F, BLASModule &M, float a, const DenseMatrix<float> &A, DenseMatrix<float> &B,
				       TriangularMatrixType type, bool diagIsOne)
	{
		linbox_check (A.rowdim () == B.rowdim ());
		linbox_check (A.coldim () == B.rowdim ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}
};

template <>
class _trsm<UnparametricRing<float>, BLASModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const UnparametricRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trsm<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	static DenseMatrix<float> &trsm_impl (const UnparametricRing<float> &F, BLASModule &M, float a, const DenseMatrix<float> &A, DenseMatrix<float> &B,
					      TriangularMatrixType type, bool diagIsOne)
	{
		linbox_check (A.rowdim () == B.rowdim ());
		linbox_check (A.coldim () == B.rowdim ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}
};

template <>
class _gemm<UnparametricRing<double>, BLASModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const UnparametricRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C)
		{ return _gemm<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, A, B, b, C); }

	static DenseMatrix<double> &gemm_impl (const UnparametricRing<double> &F, BLASModule &M,
					       double a, const DenseMatrix<double> &A, const DenseMatrix<double> &B, double b, DenseMatrix<double> &C)
	{
		linbox_check (A.rowdim () == C.rowdim ());
		linbox_check (B.coldim () == C.coldim ());
		linbox_check (A.coldim () == B.rowdim ());
		cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &A[0][0], A.disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
		return C;
	}
};

template <>
class _trmm<UnparametricRing<double>, BLASModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const UnparametricRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trmm<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	static DenseMatrix<double> &trmm_impl (const UnparametricRing<double> &F, BLASModule &M, double a, const DenseMatrix<double> &A, DenseMatrix<double> &B,
					       TriangularMatrixType type, bool diagIsOne)
	{
		linbox_check (A.rowdim () == B.rowdim ());
		linbox_check (A.coldim () == B.rowdim ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}
};

template <>
class _trsm<UnparametricRing<double>, BLASModule::Tag>
{
public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const UnparametricRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return _trsm<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	static DenseMatrix<double> &trsm_impl (const UnparametricRing<double> &F, BLASModule &M, double a, const DenseMatrix<double> &A, DenseMatrix<double> &B,
					       TriangularMatrixType type, bool diagIsOne)
	{
		linbox_check (A.rowdim () == B.rowdim ());
		linbox_check (A.coldim () == B.rowdim ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}
};

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
