/* lela/blas/level3-cblas.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Cblas implementations of level 3 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_CBLAS_H
#define __BLAS_LEVEL3_CBLAS_H

#include <algorithm>

#ifndef __LELA_BLAS_AVAILABLE
#  error "A working installation of BLAS is required to use this header-file."
#endif // __LELA_BLAS_AVAILABLE

#include "lela/cblas.h"

#include "lela/blas/context.h"
#include "lela/ring/type-wrapper.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/dense.h"
#include "lela/blas/level3-ll.h"

namespace LELA
{

namespace BLAS3
{

template <>
class _gemm<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _gemm<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (A.coldim () <= A.disp ());
		lela_check (B.coldim () <= B.disp ());
		lela_check (C.coldim () <= C.disp ());
		cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &A[0][0], A.disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
		return C;
	}

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (A.rowdim () <= A.parent ().disp ());
		lela_check (B.coldim () <= B.disp ());
		lela_check (C.coldim () <= C.disp ());
		cblas_sgemm (CblasRowMajor, CblasTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &(A.parent ())[0][0], A.parent ().disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
		return C;
	}

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (A.coldim () <= A.disp ());
		lela_check (B.rowdim () <= B.parent ().disp ());
		lela_check (C.coldim () <= C.disp ());
		cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &A[0][0], A.disp (), &(B.parent ())[0][0], B.parent ().disp (), b, &C[0][0], C.disp ());
		return C;
	}

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (A.rowdim () <= A.parent ().disp ());
		lela_check (B.rowdim () <= B.parent ().disp ());
		lela_check (C.coldim () <= C.disp ());
		cblas_sgemm (CblasRowMajor, CblasTrans, CblasTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &(A.parent ())[0][0], A.parent ().disp (), &(B.parent ())[0][0], B.parent ().disp (), b, &C[0][0], C.disp ());
		return C;
	}

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C, typename Matrix1::StorageType (), typename Matrix2::StorageType (), typename Matrix3::StorageType ()); }
};

template <>
class _trmm<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _trmm<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &(A.parent ())[0][0], A.parent ().disp (), &B[0][0], B.disp ());
		return B;
	}

public:

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix1 &op (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B,
			    TriangularMatrixType type, bool diagIsOne)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne, typename Matrix1::StorageType (), typename Matrix2::StorageType ()); }
};

template <>
class _trsm<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _trsm<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &(A.parent ())[0][0], A.parent ().disp (), &B[0][0], B.disp ());
		return B;
	}

public:

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix1 &op (const TypeWrapperRing<float> &F, Modules &M, float a, const Matrix1 &A, Matrix2 &B,
			    TriangularMatrixType type, bool diagIsOne)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne, typename Matrix1::StorageType (), typename Matrix2::StorageType ()); }
};

template <>
class _gemm<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _gemm<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &A[0][0], A.disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
		return C;
	}

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		cblas_dgemm (CblasRowMajor, CblasTrans, CblasNoTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &(A.parent ())[0][0], A.parent ().disp (), &B[0][0], B.disp (), b, &C[0][0], C.disp ());
		return C;
	}

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &A[0][0], A.disp (), &(B.parent ())[0][0], B.parent ().disp (), b, &C[0][0], C.disp ());
		return C;
	}

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == C.rowdim ());
		lela_check (B.coldim () == C.coldim ());
		lela_check (A.coldim () == B.rowdim ());
		cblas_dgemm (CblasRowMajor, CblasTrans, CblasTrans, C.rowdim (), C.coldim (), A.coldim (),
			     a, &(A.parent ())[0][0], A.parent ().disp (), &(B.parent ())[0][0], B.parent ().disp (), b, &C[0][0], C.disp ());
		return C;
	}

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C, typename Matrix1::StorageType (), typename Matrix2::StorageType (), typename Matrix3::StorageType ()); }
};

template <>
class _trmm<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _trmm<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trmm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrmm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &(A.parent ())[0][0], A.parent ().disp (), &B[0][0], B.disp ());
		return B;
	}

public:

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix1 &op (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B,
			    TriangularMatrixType type, bool diagIsOne)
		{ return trmm_impl (F, M, a, A, B, type, diagIsOne, typename Matrix1::StorageType (), typename Matrix2::StorageType ()); }
};

template <>
class _trsm<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _trsm<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &A[0][0], A.disp (), &B[0][0], B.disp ());
		return B;
	}

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::DenseTranspose, MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == B.rowdim ());
		lela_check (A.coldim () == B.rowdim ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrsm (CblasRowMajor, CblasLeft, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     B.rowdim (), B.coldim (), a, &(A.parent ())[0][0], A.parent ().disp (), &B[0][0], B.disp ());
		return B;
	}

public:

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix1 &op (const TypeWrapperRing<double> &F, Modules &M, double a, const Matrix1 &A, Matrix2 &B,
			    TriangularMatrixType type, bool diagIsOne)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne, typename Matrix1::StorageType (), typename Matrix2::StorageType ()); }
};

} // namespace BLAS3

} // namespace LELA

#endif // __BLAS_LEVEL3_CBLAS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
