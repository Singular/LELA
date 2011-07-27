/* lela/blas/level2-cblas.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Cblas implementations of level 2 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_CBLAS_H
#define __BLAS_LEVEL2_CBLAS_H

#include <algorithm>

#ifndef __LELA_BLAS_AVAILABLE
#  error "A working installation of BLAS is required to use this header-file."
#endif // __LELA_BLAS_AVAILABLE

#include "lela/cblas.h"

#include "lela/blas/context.h"
#include "lela/ring/type-wrapper.h"
#include "lela/vector/traits.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/dense.h"
#include "lela/blas/level2-ll.h"

namespace LELA
{

namespace BLAS2
{

template <>
class _gemv<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const TypeWrapperRing<float> &F, Modules &M,
				   float a, const Matrix &A, const Vector1 &x, float b, Vector2 &y,
				   MatrixStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic)
		{ return _gemv<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, A, x, b, y); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const TypeWrapperRing<float> &F, Modules &M,
				   float a, const Matrix &A, const Vector1 &x, float b, Vector2 &y,
				   MatrixStorageTypes::Dense,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real)
	{
		lela_check (A.coldim () == x.size ());
		lela_check (A.rowdim () == y.size ());
		cblas_sgemv (CblasRowMajor, CblasNoTrans, A.rowdim (), A.coldim (), a, &A[0][0], A.disp (),
			     &x[0], &x[1] - &x[0], b, &y[0], &y[1] - &y[0]);
		return y;
	}

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const TypeWrapperRing<float> &F, Modules &M,
				   float a, const Matrix &A, const Vector1 &x, float b, Vector2 &y,
				   MatrixStorageTypes::DenseTranspose,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real)
	{
		lela_check (A.coldim () == x.size ());
		lela_check (A.rowdim () == y.size ());
		cblas_sgemv (CblasRowMajor, CblasTrans, A.coldim (), A.rowdim (), a, &(A.parent ())[0][0], A.parent ().disp (),
			     &x[0], &x[1] - &x[0], b, &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const TypeWrapperRing<float> &F,
			    Modules              &M,
			    uint8                 a,
			    const Matrix         &A,
			    const Vector1        &x,
			    uint8                 b,
			    Vector2              &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::StorageType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector1>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector1>::StorageType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector2>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector2>::StorageType ()); }
};

template <>
class _trmv<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trmv<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::DenseTranspose,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &(A.parent ())[0][0], A.parent ().disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trmv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector>::StorageType ()); }
};

template <>
class _trsv<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trsv<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::DenseTranspose,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &(A.parent ())[0][0], A.parent ().disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const TypeWrapperRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trsv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<float>, Vector>::StorageType ()); }
};

template <>
class _ger<TypeWrapperRing<float>, BLASModule<float>::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 MatrixStorageTypes::Generic)
		{ return _ger<TypeWrapperRing<float>, BLASModule<float>::Tag::Parent>::op (F, M, a, x, y, A); }

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense,
				 VectorStorageTypes::Real,
				 VectorRepresentationTypes::Dense,
				 VectorStorageTypes::Real,
				 MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == y.size ());
		cblas_sger (CblasRowMajor, A.rowdim (), A.coldim (), a,
			    &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0], &A[0][0], A.disp ());
		return A;
	}

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense,
				 VectorStorageTypes::Real,
				 VectorRepresentationTypes::Dense,
				 VectorStorageTypes::Real,
				 MatrixStorageTypes::DenseTranspose)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == y.size ());
		cblas_sger (CblasRowMajor, A.coldim (), A.rowdim (), a,
			    &y[0], &y[1] - &y[0], &x[0], &x[1] - &x[0],
			    &(A.parent ())[0][0], A.parent ().disp ());
		return A;
	}

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const TypeWrapperRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<TypeWrapperRing<float>, Vector1>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<float>, Vector1>::StorageType (),
				   typename VectorTraits<TypeWrapperRing<float>, Vector2>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<float>, Vector2>::StorageType (),
				   typename Matrix::StorageType ()); }
};

template <>
class _gemv<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const TypeWrapperRing<double> &F, Modules &M,
				   double a, const Matrix &A, const Vector1 &x, double b, Vector2 &y,
				   MatrixStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic)
		{ return _gemv<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, A, x, b, y); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const TypeWrapperRing<double> &F, Modules &M,
				   double a, const Matrix &A, const Vector1 &x, double b, Vector2 &y,
				   MatrixStorageTypes::Dense,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real)
	{
		lela_check (A.coldim () == x.size ());
		lela_check (A.rowdim () == y.size ());
		cblas_dgemv (CblasRowMajor, CblasNoTrans, A.rowdim (), A.coldim (), a,
			     &A[0][0], A.disp (), &x[0], &x[1] - &x[0], b, &y[0], &y[1] - &y[0]);
		return y;
	}

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const TypeWrapperRing<double> &F, Modules &M,
				   double a, const Matrix &A, const Vector1 &x, double b, Vector2 &y,
				   MatrixStorageTypes::DenseTranspose,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real)
	{
		lela_check (A.coldim () == x.size ());
		lela_check (A.rowdim () == y.size ());
		cblas_dgemv (CblasRowMajor, CblasTrans, A.coldim (), A.rowdim (), a, &(A.parent ())[0][0], A.parent ().disp (),
			     &x[0], &x[1] - &x[0], b, &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const TypeWrapperRing<double> &F,
			    Modules              &M,
			    uint8                 a,
			    const Matrix         &A,
			    const Vector1        &x,
			    uint8                 b,
			    Vector2              &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::StorageType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector1>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector1>::StorageType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector2>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector2>::StorageType ()); }
};

template <>
class _trmv<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trmv<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::DenseTranspose,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &(A.parent ())[0][0], A.parent ().disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trmv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector>::StorageType ()); }
};

template <>
class _trsv<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trsv<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::DenseTranspose,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == x.size ());
		lela_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &(A.parent ())[0][0], A.parent ().disp (), &x[0], &x[1] - &x[0]);
		return x;
	}

public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const TypeWrapperRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trsv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector>::RepresentationType (),
				    typename VectorTraits<TypeWrapperRing<double>, Vector>::StorageType ()); }
};

template <>
class _ger<TypeWrapperRing<double>, BLASModule<double>::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 MatrixStorageTypes::Generic)
		{ return _ger<TypeWrapperRing<double>, BLASModule<double>::Tag::Parent>::op (F, M, a, x, y, A); }

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, Matrix &A,
					      VectorRepresentationTypes::Dense,
					      VectorStorageTypes::Real,
					      VectorRepresentationTypes::Dense,
					      VectorStorageTypes::Real,
					      MatrixStorageTypes::Dense)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == y.size ());
		cblas_dger (CblasRowMajor, A.rowdim (), A.coldim (), a,
			    &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0], &A[0][0], A.disp ());
		return A;
	}

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, Matrix &A,
					      VectorRepresentationTypes::Dense,
					      VectorStorageTypes::Real,
					      VectorRepresentationTypes::Dense,
					      VectorStorageTypes::Real,
					      MatrixStorageTypes::DenseTranspose)
	{
		lela_check (A.rowdim () == x.size ());
		lela_check (A.coldim () == y.size ());
		cblas_dger (CblasRowMajor, A.coldim (), A.rowdim (), a,
			    &y[0], &y[1] - &y[0], &x[0], &x[1] - &x[0], &(A.parent ())[0][0], A.parent ().disp ());
		return A;
	}

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const TypeWrapperRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<TypeWrapperRing<double>, Vector1>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<double>, Vector1>::StorageType (),
				   typename VectorTraits<TypeWrapperRing<double>, Vector2>::RepresentationType (),
				   typename VectorTraits<TypeWrapperRing<double>, Vector2>::StorageType (),
				   typename Matrix::StorageType ()); }
};

} // namespace BLAS2

} // namespace LELA

#endif // __BLAS_LEVEL2_CBLAS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
