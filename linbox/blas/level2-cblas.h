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
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

namespace BLAS2
{

template <>
class _gemv<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const UnparametricRing<float> &F, Modules &M,
				   float a, const DenseMatrix<float> &A, const Vector1 &x, float b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic)
		{ return _gemv<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, A, x, b, y, start_idx, end_idx); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const UnparametricRing<float> &F, Modules &M,
				   float a, const DenseMatrix<float> &A, const Vector1 &x, float b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixStorageTypes::Dense,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real)
	{
		linbox_check (A.coldim () == x.size ());
		linbox_check (A.rowdim () == y.size ());
		cblas_sgemv (CblasRowMajor, CblasNoTrans, A.rowdim (), A.coldim (), a, &A[0][0], A.disp (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const UnparametricRing<float> &F,
			    Modules              &M,
			    uint8                 a,
			    const Matrix         &A,
			    const Vector1        &x,
			    uint8                 b,
			    Vector2              &y,
			    size_t                start_idx = 0,
			    size_t                end_idx = (size_t) -1)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    typename Matrix::StorageType (),
				    typename VectorTraits<UnparametricRing<float>, Vector1>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector1>::StorageType (),
				    typename VectorTraits<UnparametricRing<float>, Vector2>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector2>::StorageType ()); }
};

template <>
class _trmv<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const UnparametricRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trmv<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Vector>
	static Vector &trmv_impl (const UnparametricRing<float> &F, Modules &M, const DenseMatrix<float> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		linbox_check (A.rowdim () == x.size ());
		linbox_check (A.coldim () == x.size ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const UnparametricRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trmv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<UnparametricRing<float>, Vector>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector>::StorageType ()); }
};

template <>
class _trsv<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const UnparametricRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trsv<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Vector>
	static Vector &trsv_impl (const UnparametricRing<float> &F, Modules &M, const DenseMatrix<float> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		linbox_check (A.rowdim () == x.size ());
		linbox_check (A.coldim () == x.size ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_strsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const UnparametricRing<float> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trsv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<UnparametricRing<float>, Vector>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<float>, Vector>::StorageType ()); }
};

template <>
class _ger<UnparametricRing<float>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const UnparametricRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 MatrixStorageTypes::Generic)
		{ return _ger<UnparametricRing<float>, BLASModule::Tag::Parent>::op (F, M, a, x, y, A); }

	template <class Modules, class Vector1, class Vector2>
	static DenseMatrix<float> &ger_impl (const UnparametricRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, DenseMatrix<float> &A,
					     VectorRepresentationTypes::Dense,
					     VectorStorageTypes::Real,
					     VectorRepresentationTypes::Dense,
					     VectorStorageTypes::Real,
					     MatrixStorageTypes::Dense)
	{
		linbox_check (A.rowdim () == x.size ());
		linbox_check (A.coldim () == y.size ());
		cblas_sger (CblasRowMajor, A.rowdim (), A.coldim (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0], &A[0][0], A.disp ());
		return A;
	}

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const UnparametricRing<float> &F, Modules &M, float a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<UnparametricRing<float>, Vector1>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<float>, Vector1>::StorageType (),
				   typename VectorTraits<UnparametricRing<float>, Vector2>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<float>, Vector2>::StorageType (),
				   typename Matrix::StorageType ()); }
};

template <>
class _gemv<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const UnparametricRing<double> &F, Modules &M,
				   double a, const DenseMatrix<double> &A, const Vector1 &x, double b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorStorageTypes::Generic)
		{ return _gemv<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, A, x, b, y, start_idx, end_idx); }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const UnparametricRing<double> &F, Modules &M,
				   double a, const DenseMatrix<double> &A, const Vector1 &x, double b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixStorageTypes::Dense,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real,
				   VectorRepresentationTypes::Dense,
				   VectorStorageTypes::Real)
	{
		linbox_check (A.coldim () == x.size ());
		linbox_check (A.rowdim () == y.size ());
		cblas_dgemv (CblasRowMajor, CblasNoTrans, A.rowdim (), A.coldim (), a, &A[0][0], A.disp (), &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0]);
		return y;
	}

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const UnparametricRing<float> &F,
			    Modules              &M,
			    uint8                 a,
			    const Matrix         &A,
			    const Vector1        &x,
			    uint8                 b,
			    Vector2              &y,
			    size_t                start_idx = 0,
			    size_t                end_idx = (size_t) -1)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    typename Matrix::StorageType (),
				    typename VectorTraits<UnparametricRing<double>, Vector1>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector1>::StorageType (),
				    typename VectorTraits<UnparametricRing<double>, Vector2>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector2>::StorageType ()); }
};

template <>
class _trmv<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const UnparametricRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trmv<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Vector>
	static Vector &trmv_impl (const UnparametricRing<double> &F, Modules &M, const DenseMatrix<double> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		linbox_check (A.rowdim () == x.size ());
		linbox_check (A.coldim () == x.size ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrmv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const UnparametricRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trmv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<UnparametricRing<double>, Vector>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector>::StorageType ()); }
};

template <>
class _trsv<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const UnparametricRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Generic,
				  VectorRepresentationTypes::Generic,
				  VectorStorageTypes::Generic)
		{ return _trsv<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, A, x, type, diagIsOne); }

	template <class Modules, class Vector>
	static Vector &trsv_impl (const UnparametricRing<double> &F, Modules &M, const DenseMatrix<double> &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixStorageTypes::Dense,
				  VectorRepresentationTypes::Dense,
				  VectorStorageTypes::Real)
	{
		linbox_check (A.rowdim () == x.size ());
		linbox_check (A.coldim () == x.size ());
		linbox_check (type == UpperTriangular || type == LowerTriangular);
		cblas_dtrsv (CblasRowMajor, (type == UpperTriangular) ? CblasUpper : CblasLower, CblasNoTrans, diagIsOne ? CblasUnit : CblasNonUnit,
			     A.rowdim (), &A[0][0], A.disp (), &x[0], &x[1] - &x[0]);
		return x;
	}
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const UnparametricRing<double> &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trsv_impl (F, M, A, x, type, diagIsOne,
				    typename Matrix::StorageType (),
				    typename VectorTraits<UnparametricRing<double>, Vector>::RepresentationType (),
				    typename VectorTraits<UnparametricRing<double>, Vector>::StorageType ()); }
};

template <>
class _ger<UnparametricRing<double>, BLASModule::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const UnparametricRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 VectorRepresentationTypes::Generic,
				 VectorStorageTypes::Generic,
				 MatrixStorageTypes::Generic)
		{ return _ger<UnparametricRing<double>, BLASModule::Tag::Parent>::op (F, M, a, x, y, A); }

	template <class Modules, class Vector1, class Vector2>
	static DenseMatrix<double> &ger_impl (const UnparametricRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, DenseMatrix<double> &A,
					      VectorRepresentationTypes::Dense,
					      VectorStorageTypes::Real,
					      VectorRepresentationTypes::Dense,
					      VectorStorageTypes::Real,
					      MatrixStorageTypes::Dense)
	{
		linbox_check (A.rowdim () == x.size ());
		linbox_check (A.coldim () == y.size ());
		cblas_dger (CblasRowMajor, A.rowdim (), A.coldim (), a, &x[0], &x[1] - &x[0], &y[0], &y[1] - &y[0], &A[0][0], A.disp ());
		return A;
	}

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const UnparametricRing<double> &F, Modules &M, double a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<UnparametricRing<double>, Vector1>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<double>, Vector1>::StorageType (),
				   typename VectorTraits<UnparametricRing<double>, Vector2>::RepresentationType (),
				   typename VectorTraits<UnparametricRing<double>, Vector2>::StorageType (),
				   typename Matrix::StorageType ()); }
};

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
