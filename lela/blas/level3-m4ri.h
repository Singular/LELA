/* lela/blas/level3-m4ri.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_M4RI_H
#define __BLAS_LEVEL3_M4RI_H

#ifndef __LELA_HAVE_M4RI
#  error "This header file requires that LELA be configured with libm4ri enabled. Please ensure that libm4ri is properly installed and re-run configure."
#endif

#include <m4ri/m4ri.h>

#include "lela/blas/context.h"
#include "lela/blas/level3-ll.h"
#include "lela/matrix/traits.h"
#include "lela/matrix/m4ri-matrix.h"
#include "lela/matrix/submatrix.h"

namespace LELA
{

namespace BLAS3
{

template <>
class _copy<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const GF2 &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _copy<GF2, M4RIModule::Tag::Parent>::op (F, M, A, B); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &copy_impl (const GF2 &F, Modules &M, const Matrix1 &A, Matrix2 &B,
				   MatrixStorageTypes::M4RI, MatrixStorageTypes::M4RI)
		{ mzd_copy (B._rep, A._rep); return B; }

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const GF2 &F, Modules &M, const Matrix1 &A, Matrix2 &B)
		{ return copy_impl (F, M, A, B,
				    typename Matrix1::StorageType (),
				    typename Matrix2::StorageType ()); }
};

template <>
class _scal<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix>
	static Matrix &scal_impl (const GF2 &F, Modules &M, bool a, Matrix &A, MatrixStorageTypes::Generic)
		{ return _scal<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A); }

	template <class Modules, class Matrix>
	static Matrix &scal_impl (const GF2 &F, Modules &M, bool a, Matrix &A, MatrixStorageTypes::M4RI);

public:
	template <class Modules, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, bool a, Matrix &A)
		{ return scal_impl (F, M, a, A, typename Matrix::StorageType ()); }
};

template <>
class _axpy<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _axpy<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A, B); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B,
				   MatrixStorageTypes::M4RI, MatrixStorageTypes::M4RI);

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B)
		{ return axpy_impl (F, M, a, A, B,
				    typename Matrix1::StorageType (),
				    typename Matrix2::StorageType ()); }
};

template <>
class _gemm<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const GF2 &F, Modules &M, bool a, const Matrix1 &A, const Matrix2 &B, bool b, Matrix3 &C,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _gemm<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const GF2 &F, Modules &M, bool a, const Matrix1 &A, const Matrix2 &B, bool b, Matrix3 &C,
				   MatrixStorageTypes::M4RI, MatrixStorageTypes::M4RI, MatrixStorageTypes::M4RI);

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const GF2 &F, Modules &M, bool a, const Matrix1 &A, const Matrix2 &B, bool b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C,
				    typename Matrix1::StorageType (),
				    typename Matrix2::StorageType (),
				    typename Matrix3::StorageType ()); }
};

template <>
class _trsm<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _trsm<GF2, M4RIModule::Tag::Parent>::op (F, M, a, A, B, type, diagIsOne); }

	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &trsm_impl (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
				   MatrixStorageTypes::M4RI, MatrixStorageTypes::M4RI);

public:
	template <class Modules, class Matrix1, class Matrix2>
	static Matrix2 &op (const GF2 &F, Modules &M, bool a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
		{ return trsm_impl (F, M, a, A, B, type, diagIsOne,
				    typename Matrix1::StorageType (),
				    typename Matrix2::StorageType ()); }
};

template <>
class _permute_rows<GF2, M4RIModule::Tag>
{
	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_rows_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A,
					  MatrixStorageTypes::Generic)
		{ return _permute_rows<GF2, M4RIModule::Tag::Parent>::op (F, M, P_begin, P_end, A); }

	// template <class Modules, class Iterator, class Matrix>
	// static Matrix &permute_rows_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A,
	// 				  MatrixStorageTypes::M4RI);

public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return permute_rows_impl (F, M, P_begin, P_end, A, typename Matrix::StorageType ()); }
};

template <>
class _permute_cols<GF2, M4RIModule::Tag>
{
	template <class Modules, class Iterator, class Matrix>
	static Matrix &permute_cols_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A,
					  MatrixStorageTypes::Generic)
		{ return _permute_cols<GF2, M4RIModule::Tag::Parent>::op (F, M, P_begin, P_end, A); }

	// template <class Modules, class Iterator, class Matrix>
	// static Matrix &permute_cols_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A,
	// 				  MatrixStorageTypes::M4RI);

public:
	template <class Modules, class Iterator, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
		{ return permute_cols_impl (F, M, P_begin, P_end, A, typename Matrix::StorageType ()); }
};

template <>
class _equal<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix1, class Matrix2>
	static bool equal_impl (const GF2 &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
				MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _equal<GF2, M4RIModule::Tag::Parent>::op (F, M, A, B); }

	template <class Modules, class Matrix1, class Matrix2>
	static bool equal_impl (const GF2 &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
				MatrixStorageTypes::M4RI, MatrixStorageTypes::M4RI)
	{
		if (A._rep->offset == 0 && B._rep->offset == 0)
			return mzd_equal (A._rep, B._rep);
		else
			return _equal<GF2, M4RIModule::Tag::Parent>::op (F, M, A, B);
	}

public:
	template <class Modules, class Matrix1, class Matrix2>
	static bool op (const GF2 &F, Modules &M, const Matrix1 &A, const Matrix2 &B)
		{ return equal_impl (F, M, A, B, typename Matrix1::StorageType (), typename Matrix2::StorageType ()); }
};

template <>
class _is_zero<GF2, M4RIModule::Tag>
{
	template <class Modules, class Matrix>
	static bool is_zero_impl (const GF2 &F, Modules &M, const Matrix &A, MatrixStorageTypes::Generic)
		{ return _is_zero<GF2, M4RIModule::Tag::Parent>::op (F, M, A); }

	template <class Modules, class Matrix>
	static bool is_zero_impl (const GF2 &F, Modules &M, const Matrix &A, MatrixStorageTypes::M4RI)
	{
		if (A._rep->offset == 0)
			return mzd_is_zero (A._rep);
		else
			return _is_zero<GF2, M4RIModule::Tag::Parent>::op (F, M, A);
	}

public:
	template <class Modules, class Matrix>
	static bool op (const GF2 &F, Modules &M, const Matrix &A)
		{ return is_zero_impl (F, M, A, typename Matrix::StorageType ()); }
};

} // namespace BLAS3

} // namespace LELA

#endif // __BLAS_LEVEL3_M4RI_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
