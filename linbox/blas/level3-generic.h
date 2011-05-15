/* linbox/blas/level3-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_GENERIC_H
#define __BLAS_LEVEL3_GENERIC_H

#include <algorithm>

#include "linbox/blas/context.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

namespace BLAS3
{

template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

// FIXME: Not yet implemented (need generic way to attach entry to a vector first)
template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag);

// FIXME: Not yet implemented (need generic way to attach entry to a vector first)
template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return copy (F, M, A, B, MatrixCategories::RowMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Field, class Modules, class Matrix>
bool scal_impl (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A, MatrixCategories::RowMatrixTag);

template <class Field, class Modules, class Matrix>
bool scal_impl (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A, MatrixCategories::ColMatrixTag);

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return gemm (F, M, a, A, B, b, C, MatrixCategories::RowMatrixTag (), MatrixCategories::ColMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag, MatrixCategories::RowColMatrixTag)
	{ return gemm (F, M, a, A, B, b, C, MatrixCategories::RowMatrixTag (), MatrixCategories::ColMatrixTag (), MatrixCategories::RowMatrixTag ()); }

template <class Field, class Matrix1, class Matrix2>
Matrix2 &trmm_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix1, class Matrix2>
Matrix2 &trsm_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix1, class Matrix2>
bool equal_impl (const Field &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix1, class Matrix2>
bool equal_impl (const Field &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag);

template <class Field, class Matrix>
bool is_zero_impl (const Field &F, GenericModule &M, const Matrix &A, MatrixCategories::RowMatrixTag);

template <class Field, class Matrix>
bool is_zero_impl (const Field &F, GenericModule &M, const Matrix &A, MatrixCategories::ColMatrixTag);

} // namespace BLAS3

} // namespace LinBox

#include "linbox/blas/level3-generic.tcc"

#endif // __BLAS_LEVEL3_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
