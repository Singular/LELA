/* linbox/blas/level2-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 2 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_GENERIC_H
#define __BLAS_LEVEL2_GENERIC_H

#include <algorithm>

#include "linbox/blas/context.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

namespace BLAS2
{

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::DenseVectorTag);

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::SparseVectorTag);

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::GenericVectorTag);

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::SparseVectorTag,
		    VectorCategories::GenericVectorTag);

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::DenseVectorTag)
	{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
			    MatrixCategories::ColMatrixTag (),
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::SparseVectorTag)
	{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
			    MatrixCategories::ColMatrixTag (),
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowColMatrixTag,
		    VectorCategories::SparseVectorTag,
		    VectorCategories::DenseVectorTag)
	{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
			    MatrixCategories::ColMatrixTag (),
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowColMatrixTag,
		    VectorCategories::SparseVectorTag,
		    VectorCategories::SparseVectorTag)
	{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
			    MatrixCategories::ColMatrixTag (),
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Matrix, class Vector>
Vector &trmv_impl (const Field &F, GenericModule &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowMatrixTag,
		   VectorCategories::DenseVectorTag);

template <class Field, class Matrix, class Vector>
Vector &trsv_impl (const Field &F, GenericModule &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowMatrixTag,
		   VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseVectorTag,
		  VectorCategories::DenseVectorTag,
		  MatrixCategories::RowMatrixTag);

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseVectorTag,
		  VectorCategories::SparseVectorTag,
		  MatrixCategories::RowMatrixTag);

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseVectorTag,
		  VectorCategories::DenseVectorTag,
		  MatrixCategories::ColMatrixTag);

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseVectorTag,
		  VectorCategories::SparseVectorTag,
		  MatrixCategories::ColMatrixTag);

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseVectorTag,
		  VectorCategories::DenseVectorTag,
		  MatrixCategories::RowColMatrixTag)
	{ return ger_impl (F, M, a, x, y, A,
			   typename VectorTraits<Field, Vector1>::VectorCategory (),
			   typename VectorTraits<Field, Vector2>::VectorCategory (),
			   MatrixCategories::RowMatrixTag ()); }

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseVectorTag,
		  VectorCategories::SparseVectorTag,
		  MatrixCategories::RowColMatrixTag)
	{ return ger_impl (F, M, a, x, y, A,
			   typename VectorTraits<Field, Vector1>::VectorCategory (),
			   typename VectorTraits<Field, Vector2>::VectorCategory (),
			   MatrixCategories::RowMatrixTag ()); }

} // namespace BLAS2

} // namespace LinBox

#include "linbox/blas/level2-generic.tcc"

#endif // __BLAS_LEVEL2_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
