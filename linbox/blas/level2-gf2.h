/* linbox/blas/level2-gf2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for GF2
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_GF2_H
#define __BLAS_LEVEL2_GF2_H

#include <algorithm>
#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/blas/level2-generic.h"
#include "linbox/field/gf2.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/matrix/matrix-traits.h"

namespace LinBox
{

namespace BLAS2
{

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const GF2 &F, GenericModule &M,
		    bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::DenseZeroOneVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const GF2 &F, GenericModule &M,
		    bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::SparseZeroOneVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const GF2 &F, GenericModule &M,
		    bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::HybridZeroOneVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const GF2 &F, GenericModule &M,
		    bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::DenseZeroOneVectorTag,
		    VectorCategories::GenericVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const GF2 &F, GenericModule &M,
		    bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::SparseZeroOneVectorTag,
		    VectorCategories::GenericVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const GF2 &F, GenericModule &M,
		    bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::HybridZeroOneVectorTag,
		    VectorCategories::GenericVectorTag);

template <class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseZeroOneVectorTag,
		  VectorCategories::DenseZeroOneVectorTag,
		  MatrixCategories::RowMatrixTag);

template <class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseZeroOneVectorTag,
		  VectorCategories::SparseZeroOneVectorTag,
		  MatrixCategories::RowMatrixTag);

template <class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::HybridZeroOneVectorTag,
		  VectorCategories::HybridZeroOneVectorTag,
		  MatrixCategories::RowMatrixTag);

template <class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseZeroOneVectorTag,
		  VectorCategories::DenseZeroOneVectorTag,
		  MatrixCategories::ColMatrixTag);

template <class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseZeroOneVectorTag,
		  VectorCategories::SparseZeroOneVectorTag,
		  MatrixCategories::ColMatrixTag);

template <class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::HybridZeroOneVectorTag,
		  VectorCategories::HybridZeroOneVectorTag,
		  MatrixCategories::ColMatrixTag);

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
