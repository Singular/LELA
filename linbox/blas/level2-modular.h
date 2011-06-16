/* linbox/blas/level2-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for Z/p
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_MODULAR_H
#define __BLAS_LEVEL2_MODULAR_H

#include "linbox/ring/modular.h"
#include "linbox/blas/context.h"
#include "linbox/matrix/matrix-traits.h"

namespace LinBox
{

namespace BLAS2
{

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
			 uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::DenseVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
			 uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::SparseVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Modular<uint8> &F, ZpModule<uint8> &M,
		    uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::GenericVectorTag)
	{ return gemv_col_dense (F, M, a, A, x, b, y, start_idx, end_idx, typename VectorTraits<Modular<uint8>, typename Matrix::Column>::VectorCategory ()); }

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
			 uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::DenseVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
			 uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::SparseVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Modular<uint16> &F, ZpModule<uint16> &M,
		    uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::GenericVectorTag)
	{ return gemv_col_dense (F, M, a, A, x, b, y, start_idx, end_idx, typename VectorTraits<Modular<uint16>, typename Matrix::Column>::VectorCategory ()); }

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
			 uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::DenseVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
			 uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::SparseVectorTag);

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Modular<uint32> &F, ZpModule<uint32> &M,
		    uint32 &a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::GenericVectorTag)
	{ return gemv_col_dense (F, M, a, A, x, b, y, start_idx, end_idx, typename VectorTraits<Modular<uint32>, typename Matrix::Column>::VectorCategory ()); }

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
