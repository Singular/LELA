/* linbox/blas/level1-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Modular implementations of level 1 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_MODULAR_H
#define __BLAS_LEVEL1_MODULAR_H

#include "linbox/field/modular.h"
#include "linbox/blas/context.h"
#include "linbox/blas/level1-generic.h"

namespace LinBox
{

namespace BLAS1 
{

template <class Vector1, class Vector2>
uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
		 size_t start_idx, size_t end_idx,
		 VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag);

template <class Vector1, class Vector2>
uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
		 size_t start_idx, size_t end_idx,
		 VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Vector1, class Vector2>
uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
		 size_t start_idx, size_t end_idx,
		 VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

template <class Vector1, class Vector2>
uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
		  size_t start_idx, size_t end_idx,
		  VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag);

template <class Vector1, class Vector2>
uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
		  size_t start_idx, size_t end_idx,
		  VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Vector1, class Vector2>
uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
		  size_t start_idx, size_t end_idx,
		  VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

template <class Vector1, class Vector2>
uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
		  size_t start_idx, size_t end_idx,
		  VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag);

template <class Vector1, class Vector2>
uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
		  size_t start_idx, size_t end_idx,
		  VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Vector1, class Vector2>
uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
		  size_t start_idx, size_t end_idx,
		  VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

} // namespace BLAS1

} // namespace LinBox

#include "linbox/blas/level1-modular.tcc"

#endif // __BLAS_LEVEL1_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
