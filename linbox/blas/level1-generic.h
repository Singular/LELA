/* linbox/blas/level1-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 1 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_GENERIC_H
#define __BLAS_LEVEL1_GENERIC_H

#include <algorithm>

#include "linbox/blas/context.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

namespace BLAS1 
{

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag)
	{ return _dot (F, M, res, y, x, start_idx, end_idx); }

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

template <class Field, class Vector>
void swap_impl (const Field &F, GenericModule &M, Vector &x, Vector &y, VectorCategories::GenericVectorTag)
	{ std::swap (x, y); }

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (x.size () == y.size ());

	std::copy (x.begin (), x.end (), y.begin ());
	return y;
}

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag);

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag);

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

template <class Field, class Vector>
Vector &scal_impl (const Field &F, GenericModule &M, const typename Field::Element &a, Vector &x, VectorCategories::DenseVectorTag);

template <class Field, class Vector>
Vector &scal_impl (const Field &F, GenericModule &M, const typename Field::Element &a, Vector &x, VectorCategories::SparseVectorTag);

template <class Field, class Modules, class Iterator, class Vector>
Vector &permute_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::DenseVectorTag);

template <class Field, class Modules, class Iterator, class Vector>
Vector &permute_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::SparseVectorTag);

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag);

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag)
	{ return _equal (F, M, y, x); }

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag);

template <class Field, class Vector>
bool is_zero_impl (const Field &F, GenericModule &M, const Vector &x, VectorCategories::DenseVectorTag);

template <class Field, class Vector>
bool is_zero_impl (const Field &F, GenericModule &M, const Vector &x, VectorCategories::SparseVectorTag);

template <class Field, class Vector>
int head_impl (const Field &F, GenericModule &M, typename Field::Element &a, const Vector &x, VectorCategories::DenseVectorTag);

template <class Field, class Vector>
int head_impl (const Field &F, GenericModule &M, typename Field::Element &a, const Vector &x, VectorCategories::SparseVectorTag);

} // namespace BLAS1

} // namespace LinBox

#include "linbox/blas/level1-generic.tcc"

#endif // __BLAS_LEVEL1_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
