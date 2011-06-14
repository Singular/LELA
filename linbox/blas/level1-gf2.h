/* linbox/blas/level1-gf2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for GF2
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_GF2_H
#define __BLAS_LEVEL1_GF2_H

#include <algorithm>
#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/blas/level1-generic.h"
#include "linbox/field/gf2.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

namespace BLAS1 
{

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
	{ return _dot (F, M, res, y, x, start_idx, end_idx); }

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
	{ return _dot (F, M, res, y, x, start_idx, end_idx); }

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
	{ std::copy (x.word_begin (), x.word_end (), y.word_begin ()); y.back_word () = x.back_word (); return y; }

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
	{ y.assign (x.begin (), x.end ()); return y; }

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
	{ y.assign (x.begin (), x.end ()); return y; }

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

template <class Vector>
Vector &scal_impl (const GF2 &F, GenericModule &M, bool a, Vector &x, VectorCategories::DenseZeroOneVectorTag)
	{ if (!a) { std::fill (x.word_begin (), x.word_end (), 0); x.back_word () = 0; } return x; }

template <class Vector>
Vector &scal_impl (const GF2 &F, GenericModule &M, bool a, Vector &x, VectorCategories::SparseZeroOneVectorTag)
	{ if (!a) x.clear (); return x; }

template <class Vector>
Vector &scal_impl (const GF2 &F, GenericModule &M, bool a, Vector &x, VectorCategories::HybridZeroOneVectorTag)
	{ if (!a) x.clear (); return x; }

template <class Modules, class Iterator, class Vector>
Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::SparseZeroOneVectorTag);

template <class Modules, class Iterator, class Vector>
Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::HybridZeroOneVectorTag);

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag);

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
	{ return _equal (F, M, y, x); }

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
	{ return _equal (F, M, y, x); }

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::HybridZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag);

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
	{ return _equal (F, M, y, x); }

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag);

template <class Vector>
bool is_zero_impl (const GF2 &F, GenericModule &M, const Vector &x, VectorCategories::DenseZeroOneVectorTag);

template <class Vector>
bool is_zero_impl (const GF2 &F, GenericModule &M, const Vector &x, VectorCategories::SparseZeroOneVectorTag)
	{ return x.empty (); }

template <class Vector>
bool is_zero_impl (const GF2 &F, GenericModule &M, const Vector &x, VectorCategories::HybridZeroOneVectorTag)
	{ return x.empty (); }

template <class reference, class Vector>
int head_impl (const GF2 &F, GenericModule &M, reference &a, const Vector &x, VectorCategories::DenseZeroOneVectorTag);

template <class reference, class Vector>
int head_impl (const GF2 &F, GenericModule &M, reference &a, const Vector &x, VectorCategories::SparseZeroOneVectorTag)
	{ if (x.empty ()) return -1; a = true; return x.front (); }

template <class reference, class Vector>
int head_impl (const GF2 &F, GenericModule &M, reference &a, const Vector &x, VectorCategories::HybridZeroOneVectorTag);

template <class Modules, class Vector>
std::istream &read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &v, VectorCategories::DenseZeroOneVectorTag);

template <class Modules, class Vector>
std::istream &read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &v, VectorCategories::SparseZeroOneVectorTag);

template <class Modules, class Vector>
std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &v, VectorCategories::DenseZeroOneVectorTag);

template <class Modules, class Vector>
std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &v, VectorCategories::SparseZeroOneVectorTag);

template <class Modules, class Vector>
std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &v, VectorCategories::HybridZeroOneVectorTag);

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
