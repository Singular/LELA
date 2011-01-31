/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/sparse-subvector.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from sparse-vector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SPARSE_SUBVECTOR_H
#define __VECTOR_SPARSE_SUBVECTOR_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/shifted-vector.h"

namespace LinBox
{

/** A subvector of a sparse vector. It provides an interface which
 * mimics that of the containing vector, but does not allow any
 * operations which modify the vector's contents.
  \ingroup vector
 */

template <class Vector, class Trait = typename VectorTraits<Vector>::VectorCategory>
class SparseSubvector
{
    public:
	SparseSubvector ();
	SparseSubvector (Vector &v, size_t start, size_t end);
	SparseSubvector (SparseSubvector &v, size_t start, size_t end);
	SparseSubvector (const SparseSubvector &v);

	~SparseSubvector () {}
}; // template <class Vector, Trait> class SparseSubvector

// Specialisation of SparseSubvector to sparse parallel format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseParallelVectorTag>
{
    public:
	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::first_type::value_type start, typename Vector::first_type::value_type end)
		: first (v.first, start, end), second (v.second.begin () + (first._start - v.first.begin ()), v.second.begin () + (first._end - v.first.begin ()))
		{}
	SparseSubvector (SparseSubvector &v, typename Vector::first_type::value_type start, typename Vector::first_type::value_type end)
		: first (v.first, start, end), second (v.second.begin () + (first._start - v.first.begin ()), v.second.begin () + (first._end - v.first.begin ()))
		{}
	SparseSubvector (const SparseSubvector &v)
		: first (v.first), second (v.second)
		{}

	~SparseSubvector () {}

	const ShiftedVector<typename Vector::first_type> first;
	const Subvector<typename Vector::second_type::const_iterator> second;
}; // template <class Vector> class SparseSubvector<Vector, SparseParallelVectorTag>

// Specialisation of SparseSubvector to sparse zero-one format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseZeroOneVectorTag> : public ShiftedVector<Vector>
{
    public:
	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::value_type start, typename Vector::value_type end)
		: ShiftedVector<Vector> (v, start, end)
		{}
	SparseSubvector (SparseSubvector &v, typename Vector::value_type start, typename Vector::value_type end)
		: ShiftedVector<Vector> (v, start, end)
		{}
	SparseSubvector (const SparseSubvector &v)
		: ShiftedVector<Vector> (v)
		{}

	~SparseSubvector () {}
}; // template <class Vector> class SparseSubvector<Vector, SparseZeroOneVectorTag>

// Vector traits for SparseVector wrapper
template <class Vector, class Trait>
struct VectorTraits<SparseSubvector<Vector, Trait> >
{ 
	typedef SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

} // namespace LinBox

#endif // __VECTOR_SPARSE_SUBVECTOR_H
