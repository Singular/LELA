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
#include "linbox/vector/sparse.h"

namespace LinBox
{

/** A subvector of a sparse vector. It provides an interface which
 * mimics that of the containing vector, but does not allow any
 * operations which modify the vector's contents.
  \ingroup vector
 */

template <class Vector, class Trait> // N.B. default argument in forward-declaration in linbox/vector/sparse.h
class SparseSubvector
{
    public:
	SparseSubvector ();
	SparseSubvector (Vector &v, size_t start, size_t end);
	SparseSubvector (SparseSubvector &v, size_t start, size_t end);
	SparseSubvector (const SparseSubvector &v);

	~SparseSubvector () {}
}; // template <class Vector, Trait> class SparseSubvector

// Specialisation of SparseSubvector to sparse sequence format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseSequenceVectorTag>
	: public ConstSparseVector<typename ShiftedVector<typename Vector::const_index_iterator>::const_iterator, typename Vector::const_element_iterator>
{
    public:
	typedef ConstSparseVector<typename ShiftedVector<typename Vector::const_index_iterator>::const_iterator, typename Vector::const_element_iterator> parent_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		typename Vector::const_index_iterator begin = std::lower_bound (v.index_begin (), v.index_end (), start);
		typename Vector::const_index_iterator end = std::lower_bound (v.index_begin (), v.index_end (), finish);
		
		_idx = ShiftedVector<typename Vector::const_index_iterator> (begin, end, start);
		parent_type::_idx_begin = _idx.begin ();
		parent_type::_idx_end = _idx.end ();
		parent_type::_elt_begin = v.element_begin () + (begin - v.index_begin ());
	}

	SparseSubvector (SparseSubvector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		typename Vector::const_index_iterator begin = std::lower_bound (v._idx._start, v._idx._end, v._idx._shift + start);
		typename Vector::const_index_iterator end = std::lower_bound (v._idx._start, v._idx._end, v._idx._shift + finish);
		
		_idx = ShiftedVector<typename Vector::const_index_iterator> (begin, end, v._idx._shift + start);
		parent_type::_idx_begin = _idx.begin ();
		parent_type::_idx_end = _idx.end ();
		parent_type::_elt_begin = v.element_begin () + (begin - v.index_begin ());
	}

	~SparseSubvector () {}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ _idx = v._idx; this->parent_type::operator = (v); return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<V> &v)
		{ _idx = v._idx; this->parent_type::operator = (v); return *this; }

    private:

	ShiftedVector<typename Vector::const_index_iterator> _idx;

}; // template <class Vector> class SparseSubvector<Vector, SparseParallelVectorTag>

// Specialisation of SparseSubvector to sparse parallel format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseParallelVectorTag>
{
    public:
	typedef const ShiftedVector<typename Vector::first_type> first_type;
	typedef const Subvector<typename Vector::second_type::const_iterator> second_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::first_type::value_type start, typename Vector::first_type::value_type finish)
	{
		typename Vector::const_index_iterator begin = std::lower_bound (v.first.begin (), v.first.end (), start);
		typename Vector::const_index_iterator end = std::lower_bound (v.first.begin (), v.first.end (), finish);
		
		first = first_type (begin, end, start);
		second = second_type (v.second.begin () + (begin - v.first.begin ()), v.second.begin () + (end - v.first.begin ()));
	}

	SparseSubvector (SparseSubvector &v, typename Vector::first_type::value_type start, typename Vector::first_type::value_type finish)
	{
		typename Vector::const_index_iterator begin = std::lower_bound (v.first._v.begin (), v.first._v.end (), start + v.first._shift);
		typename Vector::const_index_iterator end = std::lower_bound (v.first._v.begin (), v.first._v.end (), finish + v.first._shift);
		
		first = first_type (begin, end, start + v.first._shift);
		second = second_type (v.second.begin () + (begin - v.first.begin ()), v.second.begin () + (end - v.first.begin ()));
	}

	~SparseSubvector () {}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ first = v.first; second = v.second; return *this; }

	first_type first;
	second_type second;
}; // template <class Vector> class SparseSubvector<Vector, SparseParallelVectorTag>

// Specialisation of SparseSubvector to sparse zero-one format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseZeroOneVectorTag> : public ShiftedVector<typename Vector::const_iterator>
{
    public:
	typedef ShiftedVector<typename Vector::const_iterator> parent_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: parent_type (std::lower_bound (v.begin (), v.end (), start), std::lower_bound (v.begin (), v.end (), finish), start)
		{}
	SparseSubvector (SparseSubvector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: parent_type (std::lower_bound (v._v.begin (), v._v.end (), start + v._shift), std::lower_bound (v._v.begin (), v._v.end (), finish + v._shift), start + v._shift)
		{}
	SparseSubvector (const SparseSubvector &v)
		: ShiftedVector<typename Vector::const_iterator> (v)
		{}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ parent_type::operator = (v); return *this; }

	~SparseSubvector () {}
}; // template <class Vector> class SparseSubvector<Vector, SparseZeroOneVectorTag>

// Vector traits for SparseVector wrapper
template <class Vector, class Trait>
struct VectorTraits<SparseSubvector<Vector, Trait> >
{ 
	typedef SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

template <class Vector, class Trait>
struct VectorTraits<const SparseSubvector<Vector, Trait> >
{ 
	typedef const SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

} // namespace LinBox

#endif // __VECTOR_SPARSE_SUBVECTOR_H
