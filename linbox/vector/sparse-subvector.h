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

/** Empty struct to be used when constructing a sparse subvector of a
 * hybrid vector with word-aligned boundaries. When SparseSubvector is
 * instantiated with this tag, the meaning of the constructor-
 * parameters start and end refers to the block-index (that is,
 * column_index = bitsof (word_type) * block_index), unlike the
 * column-index which is the case with the non-word-aligned
 * version. This version is also much faster.
 */
struct HybridSubvectorWordAlignedTag {};

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

// Specialisation of SparseSubvector to sparse format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseVectorTag>
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
	SparseSubvector &operator = (const SparseSubvector<V, VectorCategories::SparseVectorTag> &v)
		{ _idx = v._idx; this->parent_type::operator = (v); return *this; }

    private:

	ShiftedVector<typename Vector::const_index_iterator> _idx;

}; // template <class Vector> class SparseSubvector<Vector, SparseVectorTag>

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

// Specialisation of SparseSubvector to word-aligned hybrid format

template <class Vector>
class SparseSubvector<Vector, HybridSubvectorWordAlignedTag>
	: public SparseSubvector<Vector, VectorCategories::SparseVectorTag>
{
    public:
	typedef typename Vector::Endianness Endianness;
	typedef typename Vector::index_type index_type;
	typedef typename Vector::word_type word_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, index_type start, index_type finish)
		: SparseSubvector<Vector, VectorCategories::SparseVectorTag> (v, start, finish)
	{}

	SparseSubvector (SparseSubvector &v, index_type start, index_type finish)
		: SparseSubvector<Vector, VectorCategories::SparseVectorTag> (v, start, finish)
	{}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ this->SparseSubvector<Vector, VectorCategories::SparseVectorTag>::operator = (v); return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<V, VectorCategories::SparseVectorTag> &v)
		{ this->SparseSubvector<Vector, VectorCategories::SparseVectorTag>::operator = (v); return *this; }

}; // template <class Vector> class SparseSubvector<Vector, HybridSubvectorWordAlignedTag>

// Vector traits for SparseVector wrapper
template <class Vector, class Trait>
struct DefaultVectorTraits<SparseSubvector<Vector, Trait> >
{ 
	typedef SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

template <class Vector, class Trait>
struct DefaultVectorTraits<const SparseSubvector<Vector, Trait> >
{ 
	typedef const SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

template <class Vector, class Trait>
struct GF2VectorTraits<SparseSubvector<Vector, Trait> >
{ 
	typedef SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

template <class Vector, class Trait>
struct GF2VectorTraits<const SparseSubvector<Vector, Trait> >
{ 
	typedef const SparseSubvector<Vector, Trait> VectorType;
	typedef Trait VectorCategory; 
};

template <class Vector>
struct GF2VectorTraits<SparseSubvector<Vector, HybridSubvectorWordAlignedTag> >
{ 
	typedef SparseSubvector<Vector, HybridSubvectorWordAlignedTag> VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

template <class Vector>
struct GF2VectorTraits<const SparseSubvector<Vector, HybridSubvectorWordAlignedTag> >
{ 
	typedef const SparseSubvector<Vector, HybridSubvectorWordAlignedTag> VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

} // namespace LinBox

#endif // __VECTOR_SPARSE_SUBVECTOR_H
