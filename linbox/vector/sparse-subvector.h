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

#include "linbox/vector/traits.h"
#include "linbox/vector/shifted-vector.h"
#include "linbox/vector/sparse.h"
#include "linbox/vector/subvector.h"

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
	typedef Trait VectorCategory; 

	SparseSubvector ();
	SparseSubvector (Vector &v, size_t start, size_t end);
	SparseSubvector (SparseSubvector &v, size_t start, size_t end);
	SparseSubvector (const SparseSubvector &v);

	~SparseSubvector () {}
}; // template <class Vector, Trait> class SparseSubvector

// Specialisation of SparseSubvector to const vector in sparse format

template <class Vector>
class SparseSubvector<const Vector, VectorCategories::SparseVectorTag>
	: public ConstSparseVector<typename ConstShiftedVector<typename Vector::const_index_iterator>::const_iterator, typename Vector::const_element_iterator>
{
    public:
	typedef VectorCategories::SparseVectorTag VectorCategory; 
	typedef ConstSparseVector<typename ConstShiftedVector<typename Vector::const_index_iterator>::const_iterator, typename Vector::const_element_iterator> parent_type;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		typename Vector::const_index_iterator begin = std::lower_bound (v.index_begin (), v.index_end (), start);
		typename Vector::const_index_iterator end = std::lower_bound (v.index_begin (), v.index_end (), finish);
		
		_idx = ConstShiftedVector<typename Vector::const_index_iterator> (begin, end, start);
		parent_type::_idx_begin = _idx.begin ();
		parent_type::_idx_end = _idx.end ();
		parent_type::_elt_begin = v.element_begin () + (begin - v.index_begin ());
	}

	SparseSubvector (const SparseSubvector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		typename Vector::const_index_iterator begin = std::lower_bound (v._idx._start, v._idx._end, v._idx._shift + start);
		typename Vector::const_index_iterator end = std::lower_bound (v._idx._start, v._idx._end, v._idx._shift + finish);
		
		_idx = ConstShiftedVector<typename Vector::const_index_iterator> (begin, end, v._idx._shift + start);
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

    protected:

	ConstShiftedVector<typename Vector::const_index_iterator> _idx;

}; // template <class Vector> class SparseSubvector<const Vector, SparseVectorTag>

// Specialisation of SparseSubvector to vector in sparse format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseVectorTag>
	: public SparseVector<typename Vector::element_type, ShiftedVector<MutableSubvector<typename Vector::index_vector> >, MutableSubvector<typename Vector::element_vector> >
{
	MutableSubvector<typename Vector::index_vector> _idx_v;

    public:
	typedef VectorCategories::SparseVectorTag VectorCategory; 
	typedef SparseVector<typename Vector::element_type, ShiftedVector<MutableSubvector<typename Vector::index_vector> >, MutableSubvector<typename Vector::element_vector> > parent_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::index_type start, typename Vector::index_type finish)
		: _idx_v (v._idx, std::lower_bound (v.index_begin (), v.index_end (), start), std::lower_bound (v.index_begin (), v.index_end (), finish))
	{
		parent_type::_elt.set_parent (v._elt);
		parent_type::_elt.set_idx_begin (_idx_v.begin () - v.index_begin ());
		parent_type::_elt.set_idx_end (_idx_v.end () - v.index_begin ());
		parent_type::_idx.set_parent (_idx_v);
		parent_type::_idx.set_shift (start);
	}

	SparseSubvector (SparseSubvector &v, typename Vector::index_type start, typename Vector::index_type finish)
		: _idx_v (v._idx_v.parent (), std::lower_bound (v._idx_v.parent ().begin (), v._idx_v.parent ().end (), v._idx.shift () + start),
			  std::lower_bound (v._idx_v.parent ().begin (), v._idx_v.parent ().end (), v._idx.shift () + finish))
	{
		parent_type::_elt.set_parent (v._elt);
		parent_type::_elt.set_idx_begin (_idx_v.begin () - v._idx_v.parent ().begin ());
		parent_type::_elt.set_idx_end (_idx_v.end () - v._idx_v.parent ().begin ());
		parent_type::_idx.set_parent (_idx_v);
		parent_type::_idx.set_shift (start);
	}

	~SparseSubvector () {}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ _idx_v = v._idx_v; this->parent_type::operator = (v); parent_type::_idx.set_parent (_idx_v); return *this; }

    protected:

	SparseSubvector (Vector &v, typename Vector::index_vector::iterator idx_begin, typename Vector::index_vector::iterator idx_end, typename Vector::index_type start)
		: _idx_v (v._idx, idx_begin, idx_end)
	{
		parent_type::_elt.set_parent (v._elt);
		parent_type::_elt.set_idx_begin (_idx_v.begin () - v.index_begin ());
		parent_type::_elt.set_idx_end (_idx_v.end () - v.index_begin ());
		parent_type::_idx.set_parent (_idx_v);
		parent_type::_idx.set_shift (start);
	}

}; // template <class Vector> class SparseSubvector<Vector, SparseVectorTag>

// Specialisation of SparseSubvector to const vector in sparse zero-one format

template <class Vector>
class SparseSubvector<const Vector, VectorCategories::SparseZeroOneVectorTag> : public ConstShiftedVector<typename Vector::const_iterator>
{
    public:
	typedef VectorCategories::SparseZeroOneVectorTag VectorCategory; 
	typedef ConstShiftedVector<typename Vector::const_iterator> parent_type;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: parent_type (std::lower_bound (v.begin (), v.end (), start), std::lower_bound (v.begin (), v.end (), finish), start)
		{}
	SparseSubvector (const SparseSubvector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: parent_type (std::lower_bound (v._v.begin (), v._v.end (), start + v._shift), std::lower_bound (v._v.begin (), v._v.end (), finish + v._shift), start + v._shift)
		{}
	SparseSubvector (const SparseSubvector &v)
		: parent_type (v)
		{}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ parent_type::operator = (v); return *this; }

	~SparseSubvector () {}
}; // template <class Vector> class SparseSubvector<const Vector, SparseZeroOneVectorTag>

// Specialisation of SparseSubvector to vector in sparse zero-one format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::SparseZeroOneVectorTag> : public ShiftedVector<MutableSubvector<Vector> >
{
	MutableSubvector<Vector> _v;

    public:
	typedef VectorCategories::SparseZeroOneVectorTag VectorCategory; 
	typedef ShiftedVector<MutableSubvector<Vector> > parent_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: _v (v, std::lower_bound (v.begin (), v.end (), start), std::lower_bound (v.begin (), v.end (), finish))
		{ parent_type::set_parent (_v); parent_type::set_shift (start); }
	SparseSubvector (SparseSubvector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: _v (v._v.parent (), std::lower_bound (v._v.begin (), v._v.end (), start + v._shift), std::lower_bound (v._v.begin (), v._v.end (), finish + v._shift))
		{ parent_type::set_parent (_v); parent_type::set_shift (start + v._shift); }

	SparseSubvector (const SparseSubvector &v)
		: _v (v._v)
		{ parent_type::set_parent (_v); parent_type::set_shift (v._shift); }

	SparseSubvector &operator = (const SparseSubvector &v)
		{ _v = v._v; parent_type::set_shift (v._shift); parent_type::set_parent (_v); return *this; }

	~SparseSubvector () {}
}; // template <class Vector> class SparseSubvector<Vector, SparseZeroOneVectorTag>

// Specialisation of SparseSubvector to (non-const) hybrid vector

template <class Vector>
class SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>
	: public SparseSubvector<Vector, VectorCategories::SparseVectorTag>
{
    public:
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
	typedef SparseSubvector<Vector, VectorCategories::SparseVectorTag> parent_type;
	typedef typename parent_type::parent_type grandparent_type;
	typedef typename Vector::Endianness Endianness;
	typedef typename Vector::index_type index_type;
	typedef typename Vector::word_type word_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, size_t start, size_t finish)
		: parent_type (v, std::lower_bound (v.index_begin () + 1, v.index_end () - 1, start >> WordTraits<word_type>::logof_size),
			       std::lower_bound (v.index_begin () + 1, v.index_end () - 1, (finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size),
			       start >> WordTraits<word_type>::logof_size)
	{
		linbox_check ((start & WordTraits<word_type>::pos_mask) == 0);
		linbox_check ((finish & WordTraits<word_type>::pos_mask) == 0 || v.empty () || finish > (size_t) (v.back ().first << WordTraits<word_type>::logof_size));
	}

	SparseSubvector (SparseSubvector &v, size_t start, size_t finish)
		: parent_type (v, start >> WordTraits<word_type>::logof_size, (finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size)
	{
		linbox_check ((start & WordTraits<word_type>::pos_mask) == 0);
		linbox_check ((finish & WordTraits<word_type>::pos_mask) == 0 || v.empty () || finish > (size_t) (v.back ().first << WordTraits<word_type>::logof_size));
	}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ this->SparseSubvector<Vector, VectorCategories::SparseVectorTag>::operator = (v); return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<V, VectorCategories::SparseVectorTag> &v)
		{ this->SparseSubvector<Vector, VectorCategories::SparseVectorTag>::operator = (v); return *this; }

}; // template <class Vector> class SparseSubvector<Vector, HybridSubvectorWordAlignedTag>

} // namespace LinBox

#endif // __VECTOR_SPARSE_SUBVECTOR_H
