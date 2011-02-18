/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/sparse-subvector-hybrid.h
 * Copyright 2011 Bradford Hovinen
 *
 * Evolved from sparse-subvector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SPARSE_SUBVECTOR_HYBRID_H
#define __VECTOR_SPARSE_SUBVECTOR_HYBRID_H

#include <vector>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/sparse-subvector.h"

namespace LinBox
{

template <class Vector>
class HybridSubvectorConstIterator;

// Specialisation of SparseSubvector to hybrid zero-one format

template <class Vector>
class SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>
{
    public:
	typedef std::pair<typename Vector::first_type::value_type, typename Vector::second_type::const_word_iterator::value_type> value_type;
	typedef typename Vector::first_type::size_type             size_type;
	typedef typename Vector::first_type::difference_type       difference_type;

	typedef typename Vector::second_type::Endianness Endianness;
	typedef typename std::iterator_traits<typename Vector::second_type::const_word_iterator>::value_type word_type;

	typedef HybridSubvectorConstIterator<Vector> const_iterator;

	friend class HybridSubvectorConstIterator<Vector>;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, size_t start, size_t end)
		: _start (start), _end (end)
		{ set_start_end (v.first.begin (), v.first.end (), v.second.wordBegin ()); _end_is_end = (_end_idx == v.first.end ()); }
	SparseSubvector (const SparseSubvector &v, size_t start, size_t end)
		: _start (v._start + start), _end (v._start + end)
		{ set_start_end (v._start_idx, v._end_idx, v._start_elt); _end_is_end = v._end_is_end && (_end_idx == v._end_idx); }
	SparseSubvector (const SparseSubvector &v)
		: _start (v._start), _end (v._end), _start_idx (v._start_idx), _end_idx (v._end_idx), _end_marker (v._end_marker), _start_elt (v._start_elt), _end_is_end (v._end_is_end) {}

	~SparseSubvector () {}

	inline const_iterator begin () const;
	inline const_iterator end   () const;

	inline size_t         size  () const { return _end_marker - _start_idx; }
	inline bool           empty () const { return _end_idx == _start_idx; }

    private:
	size_t _start, _end;
	typename Vector::first_type::const_iterator _start_idx, _end_idx, _end_marker;
	typename Vector::second_type::const_word_iterator _start_elt;
	bool _end_is_end;

	void set_start_end (typename Vector::first_type::const_iterator idx_begin,
			    typename Vector::first_type::const_iterator idx_end,
			    typename Vector::second_type::const_word_iterator elt_begin)
	{
		_start_idx = std::lower_bound (idx_begin, idx_end, _start >> WordTraits<word_type>::logof_size);
		_end_idx = std::upper_bound (idx_begin, idx_end, _end >> WordTraits<word_type>::logof_size);
		_start_elt = elt_begin + (_start_idx - idx_begin);

		if (_end_idx != _start_idx && (static_cast<size_t> (*(_end_idx - 1)) << WordTraits<word_type>::logof_size) + (_start & WordTraits<word_type>::pos_mask) >= _end)
			_end_marker = _end_idx - 1;
		else
			_end_marker = _end_idx;
	}
}; // template <class Vector> class SparseSubvector<Vector, HybridZeroOneVectorTag>

// Vector traits for SparseVector wrapper
template <class Vector>
struct VectorTraits<SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag> >
{ 
	typedef SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag> VectorType;
	typedef VectorCategories::HybridZeroOneSequenceVectorTag VectorCategory; 
};

template <class Vector>
struct VectorTraits<const SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag> >
{ 
	typedef SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag> VectorType;
	typedef VectorCategories::HybridZeroOneSequenceVectorTag VectorCategory; 
};

} // namespace LinBox

#include "linbox/vector/sparse-subvector-hybrid.inl"

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_H
