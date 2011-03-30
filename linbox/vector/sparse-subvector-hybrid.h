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
	typedef HybridSubvectorConstIterator<Vector> const_iterator;

	typedef typename Vector::index_type index_type;
	typedef typename Vector::word_type word_type;
	typedef std::pair<index_type, word_type> value_type;
	typedef typename Vector::Endianness Endianness;

	friend class HybridSubvectorConstIterator<Vector>;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, size_t start, size_t finish)
		: _start (start), _finish (finish)
		{ set_start_end (v.begin (), v.end ()); _end_is_end = (_end == v.end ()); }
	SparseSubvector (const SparseSubvector &v, size_t start, size_t finish)
		: _start (v._start + start), _finish (v._start + finish)
		{ set_start_end (v._begin, v._end); _end_is_end = v._end_is_end && (_end == v._end); }

	~SparseSubvector () {}

	inline const_iterator begin () const;
	inline const_iterator end   () const;

	inline size_t         size  () const { return _end_marker - _begin; }
	inline bool           empty () const { return _end_marker == _begin; }

    private:
	size_t _start, _finish;
	typename Vector::const_iterator _begin, _end, _end_marker;
	bool _end_is_end;

	void set_start_end (typename Vector::const_iterator begin,
			    typename Vector::const_iterator end)
	{
		_begin = std::lower_bound (begin, end, _start >> WordTraits<word_type>::logof_size, VectorWrapper::CompareSparseEntries ());
		_end = std::upper_bound (begin, end, _finish >> WordTraits<word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

		if (_end != _begin && (static_cast<size_t> ((_end - 1)->first) << WordTraits<word_type>::logof_size) + (_start & WordTraits<word_type>::pos_mask) >= _finish)
			_end_marker = _end - 1;
		else
			_end_marker = _end;
	}
}; // template <class Vector> class SparseSubvector<Vector, HybridZeroOneVectorTag>

} // namespace LinBox

#include "linbox/vector/sparse-subvector-hybrid.inl"

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_H
