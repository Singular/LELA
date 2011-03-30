/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/sparse-subvector-hybrid.inl
 * Copyright 2011 Bradford Hovinen
 *
 * Evolved from sparse-subvector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SPARSE_SUBVECTOR_HYBRID_INL
#define __VECTOR_SPARSE_SUBVECTOR_HYBRID_INL

#include <utility>

#include "linbox/vector/sparse-subvector-hybrid.h"
#include "linbox/vector/bit-subvector.h"

namespace LinBox {

template <class Vector>
class HybridSubvectorConstIterator
{
    public:
	typedef std::forward_iterator_tag iterator_category;
	typedef typename Vector::index_type index_type;
	typedef typename Vector::word_type word_type;
	typedef typename std::pair<index_type, word_type> value_type;
	typedef typename std::pair<index_type, word_type> *pointer;
	typedef typename std::iterator_traits<typename Vector::const_iterator>::difference_type difference_type;
	typedef typename Vector::Endianness Endianness;

	typedef SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag> container_type;

	class const_reference {
	public:
		const_reference () {}
		const_reference (const const_reference &r)
			: first (r.first), second (r.second) {}

		index_type first;
		word_type second;
	};

	typedef const_reference reference;

	HybridSubvectorConstIterator () {}
	HybridSubvectorConstIterator (const container_type &v, const typename Vector::const_iterator &pos)
		: _v (&v), _pos (pos), _ref_first_valid (false), _ref_second_valid (false), _start (true)
		{}

	template <class Vector1>
	HybridSubvectorConstIterator (const HybridSubvectorConstIterator<Vector1> &i)
		: _v (i._v), _pos (i._pos),
		  _ref (i._ref), _ref_first_valid (i._ref_first_valid), _ref_second_valid (i._ref_second_valid), _start (i._start) {}

	template <class Vector1>
	HybridSubvectorConstIterator &operator = (const HybridSubvectorConstIterator<Vector1> &i) {
		_v = reinterpret_cast<const container_type *> (i._v);
		_pos = i._pos;
		_ref.first = i._ref.first;
		_ref.second = i._ref.second;
		_ref_first_valid = i._ref_first_valid;
		_ref_second_valid = i._ref_second_valid;
		_start = i._start;
		return *this;
	}

	HybridSubvectorConstIterator &operator ++ () 
	{
		size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
		size_t col_index = static_cast<size_t> (_pos->first) << WordTraits<word_type>::logof_size;

		if (virtual_index < col_index) {
			++_ref.first;
			_ref_second_valid = false;
		} else {
			++_pos;
			_ref_first_valid = false;
			_ref_second_valid = false;
		}

		return *this;
	}

	HybridSubvectorConstIterator operator ++ (int) 
	{
		HybridSubvectorConstIterator tmp (*this);
		++*this;
		return tmp;
	}

	difference_type operator - (HybridSubvectorConstIterator &i) const 
		{ return _pos - i._pos; }

	const const_reference &operator * ()
		{ update_ref (); return _ref; }

	const const_reference *operator -> ()
		{ update_ref (); return &_ref; }

	template <class Vector1>
	bool operator == (const HybridSubvectorConstIterator<Vector1> &c) const 
		{ return (_pos == c._pos); }

	template <class Vector1>
	bool operator != (const HybridSubvectorConstIterator<Vector1> &c)
	{
		if (c._pos == _v->_end_marker) {
			if (_pos == _v->_end)
				return false;
			else {
				update_ref ();
				size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
				return virtual_index < _v->_finish;
			}
		}
		else
			return _pos != c._pos;
	}

    private:
	template <class Vector1>
	friend class HybridSubvectorConstIterator;

	const container_type *_v;
	typename Vector::const_iterator _pos;
	const_reference _ref;
	bool _ref_first_valid;
	bool _ref_second_valid;
	bool _start;

	inline void update_ref () {
		size_t shift, virtual_index, col_index;

		if (!_ref_first_valid || !_ref_second_valid) {
			shift = _v->_start & WordTraits<word_type>::pos_mask;

			virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
			col_index = static_cast<size_t> (_pos->first) << WordTraits<word_type>::logof_size;
		}

		if (!_ref_first_valid) {
			if (_start) {
				if (col_index < _v->_start)
					_ref.first = 0;
				else
					_ref.first = (col_index - _v->_start) >> WordTraits<word_type>::logof_size;

				_start = false;
			}
			else if (virtual_index + WordTraits<word_type>::bits < col_index)
				_ref.first = (col_index - _v->_start) >> WordTraits<word_type>::logof_size;
			else
				++_ref.first;

			virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;

			_ref_first_valid = true;
		}

		if (!_ref_second_valid) {
			if (shift == 0)
				_ref.second = _pos->second;
			else if (virtual_index < col_index)
				_ref.second = Endianness::shift_right (_pos->second, WordTraits<word_type>::bits - shift);
			else if ((_pos + 1 == _v->_end && _v->_end_is_end) || (_pos + 1)->first > _pos->first + 1)
				_ref.second = Endianness::shift_left (_pos->second, shift);
			else
				_ref.second = Endianness::shift_left (_pos->second, shift) | Endianness::shift_right ((_pos + 1)->second, WordTraits<word_type>::bits - shift);

			if ((static_cast<size_t> (_ref.first) + 1) << WordTraits<word_type>::logof_size > _v->_finish - _v->_start)
				_ref.second &= Endianness::mask_left ((_v->_finish - _v->_start) & WordTraits<word_type>::pos_mask);

			_ref_second_valid = true;
		}
	}
};

template <class Vector>
inline typename SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator
SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::begin () const
{
	return const_iterator (*this, _begin);
}

template <class Vector>
inline typename SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator
SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::end () const
{
	return const_iterator (*this, _end_marker);
}

} // namespace LinBox

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_INL
