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
	typedef typename std::iterator_traits<typename Vector::first_type::const_iterator>::value_type index_type;
	typedef typename std::iterator_traits<typename Vector::second_type::const_word_iterator>::value_type word_type;
	typedef typename std::pair<index_type, word_type> value_type;
	typedef typename Vector::first_type::difference_type difference_type;
	typedef typename Vector::second_type::const_iterator::Endianness Endianness;

	typedef SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag> container_type;

	class const_reference {
	public:
		const_reference () {}
		const_reference (const const_reference &r)
			: first (r.first), second (r.second) {}

		index_type first;
		word_type second;
	};

	HybridSubvectorConstIterator () {}
	HybridSubvectorConstIterator (const container_type &v,
				      const typename Vector::first_type::const_iterator &index_pos,
				      const typename Vector::second_type::const_word_iterator &value_pos)
		: _v (&v), _index_pos (index_pos), _value_pos (value_pos), _ref_first_valid (false), _ref_second_valid (false), _start (true)
		{}

	template <class Vector1>
	HybridSubvectorConstIterator (const HybridSubvectorConstIterator<Vector1> &i)
		: _v (i._v), _index_pos (i._index_pos), _value_pos (i._value_pos),
		  _ref (i._ref), _ref_first_valid (i._ref_first_valid), _ref_second_valid (i._ref_second_valid), _start (i._start) {}

	template <class Vector1>
	HybridSubvectorConstIterator &operator = (const HybridSubvectorConstIterator<Vector1> &i) {
		_v = reinterpret_cast<const container_type *> (i._v);
		_index_pos = i._index_pos;
		_value_pos = i._value_pos;
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
		size_t col_index = static_cast<size_t> (*_index_pos) << WordTraits<word_type>::logof_size;

		if (virtual_index < col_index) {
			++_ref.first;
			_ref_second_valid = false;
		} else {
			++_index_pos;
			++_value_pos;
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
		{ return _index_pos - i._index_pos; }

	const const_reference &operator * ()
		{ update_ref (); return _ref; }

	const const_reference *operator -> ()
		{ update_ref (); return &_ref; }

	template <class Vector1>
	bool operator == (const HybridSubvectorConstIterator<Vector1> &c) const 
		{ return (_index_pos == c._index_pos); }

	template <class Vector1>
	bool operator != (const HybridSubvectorConstIterator<Vector1> &c)
	{
		if (c._index_pos == _v->_end_marker) {
			if (_index_pos == _v->_end_idx)
				return false;
			else {
				update_ref ();
				size_t virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
				return virtual_index < _v->_end;
			}
		}
		else
			return _index_pos != c._index_pos;
	}

    private:
	template <class Vector1>
	friend class HybridSubvectorConstIterator;

	const container_type *_v;
	typename Vector::first_type::const_iterator _index_pos;
	typename Vector::second_type::const_word_iterator _value_pos;
	const_reference _ref;
	bool _ref_first_valid;
	bool _ref_second_valid;
	bool _start;

	inline void update_ref () {
		size_t shift, virtual_index, col_index;

		if (!_ref_first_valid || !_ref_second_valid) {
			shift = _v->_start & WordTraits<word_type>::pos_mask;

			virtual_index = (static_cast<size_t> (_ref.first) << WordTraits<word_type>::logof_size) + _v->_start;
			col_index = static_cast<size_t> (*_index_pos) << WordTraits<word_type>::logof_size;
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
				_ref.second = *_value_pos;
			else if (virtual_index < col_index)
				_ref.second = Endianness::shift_right (*_value_pos, WordTraits<word_type>::bits - shift);
			else if ((_index_pos + 1 == _v->_end_idx && _v->_end_is_end) || _index_pos[1] > *_index_pos + 1)
				_ref.second = Endianness::shift_left (*_value_pos, shift);
			else
				_ref.second = Endianness::shift_left (*_value_pos, shift) | Endianness::shift_right (_value_pos[1], WordTraits<word_type>::bits - shift);

			if ((static_cast<size_t> (_ref.first) + 1) << WordTraits<word_type>::logof_size > _v->_end - _v->_start)
				_ref.second &= Endianness::mask_left ((_v->_end - _v->_start) & WordTraits<word_type>::pos_mask);

			_ref_second_valid = true;
		}
	}
};

template <class Vector>
inline typename SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator
SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::begin () const
{
	return const_iterator (*this, _start_idx, _start_elt);
}

template <class Vector>
inline typename SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::const_iterator
SparseSubvector<Vector, VectorCategories::HybridZeroOneVectorTag>::end () const
{
	return const_iterator (*this, _end_marker, typename Vector::second_type::const_word_iterator ());
}

} // namespace LinBox

namespace std {
	template <class Vector>
	struct iterator_traits<LinBox::HybridSubvectorConstIterator<Vector> >
	{
		typedef forward_iterator_tag iterator_category;
		typedef typename LinBox::HybridSubvectorConstIterator<Vector>::const_reference reference;
		typedef typename LinBox::HybridSubvectorConstIterator<Vector>::const_reference *pointer;
		typedef typename LinBox::HybridSubvectorConstIterator<Vector>::value_type value_type;
		typedef typename LinBox::HybridSubvectorConstIterator<Vector>::difference_type difference_type;
	};
}

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_INL
