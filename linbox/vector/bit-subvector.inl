/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/bit-subvector.inl
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from bit-vector.inl
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_BIT_SUBVECTOR_INL
#define __VECTOR_BIT_SUBVECTOR_INL

#include "linbox/vector/bit-subvector.h"

namespace std 
{
	template <class Iterator, class ConstIterator>
	struct iterator_traits<LinBox::BitSubvectorWordIterator<Iterator, ConstIterator> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::BitSubvector<Iterator, ConstIterator>::word_reference reference;
		typedef typename LinBox::BitSubvector<Iterator, ConstIterator>::word_reference *pointer;
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type value_type;
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::difference_type difference_type;
	};

	template <class Iterator, class ConstIterator>
	struct iterator_traits<LinBox::BitSubvectorConstWordIterator<Iterator, ConstIterator> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename std::iterator_traits<typename ConstIterator::const_word_iterator>::value_type &reference;
		typedef typename std::iterator_traits<typename ConstIterator::const_word_iterator>::value_type *pointer;
		typedef typename std::iterator_traits<typename ConstIterator::const_word_iterator>::value_type value_type;
		typedef typename std::iterator_traits<typename ConstIterator::const_word_iterator>::difference_type difference_type;
	};
}

namespace LinBox {

template <class Iterator, class ConstIterator>
class BitSubvector<Iterator, ConstIterator>::word_reference {
    public:
	typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type word;
	typedef typename Iterator::Endianness Endianness;

	word_reference () {}

	word_reference (typename Iterator::word_iterator position, uint8 shift)
		: _pos (position), _shift (shift)
	{}

	~word_reference () {}

	word_reference &operator = (word_reference &a) {
		*this = (word) a;
		return a;
	}

	word_reference &operator = (word v) {
		v &= _mask;

		*_pos = (*_pos & (Endianness::shift_right (~_mask, _shift) | Endianness::mask_left (_shift))) | Endianness::shift_right (v, _shift);

		if (_shift != 0 && !_just_this_word)
			_pos[1] = (_pos[1] & (~Endianness::shift_left (_mask, WordTraits<word>::bits - _shift) | Endianness::mask_right (_shift))) |
				Endianness::shift_left (v, WordTraits<word>::bits - _shift);

		return *this;
	}

	word_reference &operator &= (word_reference &a)
		{ return *this = (word) *this & (word) a; }

	word_reference &operator &= (word v) 
		{ return *this = (word) *this & v; }

	word_reference &operator |= (word_reference &a) 
		{ return *this = (word) *this | (word) a; }

	word_reference &operator |= (word v) 
		{ return *this = (word) *this | v; }

	word_reference &operator ^= (word_reference &a) 
		{ return *this = ((word) *this) ^ ((word) a); }

	word_reference &operator ^= (word v) 
		{ return *this = ((word) *this) ^ v; }

	operator word (void) const
		{ return (Endianness::shift_left (*_pos, _shift) | ((_just_this_word || _shift == 0) ? 0 : Endianness::shift_right (_pos[1], WordTraits<word>::bits - _shift))) & _mask; }

    private:
	friend class BitSubvectorWordIterator<Iterator, ConstIterator>;
	friend class BitSubvectorConstWordIterator<Iterator, ConstIterator>;
	friend class const_word_reference;

	typename Iterator::word_iterator _pos;
	uint8 _shift;

	word _mask;
	bool _just_this_word;
};

template <class Iterator, class ConstIterator>
class BitSubvectorWordIterator {
    public:
	typedef typename std::iterator_traits<BitSubvectorWordIterator<Iterator, ConstIterator> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<BitSubvectorWordIterator<Iterator, ConstIterator> >::reference reference;
	typedef typename std::iterator_traits<BitSubvectorWordIterator<Iterator, ConstIterator> >::pointer pointer;
	typedef typename std::iterator_traits<BitSubvectorWordIterator<Iterator, ConstIterator> >::value_type value_type;
	typedef typename std::iterator_traits<BitSubvectorWordIterator<Iterator, ConstIterator> >::difference_type difference_type;
	typedef typename Iterator::Endianness Endianness;

	BitSubvectorWordIterator () {}
	BitSubvectorWordIterator (BitSubvector<Iterator, ConstIterator> &v, typename Iterator::word_iterator pos)
		: _v (&v), _ref (pos, v._begin.pos ()) { init_mask (); }
	BitSubvectorWordIterator (const BitSubvectorWordIterator<Iterator, ConstIterator> &i)
		: _v (i._v), _ref (i._ref._pos, i._ref._shift) { init_mask (); }

	BitSubvectorWordIterator<Iterator, ConstIterator> &operator = (const BitSubvectorWordIterator<Iterator, ConstIterator> &i) {
		_v = i._v;
		_ref._pos = i._ref._pos;
		_ref._shift = i._ref._shift;
		_ref._mask = i._ref._mask;
		_ref._just_this_word = i._ref._just_this_word;
		return *this;
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> &operator ++ () 
	{
		++_ref._pos;
		update_mask ();
		return *this;
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> operator ++ (int) 
	{
		BitSubvectorWordIterator<Iterator, ConstIterator> tmp (*this);
		++*this;
		return tmp;
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> operator + (difference_type i) const
	{
		return BitSubvectorWordIterator (*const_cast<BitSubvector<Iterator, ConstIterator> *> (_v), _ref._pos + i);
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> &operator += (difference_type i) 
	{
		_ref._pos += i;
		init_mask ();
		return *this;
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> &operator -- () 
	{
		--_ref._pos;
		init_mask ();
		return *this;
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> operator -- (int) 
	{
		BitSubvectorWordIterator<Iterator, ConstIterator> tmp (*this);
		--*this;
		return tmp;
	}

	BitSubvectorWordIterator<Iterator, ConstIterator> operator - (difference_type i) const
		{ return *this + -i; }

	BitSubvectorWordIterator<Iterator, ConstIterator> &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (BitSubvectorWordIterator<Iterator, ConstIterator> &i) const 
		{ return _ref._pos - i._ref._pos; }

	typename BitSubvector<Iterator, ConstIterator>::word_reference operator [] (long i) 
		{ return *(*this + i); }

	typename BitSubvector<Iterator, ConstIterator>::word_reference operator * () 
		{ return _ref; }

	typename BitSubvector<Iterator, ConstIterator>::word_reference *operator -> () 
		{ return &_ref; }

	value_type operator * () const 
		{ return _ref; }

	bool operator == (const BitSubvectorWordIterator<Iterator, ConstIterator> &c) const 
		{ return (_ref._pos == c._ref._pos); }

	bool operator != (const BitSubvectorWordIterator<Iterator, ConstIterator> &c) const 
		{ return (_ref._pos != c._ref._pos); }

    private:
	friend class BitSubvectorConstWordIterator<Iterator, ConstIterator>;

	BitSubvector<Iterator, ConstIterator> *_v;
	typename BitSubvector<Iterator, ConstIterator>::word_reference _ref;

	inline void init_mask () {
		_ref._mask = (value_type) 0 - (value_type) 1;
		update_mask ();
	}

	inline void update_mask () {
		if (((_ref._pos - _v->_begin.word () + 1) << WordTraits<value_type>::logof_size) > _v->_end - _v->_begin)
			_ref._mask = Endianness::mask_left ((_v->_end - _v->_begin) & WordTraits<value_type>::pos_mask);

		if (_ref._pos + 1 == _v->_end_marker)
			_ref._just_this_word = true;
		else
			_ref._just_this_word = false;
	}
};

template <class Iterator, class ConstIterator>
class BitSubvectorConstWordIterator {
    public:
	typedef typename std::iterator_traits<BitSubvectorConstWordIterator<Iterator, ConstIterator> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<BitSubvectorConstWordIterator<Iterator, ConstIterator> >::reference reference;
	typedef typename std::iterator_traits<BitSubvectorConstWordIterator<Iterator, ConstIterator> >::pointer pointer;
	typedef typename std::iterator_traits<BitSubvectorConstWordIterator<Iterator, ConstIterator> >::value_type value_type;
	typedef typename std::iterator_traits<BitSubvectorConstWordIterator<Iterator, ConstIterator> >::difference_type difference_type;
	typedef typename Iterator::Endianness Endianness;

	BitSubvectorConstWordIterator () {}
	BitSubvectorConstWordIterator (const BitSubvector<Iterator, ConstIterator> &v, typename Iterator::const_word_iterator pos)
		: _v (&v), _pos (pos), _ref_valid (false) {}
	BitSubvectorConstWordIterator (const BitSubvectorConstWordIterator<Iterator, ConstIterator>  &i)
		: _v (i._v), _pos (i._pos), _word (i._word), _ref_valid (i._ref_valid) {}
	BitSubvectorConstWordIterator (const BitSubvectorWordIterator<Iterator, ConstIterator>  &i)
		: _v (i._v), _pos (i._ref._pos), _ref_valid (false) {}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> &operator = (const BitSubvectorConstWordIterator<Iterator, ConstIterator>  &i) {
		_v = i._v;
		_pos = i._pos;
		_word = i._word;
		_ref_valid = i._ref_valid;
		return *this;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> &operator = (const BitSubvectorWordIterator<Iterator, ConstIterator>  &i) {
		_v = i._v;
		_pos = i._ref._pos;
		_ref_valid = false;
		return *this;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> &operator ++ () 
	{
		++_pos;
		_ref_valid = false;
		return *this;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> operator ++ (int) 
	{
		BitSubvectorConstWordIterator tmp (*this);
		++*this;
		return tmp;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> operator + (difference_type i) const
	{
		return BitSubvectorConstWordIterator<Iterator, ConstIterator> (*_v, _pos + i);
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> &operator += (difference_type i) 
	{
		_pos += i;
		_ref_valid = false;
		return *this;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> &operator -- () 
	{
		--_pos;
		_ref_valid = false;
		return *this;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> operator -- (int) 
	{
		BitSubvectorConstWordIterator<Iterator, ConstIterator> tmp (*this);
		--*this;
		return tmp;
	}

	BitSubvectorConstWordIterator<Iterator, ConstIterator> operator - (difference_type i) const
		{ return *this + -i; }

	BitSubvectorConstWordIterator<Iterator, ConstIterator> &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (BitSubvectorConstWordIterator<Iterator, ConstIterator> &i) const 
		{ return _pos - i._pos; }

	value_type operator [] (long i) 
		{ return *(*this + i); }

	value_type operator * ()
		{ update_ref (); return _word; }

	bool operator == (const BitSubvectorWordIterator<Iterator, ConstIterator> &c) const 
		{ return (_pos == c._ref._pos); }

	bool operator == (const BitSubvectorConstWordIterator<Iterator, ConstIterator> &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const BitSubvectorWordIterator<Iterator, ConstIterator> &c) const 
		{ return (_pos != c._ref._pos); }

	bool operator != (const BitSubvectorConstWordIterator<Iterator, ConstIterator> &c) const 
		{ return (_pos != c._pos); }

    private:
	friend class BitSubvectorWordIterator<Iterator, ConstIterator>;

	const BitSubvector<Iterator, ConstIterator> *_v;
	typename Iterator::const_word_iterator _pos;
	value_type _word;
	bool _ref_valid;

	inline void update_ref () {
		if (!_ref_valid) {
			if (_v->_begin.pos () == 0)
				_word = *_pos;
			else if (_pos + 1 == (_v->_end.pos () ? (_v->_end.word () + 1) : _v->_end.word ()))
				_word = Endianness::shift_left (*_pos, _v->_begin.pos ());
			else
				_word = Endianness::shift_left (*_pos, _v->_begin.pos ()) | Endianness::shift_right (_pos[1], WordTraits<value_type>::bits - _v->_begin.pos ());

			if (((_pos - _v->_begin.word () + 1) << WordTraits<value_type>::logof_size) > _v->_end - _v->_begin)
				_word &= Endianness::mask_left ((_v->_end - _v->_begin) & WordTraits<value_type>::pos_mask);

			_ref_valid = true;
		}
	}
};

} // namespace LinBox

#endif // __VECTOR_BIT_SUBVECTOR_INL
