/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/shifted-vector.h
 * Copyright 2010 Bradford Hovinen
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SHIFTED_VECTOR_H
#define __VECTOR_SHIFTED_VECTOR_H

#include <iterator>

namespace LinBox
{

template <class Vector>
class ShiftedVectorConstIterator;

// This class mimics a vector but shifts all of its elements by a fixed value. Used for sparse subvectors.
template <class Iterator>
class ShiftedVector
{
public:
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	typedef size_t size_type;
	typedef long   difference_type;

	typedef ShiftedVectorConstIterator<Iterator> const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	friend class ShiftedVectorConstIterator<Iterator>;

	ShiftedVector () {}
	ShiftedVector (Iterator begin, Iterator end, value_type shift)
		: _begin (begin), _end (end), _shift (shift)
		{}

	template <class It>
	ShiftedVector (const ShiftedVector<It> &v)
		: _begin (v._begin), _end (v._end), _shift (v._shift)
		{}

	ShiftedVector &operator = (const ShiftedVector &v)
	{
		_begin = v._begin;
		_end = v._end;
		_shift = v._shift;
		return *this;
	}

	template <class It>
	ShiftedVector &operator = (const ShiftedVector<It> &v)
	{
		_begin = v._begin;
		_end = v._end;
		_shift = v._shift;
		return *this;
	}

	inline const_iterator              begin      (void) const { return const_iterator (_begin, _shift); }
	inline const_iterator              end        (void) const { return const_iterator (_end, _shift); }

	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (end ()); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (begin ()); }

	inline size_t size () const { return _end - _begin; }
	inline bool empty () const { return _end == _begin; }

protected:
	template <class It>
	friend class ShiftedVector;

	Iterator _begin, _end;
	value_type _shift;
};

template <class Iterator>
class ShiftedVectorConstIterator
{
    public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename std::iterator_traits<Iterator>::value_type reference;
	typedef typename std::iterator_traits<Iterator>::pointer pointer;
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	typedef ptrdiff_t difference_type;

	ShiftedVectorConstIterator () {}
	ShiftedVectorConstIterator (Iterator pos, value_type shift) : _pos (pos), _shift (shift) {}

	template <class It>
	ShiftedVectorConstIterator (ShiftedVectorConstIterator<It> i)
		: _pos (i._pos), _shift (i._shift)
		{}

	ShiftedVectorConstIterator &operator = (const ShiftedVectorConstIterator &i) {
		_pos = i._pos;
		_shift = i._shift;
		return *this;
	}

	template <class It>
	ShiftedVectorConstIterator &operator = (const ShiftedVectorConstIterator<It> &i) {
		_pos = i._pos;
		_shift = i._shift;
		return *this;
	}

	ShiftedVectorConstIterator &operator ++ () 
	{
		++_pos;
		return *this;
	}

	ShiftedVectorConstIterator operator ++ (int) 
	{
		ShiftedVectorConstIterator tmp (*this);
		++*this;
		return tmp;
	}

	ShiftedVectorConstIterator operator + (difference_type i) const
	{
		return ShiftedVectorConstIterator (_pos + i, _shift);
	}

	ShiftedVectorConstIterator &operator += (difference_type i) 
	{
		_pos += i;
		return *this;
	}

	ShiftedVectorConstIterator &operator -- () 
	{
		--_pos;
		return *this;
	}

	ShiftedVectorConstIterator operator -- (int) 
	{
		ShiftedVectorConstIterator tmp (*this);
		--*this;
		return tmp;
	}

	ShiftedVectorConstIterator operator - (difference_type i) const
		{ return *this + -i; }

	ShiftedVectorConstIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	template <class It>
	difference_type operator - (const ShiftedVectorConstIterator<It> &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i) - _shift; }

	value_type operator * () const
		{ return *_pos - _shift; }

	bool operator == (const ShiftedVectorConstIterator &c) const 
		{ return (_pos == c._pos) && (_shift == c._shift); }

	template <class It>
	bool operator == (const ShiftedVectorConstIterator<It> &c) const 
		{ return (_pos == c._pos) && (_shift == c._shift); }

	bool operator != (const ShiftedVectorConstIterator &c) const 
		{ return (_pos != c._pos) || (_shift != c._shift); }

	template <class It>
	bool operator != (const ShiftedVectorConstIterator<It> &c) const 
		{ return (_pos != c._pos) || (_shift != c._shift); }

    private:
	template <class It>
	friend class ShiftedVectorConstIterator;

	Iterator _pos;
	value_type _shift;
};

} // namespace LinBox

#endif // __VECTOR_SHIFTED_VECTOR_H
