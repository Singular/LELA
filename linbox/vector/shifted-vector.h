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
template <class Vector>
class ShiftedVector
{
public:
	typedef typename Vector::value_type value_type;
	typedef size_t             size_type;
	typedef long               difference_type;

	typedef ShiftedVectorConstIterator<Vector> const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	friend class ShiftedVectorConstIterator<Vector>;

	ShiftedVector () {}
	ShiftedVector (const Vector &v, value_type shift, value_type end)
		: _shift (shift)
		{ set_start_end (v.begin (), v.end (), shift, end); }

	ShiftedVector (const ShiftedVector &v, value_type shift, value_type end)
		: _shift (v._shift + shift)
		{ set_start_end (v._start, v._end, v._shift + shift, v._shift + end); }

	ShiftedVector (const ShiftedVector &v)
		: _shift (v._shift), _start (v._start), _end (v._end) {}

	inline const_iterator              begin      (void) const { return const_iterator (_start, _shift); }
	inline const_iterator              end        (void) const { return const_iterator (_end, _shift); }

	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (end ()); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (begin ()); }

	inline size_t size () const { return _end - _start; }

	value_type _shift;
	typename Vector::const_iterator _start, _end;

protected:
	void set_start_end (typename Vector::const_iterator v_begin, typename Vector::const_iterator v_end, value_type shift, value_type end)
	{
		_start = std::lower_bound (v_begin, v_end, shift);
		_end = std::upper_bound (v_begin, v_end, end);
	}
};

} // namespace LinBox

namespace std 
{
	template <class Vector>
	struct iterator_traits<LinBox::ShiftedVectorConstIterator<Vector> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename Vector::reference reference;
		typedef typename Vector::pointer pointer;
		typedef typename Vector::value_type value_type;
		typedef long difference_type;
	};
}

namespace LinBox
{

template <class Vector>
class ShiftedVectorConstIterator
{
    public:
	typedef typename std::iterator_traits<ShiftedVectorConstIterator>::iterator_category iterator_category;
	typedef typename std::iterator_traits<ShiftedVectorConstIterator>::reference reference;
	typedef typename std::iterator_traits<ShiftedVectorConstIterator>::pointer pointer;
	typedef typename std::iterator_traits<ShiftedVectorConstIterator>::value_type value_type;
	typedef typename std::iterator_traits<ShiftedVectorConstIterator>::difference_type difference_type;

	ShiftedVectorConstIterator () {}
	ShiftedVectorConstIterator (const typename Vector::const_iterator &pos, typename Vector::value_type shift) : _pos (pos), _shift (shift) {}
	ShiftedVectorConstIterator (const ShiftedVectorConstIterator &i) : _pos (i._pos), _shift (i._shift) {}

	ShiftedVectorConstIterator &operator = (const ShiftedVectorConstIterator &i) {
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

	difference_type operator - (ShiftedVectorConstIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i) - _shift; }

	typename Vector::value_type operator * () const
		{ return *_pos - _shift; }

	bool operator == (const ShiftedVectorConstIterator &c) const 
		{ return (_pos == c._pos) && (_shift == c._shift); }

	bool operator != (const ShiftedVectorConstIterator &c) const 
		{ return (_pos != c._pos) || (_shift != c._shift); }

    private:
	typename Vector::const_iterator _pos;
	typename Vector::value_type _shift;
};

} // namespace LinBox

#endif // __VECTOR_SHIFTED_VECTOR_H
