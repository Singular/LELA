/* lela/vector/shifted-vector.h
 * Copyright 2010 Bradford Hovinen
 *
 * Wrapper around std::vector which shifts entries by a fixed amount
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_SHIFTED_VECTOR_H
#define __LELA_VECTOR_SHIFTED_VECTOR_H

#include <iterator>

#include "lela/vector/sparse.h"

namespace LELA
{

template <class Iterator>
class ShiftedVectorIterator
{
    public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef ShiftedProperty<Iterator> reference;
	typedef const ShiftedProperty<Iterator> const_reference;
	typedef ShiftedProperty<Iterator> *pointer;
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

	ShiftedVectorIterator () {}
	ShiftedVectorIterator (Iterator pos, value_type shift) : _ref (pos, shift) {}

	ShiftedVectorIterator (const ShiftedVectorIterator &i)
		: _ref (i._ref._i, i._ref._shift)
		{}

	ShiftedVectorIterator &operator = (const ShiftedVectorIterator &i)
	{
		_ref._i = i._ref._i;
		_ref._shift = i._ref._shift;
		return *this;
	}

	template <class It>
	ShiftedVectorIterator &operator = (const ShiftedVectorIterator<It> &i)
	{
		_ref._i = i._ref._i;
		_ref._shift = i._ref._shift;
		return *this;
	}

	ShiftedVectorIterator &operator ++ () 
		{ ++_ref._i; return *this; }

	ShiftedVectorIterator operator ++ (int) 
	{
		ShiftedVectorIterator tmp (*this);
		++*this;
		return tmp;
	}

	ShiftedVectorIterator operator + (difference_type i) const
		{ return ShiftedVectorIterator (_ref._i + i, _ref._shift); }

	ShiftedVectorIterator &operator += (difference_type i) 
		{ _ref._i += i; return *this; }

	ShiftedVectorIterator &operator -- () 
		{ --_ref._i; return *this; }

	ShiftedVectorIterator operator -- (int) 
	{
		ShiftedVectorIterator tmp (*this);
		--*this;
		return tmp;
	}

	ShiftedVectorIterator operator - (difference_type i) const
		{ return *this + -i; }

	ShiftedVectorIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	template <class It>
	difference_type operator - (const ShiftedVectorIterator<It> &i) const 
		{ return _ref._i - i._ref._i; }

	reference operator [] (long i)
		{ return *(*this + i); }

	reference operator * ()
		{ return _ref; }

	const reference operator * () const
		{ return _ref; }

	bool operator == (const ShiftedVectorIterator &c) const 
		{ return (_ref._i == c._ref._i) && (_ref._shift == c._ref._shift); }

	template <class It>
	bool operator == (const ShiftedVectorIterator<It> &c) const 
		{ return (_ref._i == c._ref._i) && (_ref._shift == c._ref._shift); }

	bool operator != (const ShiftedVectorIterator &c) const 
		{ return (_ref._i != c._ref._i) || (_ref._shift != c._ref._shift); }

	template <class It>
	bool operator != (const ShiftedVectorIterator<It> &c) const 
		{ return (_ref._i != c._ref._i) || (_ref._shift != c._ref._shift); }

	template <class It>
	operator ShiftedVectorIterator<It> () const
		{ return ShiftedVectorIterator<It> (_ref._i, _ref._shift); }

    private:
	template <class It>
	friend class ShiftedVectorIterator;

	template <class Vector>
	friend class ShiftedVector;

	reference _ref;
};

/** Shifted vector (const version)
 *
 * This class mimics an STL-vector but shifted all elements by a fixed amount
 */
template <class Iterator>
class ConstShiftedVector
{
public:
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	typedef size_t      size_type;
	typedef ptrdiff_t   difference_type;
	typedef typename std::iterator_traits<Iterator>::reference reference;
	typedef const reference const_reference;

	typedef ShiftedVectorIterator<Iterator> const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	friend class ShiftedVectorIterator<Iterator>;

	ConstShiftedVector () {}
	ConstShiftedVector (Iterator begin, Iterator end, value_type shift)
		: _begin (begin), _end (end), _shift (shift)
		{}

	template <class It>
	ConstShiftedVector (const ConstShiftedVector<It> &v)
		: _begin (v._begin), _end (v._end), _shift (v._shift)
		{}

	ConstShiftedVector &operator = (const ConstShiftedVector &v)
	{
		_begin = v._begin;
		_end = v._end;
		_shift = v._shift;
		return *this;
	}

	template <class It>
	ConstShiftedVector &operator = (const ConstShiftedVector<It> &v)
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

	inline value_type front     () const { return *(begin ()); }
	inline value_type back      () const { return *(end () - 1); }

	inline size_t size () const { return _end - _begin; }
	inline bool empty () const { return _end == _begin; }

protected:
	template <class It>
	friend class ConstShiftedVector;

	Iterator _begin, _end;
	value_type _shift;
};

/** Shifted vector (mutable version)
 *
 * This class mimics an STL-vector but shifted all elements by a fixed amount
 */
template <class Vector>
class ShiftedVector
{
public:
	typedef typename std::iterator_traits<typename Vector::iterator>::value_type      value_type;
	typedef typename std::iterator_traits<typename Vector::iterator>::reference       reference;
	typedef const typename std::iterator_traits<typename Vector::iterator>::reference const_reference;
	typedef size_t      size_type;
	typedef ptrdiff_t   difference_type;

	typedef ShiftedVectorIterator<typename Vector::iterator> iterator;
	typedef ShiftedVectorIterator<typename Vector::const_iterator> const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	ShiftedVector () {}
	ShiftedVector (Vector &v, value_type shift)
		: _v (&v), _shift (shift)
		{}

	ShiftedVector (const ShiftedVector &v)
		: _v (v._v), _shift (v._shift) {}

	ShiftedVector &operator = (const ShiftedVector &v)
	{
		_v = v._v;
		_shift = v._shift;
		return *this;
	}

	inline iterator               begin  ()       { return iterator (_v->begin (), _shift); }
	inline iterator               end    ()       { return iterator (_v->end (), _shift); }
	inline const_iterator         begin  () const { return const_iterator (_v->begin (), _shift); }
	inline const_iterator         end    () const { return const_iterator (_v->end (), _shift); }

	inline reverse_iterator       rbegin ()       { return reverse_iterator (end ()); }
	inline reverse_iterator       rend   ()       { return reverse_iterator (begin ()); }
	inline const_reverse_iterator rbegin () const { return const_reverse_iterator (end ()); }
	inline const_reverse_iterator rend   () const { return const_reverse_iterator (begin ()); }

	inline reference       operator[] (size_type n)       { return *(begin () + n); }
	inline const_reference operator[] (size_type n) const { return *(begin () + n); }

	inline reference at (size_type n)
	{
		if (n >= size ())
			throw std::out_of_range ("n");
		else
			return (*this)[n];
	}

	inline const_reference at (size_type n) const
	{
		if (n >= size ())
			throw std::out_of_range ("n");
		else
			return (*this)[n];
	}

	inline value_type front () const { return *(begin ()); }
	inline value_type back  () const { return *(end () - 1); }

	inline void            push_back (const value_type &x)  { _v->push_back (x + _shift); }
	inline void            clear     ()            { _v->clear (); }
	inline void            resize    (size_type s) { _v->resize (s); }

	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
	{
		clear ();

		while (first != last)
			insert (end (), *first++);
	}

	inline iterator insert (iterator pos, const value_type &x)
		{ return iterator (_v->insert (pos._ref._i, x + _shift), _shift); }

	template <class It>
	void insert (iterator pos, It begin, It end)
	{
		while (begin != end) {
			pos = insert (pos, *begin++);
			++pos;
		}
	}

	inline iterator erase (iterator pos)
		{ return _v->erase (pos._pos); }

	template <class It>
	void erase (It begin, It end)
		{ _v->erase (begin, end); }

	inline size_type       size      () const { return _v->size ();  }
	inline bool            empty     () const { return _v->empty (); }
	inline size_type       max_size  () const { return _v->max_size ();  }

	inline bool operator == (const ShiftedVector &v) const
		{ return (_v == v._v) && (_shift == v._shift); }

	Vector       &parent     ()       { return *_v; }
	const Vector &parent     () const { return *_v; }
	value_type    shift      () const { return _shift; }

	void          set_parent (Vector &v)        { _v = &v; }
	void          set_shift  (value_type shift) { _shift = shift; }

protected:
	Vector *_v;
	value_type _shift;
};

} // namespace LELA

#endif // __LELA_VECTOR_SHIFTED_VECTOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
