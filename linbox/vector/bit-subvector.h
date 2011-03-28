/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/bit-subvector.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from bit-vector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_BIT_SUBVECTOR_H
#define __VECTOR_BIT_SUBVECTOR_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"

namespace LinBox
{

template <class Iterator, class ConstIterator = Iterator>
class BitSubvectorWordIterator;

template <class Iterator, class ConstIterator = Iterator>
class BitSubvectorConstWordIterator;

/** A subvector of a BitVector. The interface is identical to that of \ref BitVector.
 \ingroup vector
 */

template <class Iterator, class ConstIterator = Iterator>
class BitSubvector
{
    public:
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
	typedef typename std::iterator_traits<Iterator>::pointer	 pointer;
	typedef typename std::iterator_traits<Iterator>::reference	 reference;
	typedef Iterator iterator;
	typedef ConstIterator const_iterator;

	typedef typename std::iterator_traits<iterator>::size_type	 size_type;
	typedef typename std::iterator_traits<const_iterator>::const_reference const_reference;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	class word_reference;
	typedef BitSubvectorWordIterator<Iterator, ConstIterator> word_iterator;
	typedef BitSubvectorConstWordIterator<Iterator, ConstIterator> const_word_iterator;

	typedef std::reverse_iterator<word_iterator> reverse_word_iterator;
	typedef std::reverse_iterator<const_word_iterator> const_reverse_word_iterator;

	friend class word_reference;
	friend class BitSubvectorConstWordIterator<Iterator, ConstIterator>;
	friend class BitSubvectorWordIterator<Iterator, ConstIterator>;

	BitSubvector () {}

	template <class Endianness>
	BitSubvector (BitVector<Endianness> &v, size_t start, size_t end)
		: _begin (v.begin () + start), _end (v.begin () + end) { set_end_word (); }

	BitSubvector (BitSubvector &v, size_t start, size_t end)
		: _begin (v._begin + start), _end (v._begin + end) { set_end_word (); }

	BitSubvector (iterator begin, iterator end)
		: _begin (begin), _end (end) { set_end_word (); }
		
	BitSubvector (const BitSubvector &v) 
		: _begin (v._begin), _end (v._end), _end_word (v._end_word), _end_marker (v._end_marker) {}

	~BitSubvector () {}

	inline iterator                    begin      (void)       { return _begin; }
	inline const_iterator              begin      (void) const { return const_iterator (_begin); }
	inline iterator                    end        (void)       { return _end; }
	inline const_iterator              end        (void) const { return const_iterator (_end); }

	inline reverse_iterator            rbegin     (void)       { return reverse_iterator (_end); }
	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (_end); }
	inline reverse_iterator            rend       (void)       { return reverse_iterator (_begin); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (_begin); }

	inline word_iterator               wordBegin  (void)       { return word_iterator (*this, _begin.word ()); }
	inline const_word_iterator         wordBegin  (void) const { return const_word_iterator (*this, _begin.word ()); }
	inline word_iterator               wordEnd    (void)       { return word_iterator (*this, _end_word); }
	inline const_word_iterator         wordEnd    (void) const { return const_word_iterator (*this, _end_word); }

	inline reverse_word_iterator       wordRbegin (void)       { return reverse_word_iterator (wordEnd ()); }
	inline const_reverse_word_iterator wordRbegin (void) const { return const_reverse_word_iterator (wordEnd ()); }
	inline reverse_word_iterator       wordRend   (void)       { return reverse_word_iterator (wordBegin ()); }
	inline const_reverse_word_iterator wordRend   (void) const { return const_reverse_word_iterator (wordBegin ()); }

	// Element access

	inline reference       operator[] (size_type n)       { return *(_begin + n); }
	inline const_reference operator[] (size_type n) const { return *(_begin + n); }

	inline reference at (size_type n)
	{   
		iterator p = _begin + n;
		if ( _begin <= p && p < _end ) 
			return *p;
		else 
			throw std::out_of_range(); //out of range error message.
	}

	inline const_reference at (size_type n) const
	{
		const_iterator p = _begin + n;
		if ( _begin <= p && p < _end)
			return *p;
		else 
			throw std::out_of_range(); //out of range error message
	}

	BitSubvector &operator = (const BitSubvector &sub)
		{ _begin = sub._begin; _end = sub._end; _end_word = sub._end_word; _end_marker = sub._end_marker; return *this; }

	inline reference       front (void)       { return *_begin; }
	inline const_reference front (void) const { return *_begin; }
	inline reference       back  (void)       { return *(_end - 1); }
	inline const_reference back  (void) const { return *(_end - 1); }

	inline size_type size      (void) const { return _end - _begin;    }
	inline bool      empty     (void) const { return _end == _begin;   }
	inline size_type max_size  (void) const { return _end - _begin; }

	inline bool operator == (const BitSubvector &v) const
		{ return (_begin == v._begin) && (_end == v._end); }
	inline bool operator != (const BitSubvector &v) const
		{ return (_begin != v._begin) || (_end != v._end); }

    protected:

	Iterator     _begin;
	Iterator     _end;

	typename Iterator::word_iterator _end_word, _end_marker;

	void set_end_word ()
	{
		_end_marker = _end_word = _end.word ();

		if (_end.pos () > _begin.pos ())
			++_end_word;
		if (_end.pos () > 0)
			++_end_marker;
	}

}; // template <class Iterator> class BitSubvector

// Vector traits for BitSubvector
template <class Iterator, class ConstIterator>
struct GF2VectorTraits<BitSubvector<Iterator, ConstIterator> >
{ 
	typedef BitSubvector<Iterator, ConstIterator> VectorType;
	typedef VectorCategories::DenseZeroOneVectorTag VectorCategory; 
};

// Vector traits for BitSubvector
template <class Iterator, class ConstIterator>
struct GF2VectorTraits<const BitSubvector<Iterator, ConstIterator> >
{ 
	typedef BitSubvector<Iterator, ConstIterator> VectorType;
	typedef VectorCategories::DenseZeroOneVectorTag VectorCategory; 
};

} // namespace LinBox

#include "linbox/vector/bit-subvector.inl"

#endif // __VECTOR_BIT_SUBVECTOR_H
