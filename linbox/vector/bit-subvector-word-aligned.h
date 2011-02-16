/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/bit-subvector-word-aligned.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from bit-subvector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_BIT_SUBVECTOR_WORD_ALIGNED_H
#define __VECTOR_BIT_SUBVECTOR_WORD_ALIGNED_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-iterator.h"

namespace LinBox
{

/** A subvector of a \ref BitVector which is restricted to
 * word-aligned subvectors. More efficient than \ref BitSubvector, which
 * does not have this restriction.
 \ingroup vector
 */

template <class Iterator, class ConstIterator, class Endianness>
class BitSubvectorWordAligned
{
    public:
	typedef Iterator word_iterator;
	typedef ConstIterator const_word_iterator;

	typedef BitVectorIterator<word_iterator, const_word_iterator, Endianness> iterator;
	typedef BitVectorConstIterator<const_word_iterator, Endianness>           const_iterator;
	typedef BitVectorReference<word_iterator, Endianness>                     reference;
	typedef BitVectorConstReference<const_word_iterator, Endianness>          const_reference;
	typedef iterator                                              pointer;
	typedef typename BitVectorIterator<word_iterator, const_word_iterator>::difference_type difference_type;
	typedef typename BitVectorIterator<word_iterator, const_word_iterator>::size_type	   size_type;

	typedef typename std::iterator_traits<word_iterator>::value_type word;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef std::reverse_iterator<word_iterator> reverse_word_iterator;
	typedef std::reverse_iterator<const_word_iterator> const_reverse_word_iterator;

	BitSubvectorWordAligned () {}

	BitSubvectorWordAligned (Iterator begin, Iterator end, size_t bit_len = 0)
		: _begin (begin), _end (end)
	{
		if (bit_len == 0)
			_bit_len = (end - begin) * WordTraits<word>::bits;
		else
			_bit_len = bit_len;
	}
		
	BitSubvectorWordAligned (const BitSubvectorWordAligned &v) 
		: _begin (v._begin), _end (v._end), _bit_len (v._bit_len) {}

	~BitSubvectorWordAligned () {}

	inline iterator                    begin      (void)       { return iterator (_begin, 0); }
	inline const_iterator              begin      (void) const { return const_iterator (_begin, 0); }
	inline iterator                    end        (void)
	{
		if (!(_bit_len & WordTraits<word>::pos_mask))
			return iterator (_end, 0);
		else
			return iterator (_end - 1, _bit_len & WordTraits<word>::pos_mask);
	}
	inline const_iterator              end        (void) const
	{
		if (!(_bit_len & WordTraits<word>::pos_mask))
			return const_iterator (_end, 0);
		else
			return const_iterator (_end - 1, _bit_len & WordTraits<word>::pos_mask);
	}

	inline reverse_iterator            rbegin     (void)       { return reverse_iterator (_end); }
	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (_end); }
	inline reverse_iterator            rend       (void)       { return reverse_iterator (_begin); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (_begin); }

	inline word_iterator               wordBegin  (void)       { return _begin; }
	inline const_word_iterator         wordBegin  (void) const { return _begin; }
	inline word_iterator               wordEnd    (void)       { return _end; }
	inline const_word_iterator         wordEnd    (void) const { return _end; }

	inline reverse_word_iterator       wordRbegin (void)       { return reverse_word_iterator (wordEnd ()); }
	inline const_reverse_word_iterator wordRbegin (void) const { return const_reverse_word_iterator (wordEnd ()); }
	inline reverse_word_iterator       wordRend   (void)       { return reverse_word_iterator (wordBegin ()); }
	inline const_reverse_word_iterator wordRend   (void) const { return const_reverse_word_iterator (wordBegin ()); }

	// Element access

	inline reference       operator[] (size_type n)       { return *(iterator (_begin + n / WordTraits<word>::bits, n & WordTraits<word>::pos_mask)); }
	inline const_reference operator[] (size_type n) const { return *(const_iterator (_begin + n / WordTraits<word>::bits, n & WordTraits<word>::pos_mask)); }

	inline reference at (size_type n)
	{
		if (n < _bit_len)
			return *iterator (_begin + n / WordTraits<word>::bits, n & WordTraits<word>::pos_mask);
		else 
			throw std::out_of_range(); //out of range error message.
	}

	inline const_reference at (size_type n) const
	{
		if (n < _bit_len)
			return *const_iterator (_begin + n / WordTraits<word>::bits, n & WordTraits<word>::pos_mask);
		else 
			throw std::out_of_range(); //out of range error message
	}

	BitSubvectorWordAligned &operator = (const BitSubvectorWordAligned &sub)
		{ _begin = sub._begin; _end = sub._end; _bit_len = sub._bit_len; return *this; }

	inline reference       front (void)       { return *iterator (_begin, 0); }
	inline const_reference front (void) const { return *const_iterator (_begin, 0); }
	inline reference       back  (void)       { return *iterator (_end, _bit_len & WordTraits<word>::pos_mask); }
	inline const_reference back  (void) const { return *const_iterator (_end, _bit_len & WordTraits<word>::pos_mask); }

	inline word       &front_word (void)       { return *_begin; }
	inline const word &front_word (void) const { return *_begin; }
	inline word       &back_word  (void)       { return *_end; }
	inline const word &back_word  (void) const { return *_end; }

	inline size_type size      (void) const { return _bit_len; }
	inline bool      empty     (void) const { return _bit_len == 0; }
	inline size_type max_size  (void) const { return _bit_len; }

	inline size_type word_size (void) const { return _end - _begin; }

	inline bool operator == (const BitSubvectorWordAligned &v) const
		{ return (_begin == v._begin) && (_end == v._end) && (_bit_len == v._bit_len); }
	inline bool operator != (const BitSubvectorWordAligned &v) const
		{ return (_begin != v._begin) || (_end != v._end) || (_bit_len != v._bit_len); }

    protected:

	Iterator     _begin;
	Iterator     _end;
	size_t       _bit_len;

}; // template <class Iterator> class BitSubvectorWordAligned

// Vector traits for BitVector wrapper
template <class Iterator, class ConstIterator, class Endianness>
struct VectorTraits<BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> >
{ 
	typedef BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> VectorType;
	typedef VectorCategories::DenseZeroOneVectorTag VectorCategory; 
};

template <class Iterator, class ConstIterator, class Endianness>
struct VectorTraits< std::pair<std::vector<size_t>, BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> > >
{ 
	typedef std::pair<std::vector<size_t>, BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> > VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

} // namespace LinBox

#endif // __VECTOR_BIT_SUBVECTOR_WORD_ALIGNED_H
