/* lela/vector/bit-subvector-word-aligned.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from bit-subvector.h
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_BIT_SUBVECTOR_WORD_ALIGNED_H
#define __LELA_VECTOR_BIT_SUBVECTOR_WORD_ALIGNED_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "lela/vector/traits.h"
#include "lela/vector/bit-iterator.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/bit-subvector.h"

namespace LELA
{

/** A subvector of a \ref BitVector which is restricted to
 * word-aligned subvectors. More efficient than \ref BitSubvector, which
 * does not have this restriction.
 \ingroup vector
 */

template <class Iterator, class ConstIterator, class _Endianness>
class BitSubvectorWordAligned
{
    public:
	typedef Iterator word_iterator;
	typedef ConstIterator const_word_iterator;
	typedef _Endianness Endianness;

	typedef BitVectorIterator<word_iterator, const_word_iterator, Endianness> iterator;
	typedef BitVectorConstIterator<const_word_iterator, Endianness>           const_iterator;
	typedef BitVectorReference<word_iterator, Endianness>                     reference;
	typedef BitVectorConstReference<const_word_iterator, Endianness>          const_reference;
	typedef typename std::iterator_traits<Iterator>::pointer pointer;
	typedef typename BitVectorIterator<word_iterator, const_word_iterator>::difference_type difference_type;
	typedef typename BitVectorIterator<word_iterator, const_word_iterator>::size_type	size_type;

	typedef typename std::iterator_traits<word_iterator>::value_type word_type;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef std::reverse_iterator<word_iterator> reverse_word_iterator;
	typedef std::reverse_iterator<const_word_iterator> const_reverse_word_iterator;

	typedef VectorRepresentationTypes::Dense01 RepresentationType; 
	typedef VectorStorageTypes::Real StorageType;
	typedef BitVector<Endianness> ContainerType;
	typedef BitSubvector<iterator, const_iterator> SubvectorType;
	typedef BitSubvector<const_iterator, const_iterator> ConstSubvectorType;
	typedef BitSubvectorWordAligned<word_iterator, const_word_iterator> AlignedSubvectorType;
	typedef BitSubvectorWordAligned<const_word_iterator, const_word_iterator> ConstAlignedSubvectorType;
	static const int align = WordTraits<word_type>::bits;

	BitSubvectorWordAligned () {}

	BitSubvectorWordAligned (BitVector<Endianness> &v, size_t start, size_t finish)
		: _begin (v.word_begin () + (start >> WordTraits<word_type>::logof_size)),
		  _end (v.word_begin () + ((finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size)),
		  _bit_len (finish - start)
	{
		lela_check ((start & WordTraits<word_type>::pos_mask) == 0);
		lela_check (finish == v.size () || (finish & WordTraits<word_type>::pos_mask) == 0);
	}

	BitSubvectorWordAligned (const BitVector<Endianness> &v, size_t start, size_t finish)
		: _begin (v.word_begin () + (start >> WordTraits<word_type>::logof_size)),
		  _end (v.word_begin () + ((finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size)),
		  _bit_len (finish - start)
	{
		lela_check ((start & WordTraits<word_type>::pos_mask) == 0);
		lela_check (finish == v.size () || (finish & WordTraits<word_type>::pos_mask) == 0);
	}

	template <class It, class CIt>
	BitSubvectorWordAligned (const BitSubvectorWordAligned<It, CIt, Endianness> &v, size_t start, size_t finish)
		: _begin (v._begin + (start >> WordTraits<word_type>::logof_size)),
		  _end (v._begin + ((finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size)),
		  _bit_len (finish - start)
	{
		lela_check ((start & WordTraits<word_type>::pos_mask) == 0);
		lela_check (finish == v.size () || (finish & WordTraits<word_type>::pos_mask) == 0);
	}

	BitSubvectorWordAligned (Iterator begin, Iterator end, size_t bit_len = 0)
		: _begin (begin), _end (end)
	{
		if (bit_len == 0)
			_bit_len = (end - begin) * WordTraits<word_type>::bits;
		else
			_bit_len = bit_len;
	}
		
	BitSubvectorWordAligned (const BitSubvectorWordAligned &v) 
		: _begin (v._begin), _end (v._end), _bit_len (v._bit_len) {}

	template <class It, class CIt>
	BitSubvectorWordAligned (const BitSubvectorWordAligned<It, CIt, Endianness> &v) 
		: _begin (v._begin), _end (v._end), _bit_len (v._bit_len) {}

	~BitSubvectorWordAligned () {}

	inline iterator                    begin      (void)       { return iterator (_begin, 0); }
	inline const_iterator              begin      (void) const { return const_iterator (_begin, 0); }
	inline iterator                    end        (void)
	{
		if (!(_bit_len & WordTraits<word_type>::pos_mask))
			return iterator (_end, 0);
		else
			return iterator (_end - 1, _bit_len & WordTraits<word_type>::pos_mask);
	}
	inline const_iterator              end        (void) const
	{
		if (!(_bit_len & WordTraits<word_type>::pos_mask))
			return const_iterator (_end, 0);
		else
			return const_iterator (_end - 1, _bit_len & WordTraits<word_type>::pos_mask);
	}

	inline reverse_iterator            rbegin     (void)       { return reverse_iterator (_end); }
	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (_end); }
	inline reverse_iterator            rend       (void)       { return reverse_iterator (_begin); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (_begin); }

	inline word_iterator               word_begin  (void)       { return _begin; }
	inline const_word_iterator         word_begin  (void) const { return _begin; }
	inline word_iterator               word_end    (void)       { return (_begin == _end) ? _end : _end - 1; }
	inline const_word_iterator         word_end    (void) const { return (_begin == _end) ? _end : _end - 1; }

	inline reverse_word_iterator       word_rbegin (void)       { return reverse_word_iterator (word_end ()); }
	inline const_reverse_word_iterator word_rbegin (void) const { return const_reverse_word_iterator (word_end ()); }
	inline reverse_word_iterator       word_rend   (void)       { return reverse_word_iterator (word_begin ()); }
	inline const_reverse_word_iterator word_rend   (void) const { return const_reverse_word_iterator (word_begin ()); }

	// Element access

	inline reference       operator[] (size_type n)       { return *(iterator (_begin + n / WordTraits<word_type>::bits, n & WordTraits<word_type>::pos_mask)); }
	inline const_reference operator[] (size_type n) const { return *(const_iterator (_begin + n / WordTraits<word_type>::bits, n & WordTraits<word_type>::pos_mask)); }

	inline reference at (size_type n)
	{
		if (n < _bit_len)
			return *iterator (_begin + n / WordTraits<word_type>::bits, n & WordTraits<word_type>::pos_mask);
		else 
			throw std::out_of_range (std::string ("n"));
	}

	inline const_reference at (size_type n) const
	{
		if (n < _bit_len)
			return *const_iterator (_begin + n / WordTraits<word_type>::bits, n & WordTraits<word_type>::pos_mask);
		else 
			throw std::out_of_range (std::string ("n"));
	}

	BitSubvectorWordAligned &operator = (const BitSubvectorWordAligned &sub)
		{ _begin = sub._begin; _end = sub._end; _bit_len = sub._bit_len; return *this; }

	inline reference       front (void)       { return *iterator (_begin, 0); }
	inline const_reference front (void) const { return *const_iterator (_begin, 0); }
	inline reference       back  (void)       { return *iterator (_end, _bit_len & WordTraits<word_type>::pos_mask); }
	inline const_reference back  (void) const { return *const_iterator (_end, _bit_len & WordTraits<word_type>::pos_mask); }

	inline word_type       &front_word (void)       { return *_begin; }
	inline const word_type &front_word (void) const { return *_begin; }
	inline word_type       &back_word  (void)       { return *(_end - 1); }
	inline const word_type &back_word  (void) const { return *(_end - 1); }

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

	template <class It, class CIt, class E>
	friend class BitSubvectorWordAligned;

}; // template <class Iterator> class BitSubvectorWordAligned

} // namespace LELA

namespace std {

template<class Iterator, class ConstIterator, class Endianness>
void swap (LELA::BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> &x, LELA::BitSubvectorWordAligned<Iterator, ConstIterator, Endianness> &y)
{
	typename LELA::BitSubvectorWordAligned<Iterator, ConstIterator, Endianness>::word_iterator i_x, i_y;

	for (i_x = x.word_begin (), i_y = y.word_begin (); i_x != x.word_end () && i_y != y.word_end (); ++i_x, ++i_y)
		std::swap (*i_x, *i_y);

	std::swap (x.back_word (), y.back_word ());
}

} // namespace std

#endif // __LELA_VECTOR_BIT_SUBVECTOR_WORD_ALIGNED_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
