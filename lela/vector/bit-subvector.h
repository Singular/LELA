/* lela/vector/bit-subvector.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from bit-vector.h
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_BIT_SUBVECTOR_H
#define __LELA_VECTOR_BIT_SUBVECTOR_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "lela/util/debug.h"
#include "lela/vector/traits.h"
#include "lela/vector/bit-vector.h"

namespace LELA
{

/** A subvector of a BitVector. The interface is identical to that of \ref BitVector.
 \ingroup vector
 */

template <class Iterator, class ConstIterator>
class BitSubvector
{
    public:
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
	typedef typename std::iterator_traits<Iterator>::pointer	 pointer;
	typedef typename std::iterator_traits<Iterator>::reference	 reference;
	typedef Iterator iterator;
	typedef ConstIterator const_iterator;

	typedef VectorRepresentationTypes::Dense01 RepresentationType; 
	typedef VectorStorageTypes::Real StorageType;
	typedef BitVector<typename Iterator::Endianness> ContainerType;
	typedef BitSubvector<iterator, const_iterator> SubvectorType;
	typedef BitSubvector<const_iterator, const_iterator> ConstSubvectorType;
	typedef BitSubvector<iterator, const_iterator> AlignedSubvectorType;
	typedef BitSubvector<const_iterator, const_iterator> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef typename Iterator::size_type	 size_type;
	typedef typename Iterator::const_reference const_reference;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	class word_reference;
	class word_end_reference;
	class word_iterator;
	class const_word_iterator;

	typedef std::reverse_iterator<word_iterator> reverse_word_iterator;
	typedef std::reverse_iterator<const_word_iterator> const_reverse_word_iterator;

	friend class word_reference;
	friend class word_iterator;
	friend class const_word_iterator;

	typedef typename Iterator::Endianness Endianness;
	typedef typename std::iterator_traits<word_iterator>::value_type word_type;

	BitSubvector () {}

	template <class Endianness>
	BitSubvector (BitVector<Endianness> &v, size_t start, size_t end)
		: _begin (v.begin () + start), _end (v.begin () + end) { set_end_word (); }

	template <class It, class CIt>
	BitSubvector (BitSubvectorWordAligned<It, CIt, Endianness> &v, size_t start, size_t end)
		: _begin (v.begin () + start), _end (v.begin () + end) { set_end_word (); }

	BitSubvector (BitSubvector &v, size_t start, size_t end)
		: _begin (v._begin + start), _end (v._begin + end) { set_end_word (); }

	BitSubvector (iterator begin, iterator end)
		: _begin (begin), _end (end) { set_end_word (); }
		
	BitSubvector (const BitSubvector &v) 
		: _begin (v._begin), _end (v._end), _end_word (v._end_word), _end_mask (v._end_mask) {}

	~BitSubvector () {}

	inline iterator                    begin      (void)       { return _begin; }
	inline const_iterator              begin      (void) const { return const_iterator (_begin); }
	inline iterator                    end        (void)       { return _end; }
	inline const_iterator              end        (void) const { return const_iterator (_end); }

	inline reverse_iterator            rbegin     (void)       { return reverse_iterator (_end); }
	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (_end); }
	inline reverse_iterator            rend       (void)       { return reverse_iterator (_begin); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (_begin); }

	inline word_iterator               word_begin  (void)       { return word_iterator (_begin.word (), _begin.pos ()); }
	inline const_word_iterator         word_begin  (void) const { return const_word_iterator (_begin.word (), _begin.pos ()); }
	inline word_iterator               word_end    (void)       { return word_iterator (_end_word, _begin.pos ()); }
	inline const_word_iterator         word_end    (void) const { return const_word_iterator (_end_word, _begin.pos ()); }

	inline reverse_word_iterator       word_rbegin (void)       { return reverse_word_iterator (word_end ()); }
	inline const_reverse_word_iterator word_rbegin (void) const { return const_reverse_word_iterator (word_end ()); }
	inline reverse_word_iterator       word_rend   (void)       { return reverse_word_iterator (word_begin ()); }
	inline const_reverse_word_iterator word_rend   (void) const { return const_reverse_word_iterator (word_begin ()); }

	// Element access

	inline reference       operator[] (size_type n)       { return *(_begin + n); }
	inline const_reference operator[] (size_type n) const { return *(_begin + n); }

	inline reference at (size_type n)
	{   
		iterator p = _begin + n;
		if ( _begin <= p && p < _end ) 
			return *p;
		else 
			throw std::out_of_range (std::string ("n"));
	}

	inline const_reference at (size_type n) const
	{
		const_iterator p = _begin + n;
		if ( _begin <= p && p < _end)
			return *p;
		else 
			throw std::out_of_range (std::string ("n"));
	}

	BitSubvector &operator = (const BitSubvector &sub)
		{ _begin = sub._begin; _end = sub._end; _end_word = sub._end_word; _end_mask = sub._end_mask; return *this; }

	inline reference       front (void)       { return *_begin; }
	inline const_reference front (void) const { return *_begin; }
	inline reference       back  (void)       { return *(_end - 1); }
	inline const_reference back  (void) const { return *(_end - 1); }

	inline word_reference     front_word (void)       { return *word_begin (); }
	inline word_type          front_word (void) const { return *word_begin (); }
	inline word_end_reference back_word  (void)       { return word_end_reference (_end_word, _begin.pos (), _end_mask); }
	inline word_type          back_word  (void) const { return (word_type) word_end_reference (_end_word, _begin.pos (), _end_mask); }

	inline size_type size      (void) const { return _end - _begin;    }
	inline bool      empty     (void) const { return _end == _begin;   }
	inline size_type max_size  (void) const { return _end - _begin; }

	inline bool operator == (const BitSubvector &v) const
		{ return (_begin == v._begin) && (_end == v._end); }
	inline bool operator != (const BitSubvector &v) const
		{ return (_begin != v._begin) || (_end != v._end); }

	inline size_type word_size () const { return word_end () - word_begin () + 1; }

    protected:

	// BitVector::iterators pointing to the beginning and end of the bit-subvector, respectively
	Iterator     _begin;
	Iterator     _end;

	// Points to word which is past the end for the word-iterator and from which back_word is constructed
	typename Iterator::word_iterator _end_word;

	// Mask to apply to back_word based on position of end-iterator
	word_type _end_mask;

	// Note: We adopt the convention that the word_iterator always
	// ends *before* the last word and back_word never has a mask
	// of zero (unless the whole subvector is empty), even when
	// the size of the vector is a multiple of
	// WordTraits<word_type>::bits

	void set_end_word ()
	{
		if (_begin != _end && _begin.pos () >= _end.pos ())
			_end_word = _end.word () - 1;
		else
			_end_word = _end.word ();

		_end_mask = Endianness::mask_left ((_end.pos () - _begin.pos ()) & WordTraits<word_type>::pos_mask);

		if (_begin != _end && _end_mask == 0)
			_end_mask = (word_type) -1;
	}

    public:

	class word_reference
	{
	public:
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type word;
		typedef typename Iterator::Endianness Endianness;

		word_reference () {}

		word_reference (typename Iterator::word_iterator position, uint8 shift)
			: _pos (position), _shift (shift)
		{
			_m.parts.low = (word_type) -1;
			_m.parts.high = 0ULL;
			_m.full = Endianness::shift_right (_m.full, _shift);
		}

		word_reference (const word_reference &ref)
			: _pos (ref._pos), _shift (ref._shift)
			{ _m.full = ref._m.full; }

		~word_reference () {}

		word_reference &operator = (word_reference &a)
		{
			*this = (word) a;
			return a;
		}

		word_reference &operator = (word v)
		{
			typename Endianness::word_pair w;

			w.parts.low = v;
			w.parts.high = 0ULL;

			w.full = Endianness::shift_right (w.full, _shift);

			*_pos = (*_pos & ~_m.parts.low) | w.parts.low;
			_pos[1] = (_pos[1] & ~_m.parts.high) | w.parts.high;

			return *this;
		}

		word_reference &operator &= (word_reference &a)
			{ return *this &= (word) a; }

		word_reference &operator &= (word v) 
			{ return *this = (word) *this & v; }

		word_reference &operator |= (word_reference &a) 
			{ return *this |= (word) a; }

		word_reference &operator |= (word v) 
		{
			typename Endianness::word_pair w;

			w.parts.low = v;
			w.parts.high = 0ULL;

			w.full = Endianness::shift_right (w.full, _shift);

			*_pos |= w.parts.low;
			_pos[1] |= w.parts.high;

			return *this;
		}

		word_reference &operator ^= (word_reference &a) 
			{ return *this ^= ((word) a); }

		word_reference &operator ^= (word v) 
		{
			typename Endianness::word_pair w;

			w.parts.low = v;
			w.parts.high = 0ULL;

			w.full = Endianness::shift_right (w.full, _shift);

			*_pos ^= w.parts.low;
			_pos[1] ^= w.parts.high;

			return *this;
		}

		operator word (void) const
		{
			typename Endianness::word_pair v;

			v.parts.low = *_pos;
			v.parts.high = _pos[1];
			v.full = Endianness::shift_left (v.full, _shift);

			return v.parts.low;
		}

	private:
		friend class word_iterator;
		friend class const_word_iterator;
		friend class const_word_reference;

		/* Semantics: _pos and _pos + 1 must always point to a
		 * valid position in memory. If _pos + 1 is past the
		 * end, the corresponding word_iterator is too.
		 */
		typename Iterator::word_iterator _pos;

		/* Masks to be used when assigning to _pos */
		typename Endianness::word_pair _m;

		/* Offset of word-reference in bits */
		uint8 _shift;
	};

	class word_end_reference
	{
	public:
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type word;
		typedef typename Iterator::Endianness Endianness;

		word_end_reference () {}

		word_end_reference (typename Iterator::word_iterator position, uint8 shift, word mask)
			: _pos (position), _shift (shift), _mask (mask)
			{}

		word_end_reference (const word_end_reference &w)
			: _pos (w._pos), _shift (w._shift), _mask (w._mask)
			{}

		~word_end_reference () {}

		word_end_reference &operator = (const word_end_reference &a)
		{
			*this = (word) a;
			return *this;
		}

		word_end_reference &operator = (word v)
		{
			typename Endianness::word_pair w, m;

			w.parts.low = v & _mask;
			w.parts.high = 0ULL;
			m.parts.low = _mask;
			m.parts.high = 0ULL;

			w.full = Endianness::shift_right (w.full, _shift);
			m.full = Endianness::shift_right (m.full, _shift);

			*_pos = (*_pos & ~m.parts.low) | w.parts.low;

			if (m.parts.high != 0)
				_pos[1] = (_pos[1] & ~m.parts.high) | w.parts.high;

			return *this;
		}

		word_end_reference &operator &= (word_end_reference &a)
			{ return *this = (word) *this & (word) a; }

		word_end_reference &operator &= (word v) 
			{ return *this = (word) *this & v; }

		word_end_reference &operator |= (word_end_reference &a) 
			{ return *this = (word) *this | (word) a; }

		word_end_reference &operator |= (word v) 
			{ return *this = (word) *this | v; }

		word_end_reference &operator ^= (word_end_reference &a) 
			{ return *this = ((word) *this) ^ ((word) a); }

		word_end_reference &operator ^= (word v) 
			{ return *this = ((word) *this) ^ v; }

		operator word (void) const
		{
			typename Endianness::word_pair w, m;

			m.parts.low = _mask;
			m.parts.high = 0ULL;

			m.full = Endianness::shift_right (m.full, _shift);

			w.parts.low = *_pos;

			if (m.parts.high != 0)
				w.parts.high = _pos[1];
			else
				w.parts.high = 0ULL;

			w.full = Endianness::shift_left (w.full, _shift);

			return w.parts.low & _mask;
		}

	private:
		friend class word_iterator;
		friend class const_word_iterator;
		friend class const_word_reference;

		typename Iterator::word_iterator _pos;
		uint8 _shift;
		word _mask;
	};

	class word_iterator
	{
	public:
		typedef std::random_access_iterator_tag iterator_category;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference reference;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference *pointer;
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type value_type;
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::difference_type difference_type;
		typedef typename Iterator::Endianness Endianness;

		word_iterator () {}
		word_iterator (typename Iterator::word_iterator pos, uint8 shift)
			: _ref (pos, shift) {}
		word_iterator (const word_iterator &i)
			: _ref (i._ref._pos, i._ref._shift) {}

		word_iterator &operator = (const word_iterator &i)
		{
			_ref._pos = i._ref._pos;
			_ref._m.full = i._ref._m.full;
			_ref._shift = i._ref._shift;
			return *this;
		}

		word_iterator &operator ++ () 
		{
			++_ref._pos;
			return *this;
		}

		word_iterator operator ++ (int) 
		{
			word_iterator tmp (*this);
			++*this;
			return tmp;
		}

		word_iterator operator + (difference_type i) const
			{ return word_iterator (_ref._pos + i, _ref._shift); }

		word_iterator &operator += (difference_type i) 
		{
			_ref._pos += i;
			return *this;
		}

		word_iterator &operator -- () 
		{
			--_ref._pos;
			return *this;
		}

		word_iterator operator -- (int) 
		{
			word_iterator tmp (*this);
			--*this;
			return tmp;
		}

		word_iterator operator - (word_iterator i) const
			{ return _ref._pos - i._ref._pos; }

		word_iterator operator - (const_word_iterator i) const
			{ return _ref._pos - i._pos; }

		word_iterator operator - (difference_type i) const
			{ return *this + -i; }

		word_iterator &operator -= (difference_type i) 
			{ return *this += -i; }

		difference_type operator - (word_iterator &i) const 
			{ return _ref._pos - i._ref._pos; }

		typename BitSubvector<Iterator, ConstIterator>::word_reference operator [] (long i) 
			{ return *(*this + i); }

		typename BitSubvector<Iterator, ConstIterator>::word_reference operator * () 
			{ return _ref; }

		typename BitSubvector<Iterator, ConstIterator>::word_reference *operator -> () 
			{ return &_ref; }

		value_type operator * () const 
			{ return _ref; }

		bool operator == (const word_iterator &c) const 
			{ return (_ref._pos == c._ref._pos); }

		bool operator != (const word_iterator &c) const 
			{ return (_ref._pos != c._ref._pos); }

	private:
		friend class const_word_iterator;

		word_reference _ref;
	};

	class const_word_iterator
	{
	public:
		typedef std::random_access_iterator_tag iterator_category;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference reference;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference *pointer;
		typedef typename std::iterator_traits<typename Iterator::const_word_iterator>::value_type value_type;
		typedef typename std::iterator_traits<typename Iterator::const_word_iterator>::difference_type difference_type;
		typedef typename Iterator::Endianness Endianness;

		const_word_iterator () {}
		const_word_iterator (typename Iterator::const_word_iterator pos, uint8 shift)
			: _pos (pos), _shift (shift) {}
		const_word_iterator (const const_word_iterator  &i)
			: _pos (i._pos), _shift (i._shift) {}
		const_word_iterator (const word_iterator  &i)
			: _pos (i._ref._pos), _shift (i._ref._shift) {}

		const_word_iterator &operator = (const const_word_iterator  &i)
		{
			_pos = i._pos;
			_shift = i._shift;
			return *this;
		}

		const_word_iterator &operator = (const word_iterator  &i)
		{
			_pos = i._ref._pos;
			_shift = i._ref._shift;
			return *this;
		}

		const_word_iterator &operator ++ () 
		{
			++_pos;
			return *this;
		}

		const_word_iterator operator ++ (int) 
		{
			const_word_iterator tmp (*this);
			++*this;
			return tmp;
		}

		const_word_iterator operator + (difference_type i) const
			{ return const_word_iterator (_pos + i, _shift); }

		const_word_iterator &operator += (difference_type i) 
		{
			_pos += i;
			return *this;
		}

		const_word_iterator &operator -- () 
		{
			--_pos;
			return *this;
		}

		const_word_iterator operator -- (int) 
		{
			const_word_iterator tmp (*this);
			--*this;
			return tmp;
		}

		const_word_iterator operator - (difference_type i) const
			{ return *this + -i; }

		const_word_iterator &operator -= (difference_type i) 
			{ return *this += -i; }

		difference_type operator - (const_word_iterator i) const 
			{ return _pos - i._pos; }

		value_type operator [] (long i) 
			{ return *(*this + i); }

		value_type operator * ()
		{
			typename Endianness::word_pair w;
			w.parts.low = *_pos;
			w.parts.high = _pos[1];
			w.full = Endianness::shift_left (w.full, _shift);
			return w.parts.low;
		}

		bool operator == (const word_iterator &c) const 
			{ return (_pos == c._ref._pos); }

		bool operator == (const const_word_iterator &c) const 
			{ return (_pos == c._pos); }

		bool operator != (const word_iterator &c) const 
			{ return (_pos != c._ref._pos); }

		bool operator != (const const_word_iterator &c) const 
			{ return (_pos != c._pos); }

	private:
		friend class word_iterator;

		/* Semantics: see word_reference::_pos */
		typename Iterator::const_word_iterator _pos;
		uint8 _shift;
	};

}; // template <class Iterator> class BitSubvector

} // namespace LELA

namespace std {

template<class Iterator, class ConstIterator>
void swap (LELA::BitSubvector<Iterator, ConstIterator> &x, LELA::BitSubvector<Iterator, ConstIterator> &y)
{
	lela_check (x.size () == y.size ());

	typename LELA::BitSubvector<Iterator, ConstIterator>::word_iterator i_x, i_y;

	for (i_x = x.word_begin (), i_y = y.word_begin (); i_x != x.word_end () && i_y != y.word_end (); ++i_x, ++i_y) {
		typename Iterator::word_type word = *i_x;
		*i_x = *i_y;
		*i_y = word;
	}

	typename Iterator::word_type word = x.back_word ();
	x.back_word () = y.back_word ();
	y.back_word () = word;
}

} // namespace std

#endif // __LELA_VECTOR_BIT_SUBVECTOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
