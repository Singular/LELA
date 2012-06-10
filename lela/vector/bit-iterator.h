/* lela/vector/bit-iterator.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Reference- and iterator-clases for bit-vectors
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_BIT_ITERATOR_H
#define __LELA_VECTOR_BIT_ITERATOR_H

#include <stdexcept>
#include <vector>

#include "lela/lela-config.h"
#include "lela/integer.h"

namespace LELA
{

#ifndef __LELA_BITVECTOR_WORD_TYPE
#  if __LELA_SIZEOF_LONG == 8
#    define __LELA_BITVECTOR_WORD_TYPE uint64
#  else
#    define __LELA_BITVECTOR_WORD_TYPE uint32
#  endif
#endif

/** Structure definining information on word-sizes
 */
template <class _Word>
struct WordTraits {
	typedef _Word Word;
	static const unsigned int bits;
	static const unsigned int logof_size;
	static const unsigned int pos_mask;
	static const Word all_ones;

	static inline bool ParallelParity (Word t);
};

/** Specialisation for 8-bit words */
template <>
struct WordTraits<uint8> {
	typedef uint8 Word;
	typedef uint16 DoubleWord;
	static const unsigned int bits = 8;
	static const unsigned int logof_size = 3;
	static const unsigned int pos_mask = 0x07;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};

/** Specialisation for int */
#if __LELA_SIZEOF_INT == 4
template <>
struct WordTraits<unsigned int> {
	typedef unsigned int Word;
	typedef uint64 DoubleWord;
	static const unsigned int bits = 32;
	static const unsigned int logof_size = 5;
	static const unsigned int pos_mask = 0x1F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};
#elif __LELA_SIZEOF_INT == 8
template <>
struct WordTraits<unsigned int> {
	typedef unsigned int Word;
	typedef uint128 DoubleWord;
	static const unsigned int bits = 64;
	static const unsigned int logof_size = 6;
	static const unsigned int pos_mask = 0x3F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 32);
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};
#endif // __LELA_SIZEOF_INT

/** Specialisation for long */
#if __LELA_SIZEOF_LONG == 4
template <>
struct WordTraits<unsigned long> {
	typedef unsigned long Word;
	typedef uint64 DoubleWord;
	static const unsigned int bits = 32;
	static const unsigned int logof_size = 5;
	static const unsigned int pos_mask = 0x1F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};
#elif __LELA_SIZEOF_LONG == 8
template <>
struct WordTraits<unsigned long> {
	typedef unsigned long Word;
	typedef uint128 DoubleWord;
	static const unsigned int bits = 64;
	static const unsigned int logof_size = 6;
	static const unsigned int pos_mask = 0x3F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 32);
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};
#endif

/** Specialisation for unsigned long long */
#if __LELA_SIZEOF_LONG_LONG == 8
template <>
struct WordTraits<unsigned long long> {
	typedef unsigned long long Word;
	typedef uint128 DoubleWord;
	static const unsigned int bits = 64;
	static const unsigned int logof_size = 6;
	static const unsigned int pos_mask = 0x3F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 32);
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};
#elif __LELA_SIZEOF_LONG_LONG == 16
template <>
struct WordTraits<unsigned long long> {
	typedef unsigned long long Word;
	typedef uint256 DoubleWord;
	static const unsigned int bits = 128;
	static const unsigned int logof_size = 7;
	static const unsigned int pos_mask = 0x7F;
	static const Word all_ones = static_cast<const Word> (-1);

	static inline bool ParallelParity (Word t) {
		t ^= (t >> 64);
		t ^= (t >> 32);
		t ^= (t >> 16);
		t ^= (t >> 8);
		t ^= (t >> 4);
		t &= 0xf;
		return bool( (0x6996 >> t) & 0x1);
	}
};
#endif

// Generic routines for big endian word-order

// Big endian version: Position zero is the highest bit on the word
// (so it prints the right way around if you convert the word to
// binary)
template <class _word>
class BigEndian {
public:
	typedef _word word;

	union word_pair {
		struct {
#ifdef __LELA_HAVE_LITTLE_ENDIAN
			word high;
			word low;
#elif __LELA_HAVE_BIG_ENDIAN
			word low;
			word high;
#endif // __LELA_HAVE_LITTLE_ENDIAN
		} parts;

		typename WordTraits<word>::DoubleWord full;

		void init () { parts.low = parts.high = 0; }

		void advance () { parts.low = parts.high; }

		template <class iterator>
		void read_next (iterator i)
			{ parts.low = parts.high; parts.high = *i; }
	};

	// Constant representing a one in the position zero in the word
	static const word e_0 = 1ULL << (WordTraits<word>::bits - 1);

	// Shift the given word pos positions to the right
	template <class Word>
	static inline Word shift_right (Word w, uint8 pos) { return w >> pos; }

	// Shift the given word pos positions to the left
	template <class Word>
	static inline Word shift_left (Word w, uint8 pos) { return w << pos; }

	// Return a word with all positions from pos onwards set to one and the rest set to zero
	static inline word mask_right (uint8 pos) { return (((e_0 >> pos) << 1) - 1); }

	// Return a word with all positions up to and not including pos set to one and the rest set to zero
	static inline word mask_left (uint8 pos) { return ~(((e_0 >> pos) << 1) - 1); }

	// Return a word with only the bit at the lowest position of the input word set
	static inline word first_position (word w) {
		word v = e_0;

		while (v && !(w & v)) v >>= 1;

		return v;
	}

	// Return e_j
	static inline word e_j (uint8 j) { return shift_right (e_0, j); }
};

// Little endian version: position zero is in the lowest position in
// the word (so e_k is 2^k)
template <class _word>
class LittleEndian {
public:
	typedef _word word;

	union word_pair {
		struct {
#ifdef __LELA_HAVE_LITTLE_ENDIAN
			word low;
			word high;
#elif __LELA_HAVE_BIG_ENDIAN
			word high;
			word low;
#endif // __LELA_HAVE_LITTLE_ENDIAN
		} parts;

		typename WordTraits<word>::DoubleWord full;

		void init () { parts.low = parts.high = 0; }

		void advance () { parts.low = parts.high; }

		template <class iterator>
		void read_next (iterator i)
			{ parts.low = parts.high; parts.high = *i; }
	};

	// Constant representing a one in the position zero in the word
	static const word e_0 = 1ULL;

	// Shift the given word pos positions to the right
	template <class Word>
	static inline Word shift_right (Word w, uint8 pos) { return w << pos; }

	// Shift the given word pos positions to the left
	template <class Word>
	static inline Word shift_left (Word w, uint8 pos) { return w >> pos; }

	// Return a word with all positions from pos onwards set to one and the rest set to zero
	static inline word mask_right (uint8 pos) { return ~((e_0 << pos) - 1); }

	// Return a word with all positions up to and not including pos set to one and the rest set to zero
	static inline word mask_left (uint8 pos) { return ((e_0 << pos) - 1); }

	// Return a word with only the bit at the lowest position of the input word set
	static inline word first_position (word w) { return ((w ^ (w - 1)) >> 1) + 1; }

	// Return e_j
	static inline word e_j (uint8 j) { return shift_right (e_0, j); }
};

// Define default endianness

// We use the endianness which is given by the version of libm4ri if
// it is available, since in any case it must be compatible with what
// libm4ri uses. Otherwise we use little endian by default.
#ifdef __LELA_HAVE_M4RI
#  ifdef __LELA_HAVE_M4RI_GE_20110601
#    define DefaultEndianness LittleEndian
#  else // !__LELA_HAVE_M4RI_GE_20110601
#    define DefaultEndianness BigEndian
#  endif // __LELA_HAVE_M4RI_GE_20110601
#else // !__LELA_HAVE_M4RI
#  define DefaultEndianness LittleEndian
#endif // __LELA_HAVE_M4RI

template <class word_iterator, class const_word_iterator, class Endianness>
class BitVectorIterator;

template <class const_word_iterator, class Endianness>
class BitVectorConstIterator;

template <class word_iterator, class _Endianness = DefaultEndianness<typename std::iterator_traits<word_iterator>::value_type> >
class BitVectorReference
{
    public:
	typedef _Endianness Endianness;

	BitVectorReference (word_iterator word, uint8 position)
		: _word (word), _pos (position) {}

	BitVectorReference (const BitVectorReference &v)
		: _word (v._word), _pos (v._pos) {}

	~BitVectorReference () {}

	BitVectorReference &operator = (BitVectorReference &a) 
		{ return *this = (bool) a; }

	BitVectorReference &operator = (bool v) 
	{ 
		*_word = v ? (*_word | (Endianness::e_j (_pos))) : (*_word & ~Endianness::e_j (_pos));
		return *this;
	}

	BitVectorReference &operator &= (BitVectorReference &a) 
		{ *_word &= ~Endianness::e_j (_pos) | Endianness::shift_right (a.get_bit (), _pos - a._pos); return *this; }

	BitVectorReference &operator &= (bool v) 
		{ if (!v) *_word &= ~Endianness::e_j (_pos); return *this; }

	BitVectorReference &operator |= (BitVectorReference &a) 
		{ *_word |= Endianness::shift_right (a.get_bit (), _pos - a._pos); return *this; }

	BitVectorReference &operator |= (bool v) 
		{ if (v) *_word |= Endianness::e_j (_pos); return *this; }

	BitVectorReference &operator ^= (BitVectorReference &a) 
		{ *_word ^= Endianness::shift_right (a.get_bit (), _pos - a._pos); return *this; }

	BitVectorReference &operator ^= (bool v) 
		{ if (v) *_word ^= Endianness::e_j (_pos); return *this; }

	operator bool (void) const
		{ return *_word & Endianness::e_j (_pos); return *this; }

    private:
	template <class _word_iterator, class const_word_iterator, class Endianness>
	friend class BitVectorIterator;

	template <class const_word_iterator, class Endianness>
	friend class BitVectorConstIterator;

	template <class const_word_iterator, class Endianness>
	friend class BitVectorConstReference;

	typename std::iterator_traits<word_iterator>::value_type neg_mask_word (void) { return *_word & ~Endianness::e_j (_pos); }
	typename std::iterator_traits<word_iterator>::value_type get_bit ()           { return *_word & Endianness::e_j (_pos); }

	word_iterator _word;
	uint8         _pos;
};

template <class word_iterator, class Endianness>
inline std::istream &operator >> (std::istream &is, BitVectorReference<word_iterator, Endianness> &a) 
	{ bool v; is >> v; a = v; return is; }

template <class word_iterator, class Endianness>
inline std::ostream &operator << (std::ostream &os, BitVectorReference<word_iterator, Endianness> &a) 
	{ os << bool (a); return os; }

template <class const_word_iterator, class _Endianness = DefaultEndianness<typename std::iterator_traits<const_word_iterator>::value_type> >
class BitVectorConstReference
{
    public:
	typedef _Endianness Endianness;

	template <class Iterator>
	BitVectorConstReference (BitVectorReference<Iterator, Endianness> r)
		: _word (r._word), _pos (r._pos) {}

	template <class Iterator>
	BitVectorConstReference (Iterator word, uint8 position)
		: _word (word), _pos (position) {}

	BitVectorConstReference (const_word_iterator word, uint8 position)
		: _word (word), _pos (position) {}

	BitVectorConstReference (const BitVectorConstReference &v)
		: _word (v._word), _pos (v._pos) {}

	~BitVectorConstReference () {}

	operator bool (void) const
		{ return (*_word & Endianness::e_j (_pos)) != 0; }

    private:
	friend class BitVectorConstIterator<const_word_iterator, Endianness>;

	const_word_iterator _word;
	uint8               _pos;
};

template<class const_word_iterator, class Endianness>
inline std::ostream &operator << (std::ostream &os, BitVectorConstReference<const_word_iterator, Endianness> &a) 
	{ os << bool (a); return os; }

// class BitVectorIterator : public std::iterator <std::random_access_iterator_tag, bool>
template <class _word_iterator, class _const_word_iterator, class _Endianness = DefaultEndianness<typename std::iterator_traits<_word_iterator>::value_type> >
class BitVectorIterator : public std::_Bit_iterator
{
    public:

	typedef _word_iterator word_iterator;
	typedef _const_word_iterator const_word_iterator;
	typedef _Endianness Endianness;

	typedef std::random_access_iterator_tag iterator_category;
	typedef BitVectorReference<word_iterator, Endianness> reference;
	typedef BitVectorConstReference<const_word_iterator, Endianness> const_reference;
	typedef bool *pointer;
	typedef bool value_type;
	typedef long difference_type;
	typedef size_t size_type;
	typedef typename std::iterator_traits<const_word_iterator>::value_type word_type;

	BitVectorIterator () : _ref (word_iterator (), 0UL) {}
	BitVectorIterator (word_iterator word, uint8 position) : _ref (word, position) {}
	BitVectorIterator (const BitVectorIterator &i) : _ref (i._ref._word, i._ref._pos) {}

	BitVectorIterator &operator = (const BitVectorIterator &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	BitVectorIterator &operator ++ () 
	{
		if (++_ref._pos > (WordTraits<word_type>::bits - 1)) {
			++_ref._word;
			_ref._pos = 0UL;
		}

		return *this;
	}

	BitVectorIterator operator ++ (int) 
	{
		BitVectorIterator tmp (*this);
		++*this;
		return tmp;
	}

	BitVectorIterator operator + (difference_type i) const
	{
		word_iterator new_word = _ref._word + (i >> WordTraits<word_type>::logof_size);
		uint8 new_pos = _ref._pos + (i & WordTraits<word_type>::pos_mask);

		new_word += new_pos >> WordTraits<word_type>::logof_size;
		new_pos &= WordTraits<word_type>::pos_mask;

		return BitVectorIterator (new_word, new_pos);
	}

	BitVectorIterator &operator += (difference_type i) 
	{
		_ref._word += i >> WordTraits<word_type>::logof_size;
		_ref._pos  += i & WordTraits<word_type>::pos_mask;
		_ref._word += _ref._pos >> WordTraits<word_type>::logof_size;
		_ref._pos  &= WordTraits<word_type>::pos_mask;
		return *this;
	}

	BitVectorIterator &operator -- () 
	{
		if (--_ref._pos > (WordTraits<word_type>::bits - 1)) {
			--_ref._word;
			_ref._pos = (WordTraits<word_type>::bits - 1);
		}

		return *this;
	}

	BitVectorIterator operator -- (int) 
	{
		BitVectorIterator tmp (*this);
		--*this;
		return tmp;
	}

	BitVectorIterator operator - (difference_type i) const
		{ return *this + -i; }

	BitVectorIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (BitVectorIterator i) const 
		{ return (_ref._word - i._ref._word) * WordTraits<word_type>::bits + (_ref._pos - i._ref._pos); }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ return _ref; }

	const_reference operator * () const 
		{ return _ref; }

	bool operator == (const BitVectorIterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator != (const BitVectorIterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

	word_iterator word () const { return _ref._word; }
	uint8 pos () const { return _ref._pos; }

    private:
	friend class BitVectorConstIterator<const_word_iterator, Endianness>;

	reference _ref;
};

template <class _const_word_iterator, class _Endianness = DefaultEndianness<typename std::iterator_traits<_const_word_iterator>::value_type> >
class BitVectorConstIterator : public std::iterator <std::random_access_iterator_tag, bool>
{
    public:

	typedef _const_word_iterator word_iterator;
	typedef _const_word_iterator const_word_iterator;
	typedef _Endianness Endianness;

	typedef std::random_access_iterator_tag iterator_category;
	typedef BitVectorConstReference<const_word_iterator, Endianness> reference;
	typedef reference const_reference;
	typedef const bool *pointer;
	typedef bool value_type;
	typedef long difference_type;
	typedef size_t size_type;
	typedef typename std::iterator_traits<const_word_iterator>::value_type word_type;

	BitVectorConstIterator () : _ref (const_word_iterator (), 0UL) {}
	BitVectorConstIterator (const_word_iterator word, uint8 position) : _ref (word, position) {}
	BitVectorConstIterator (const BitVectorConstIterator &i) : _ref (i._ref._word, i._ref._pos) {}
	template <class Iterator>
	BitVectorConstIterator (const BitVectorIterator<Iterator, const_word_iterator, Endianness> &i) : _ref (i._ref._word, i._ref._pos) {}

	BitVectorConstIterator &operator = (const BitVectorConstIterator &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	template <class Iterator>
	BitVectorConstIterator &operator = (const BitVectorIterator<Iterator, const_word_iterator, Endianness> &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	BitVectorConstIterator &operator ++ () 
	{
		if (++_ref._pos > WordTraits<word_type>::bits - 1) {
			++_ref._word;
			_ref._pos = 0UL;
		}

		return *this;
	}

	BitVectorConstIterator operator ++ (int) 
	{
		BitVectorConstIterator tmp (*this);
		++*this;
		return tmp;
	}

	BitVectorConstIterator operator + (long i) const
	{
		const_word_iterator new_word = _ref._word + (i >> WordTraits<word_type>::logof_size);
		uint8 new_pos = _ref._pos + (i & WordTraits<word_type>::pos_mask);

		new_word += new_pos >> WordTraits<word_type>::logof_size;
		new_pos &= WordTraits<word_type>::pos_mask;

		return BitVectorConstIterator (new_word, new_pos);
	}

	BitVectorConstIterator &operator += (long i) 
	{
		_ref._word += i >> WordTraits<word_type>::logof_size;
		_ref._pos  += i & WordTraits<word_type>::pos_mask;
		_ref._word += _ref._pos >> WordTraits<word_type>::logof_size;
		_ref._pos  &= WordTraits<word_type>::pos_mask;
		return *this;
	}

	BitVectorConstIterator &operator -- () 
	{
		if (--_ref._pos > WordTraits<word_type>::bits - 1) {
			--_ref._word;
			_ref._pos = WordTraits<word_type>::bits - 1;
		}

		return *this;
	}

	BitVectorConstIterator operator -- (int) 
	{
		BitVectorConstIterator tmp (*this);
		--*this;
		return tmp;
	}

	BitVectorConstIterator operator - (difference_type i) const 
		{ return *this + -i; }

	BitVectorConstIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (BitVectorConstIterator i) const 
		{ return (_ref._word - i._ref._word) * WordTraits<word_type>::bits + (_ref._pos - i._ref._pos); }

	reference operator [] (difference_type i) const
		{ return *(*this + i); }

	reference operator * () const
		{ return _ref; }

	bool operator == (const BitVectorConstIterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	template <class Iterator>
	bool operator == (const BitVectorIterator<Iterator, const_word_iterator, Endianness> &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator != (const BitVectorConstIterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

	template <class Iterator>
	bool operator != (const BitVectorIterator<Iterator, const_word_iterator, Endianness> &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

	const_word_iterator word () const { return _ref._word; }
	uint8 pos () const { return _ref._pos; }

    private:

	const_reference _ref;
};

} // namespace LELA

namespace std 
{

template <class word_iterator, class Endianness>
void swap (LELA::BitVectorReference<word_iterator, Endianness> x, LELA::BitVectorReference<word_iterator, Endianness> y)
	{ bool t = x; x = y; y = t; }

} // namespace std

#endif // __LELA_VECTOR_BIT_ITERATOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
