/* linbox/vector/bit-vector.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_bit_vector_H
#define __LINBOX_bit_vector_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-iterator.h"

namespace LinBox
{

/** A vector of boolean 0-1 values, stored compactly to save space. 
 *
 * BitVector provides an additional iterator, word_iterator, that gives
 * the bits in compact 32-bit words, so that vector operations may be done in
 * parallel. It is similar to the STL bit_vector except that it provides the
 * aforementioned additional iterator.
 \ingroup vector
 */

template <class _Endianness = DefaultEndianness>
class BitVector
{
    public:
	typedef bool        value_type;
	typedef size_t      size_type;
	typedef long         difference_type;
	typedef typename _Endianness::word word_type;
	typedef typename std::vector<word_type>::iterator               word_iterator;
	typedef typename std::vector<word_type>::const_iterator         const_word_iterator;
	typedef typename std::vector<word_type>::reverse_iterator       reverse_word_iterator;
	typedef typename std::vector<word_type>::const_reverse_iterator const_reverse_word_iterator;

	typedef _Endianness Endianness;

	BitVector () {}
	BitVector (std::vector<bool> &v)
		{ *this = v; }
	BitVector (std::vector<word_type> &v)
		: _v (v), _size (_v.size () * WordTraits<word_type>::bits) {}
	BitVector (size_t n, bool val = false)
		{ resize (n, val); }
		
	// Copy constructor
	BitVector (const BitVector<> &v) 
		: _v (v._v), _size (v._size) {}

	~BitVector () {}

	// Reference

	typedef BitVectorReference<word_iterator, Endianness> reference;
	typedef BitVectorConstReference<const_word_iterator, Endianness> const_reference;

	// Iterators

	typedef BitVectorIterator<word_iterator, const_word_iterator, Endianness> iterator;
	typedef BitVectorConstIterator<const_word_iterator, Endianness> const_iterator;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef iterator    pointer;
	typedef const_iterator    const_pointer;

	inline iterator                    begin      (void);
	inline const_iterator              begin      (void) const;
	inline iterator                    end        (void);
	inline const_iterator              end        (void) const;

	inline reverse_iterator            rbegin     (void);
	inline const_reverse_iterator      rbegin     (void) const;
	inline reverse_iterator            rend       (void);
	inline const_reverse_iterator      rend       (void) const;

	inline word_iterator               wordBegin  (void)       { return _v.begin (); }
	inline const_word_iterator         wordBegin  (void) const { return _v.begin (); }
	inline word_iterator               wordEnd    (void)       { return _v.end (); }
	inline const_word_iterator         wordEnd    (void) const { return _v.end (); }

	inline reverse_word_iterator       wordRbegin (void)       { return _v.rbegin (); }
	inline const_reverse_word_iterator wordRbegin (void) const { return _v.rbegin (); }
	inline reverse_word_iterator       wordRend   (void)       { return _v.rend (); }
	inline const_reverse_word_iterator wordRend   (void) const { return _v.rend (); }

	// Element access

	inline reference       operator[] (size_type n);
	inline const_reference operator[] (size_type n) const;

	inline reference at       (size_type n);
	inline const_reference at (size_type n) const;

	inline reference       front (void);
	inline const_reference front (void) const;
	inline reference       back  (void);
	inline const_reference back  (void) const;

	// Word access

	inline word_type       &front_word (void)       { return _v.front (); }
	inline const word_type &front_word (void) const { return _v.front (); }
	inline word_type       &back_word  (void)       { return _v.back (); }
	inline const word_type &back_word  (void) const { return _v.back (); }

	// Appending to end of vector

	inline void push_back (bool v);

	// Push a whole word onto the back of the bit-vector
	// WARNING: Does not change size! Call resize afterwards to set the correct size.
	inline void push_word_back (word_type w)
		{ _v.push_back (w); }

	// Inserting a word into the vector
	inline void insertWord (word_iterator position, word_type word)
		{ _v.insert (position, word); }

	// Erasing a word from the vector
	inline void eraseWord (word_iterator position)
		{ _v.erase (position); }

	template<class Container>
	inline BitVector &operator = (const Container& x);

	inline void resize (size_type new_size, bool val = false);

	inline void clear () { _v.clear (); _size = 0; }

	inline size_type size      (void) const { return _size;            }
	inline bool      empty     (void) const { return _v.empty ();      }

	inline size_type word_size (void) const { return _v.size ();       }

	// FIXME: This should be _v.max_size (), not _v.size (), or am I mistaken?
	inline size_type max_size  (void) const { return _v.size  () * WordTraits<word_type>::bits; }

	inline bool operator == (const BitVector &v) const;

	// Used for hybrid vectors to make sure _size is correctly set after setting the vector up
	void fix_size (size_type size_mod_word) {
		if (size_mod_word & WordTraits<word_type>::pos_mask)
			_size = ((_v.size () - 1) << WordTraits<word_type>::logof_size) + (size_mod_word & WordTraits<word_type>::pos_mask);
		else
			_size = (_v.size () << WordTraits<word_type>::logof_size);
	}

    protected:

	std::vector<word_type> _v;
	size_t              _size;

}; // template <class Vector> class ReverseVector

// Vector traits for BitVector wrapper
template <class Endianness>
struct GF2VectorTraits<BitVector<Endianness> >
{ 
	typedef BitVector<Endianness> VectorType;
	typedef VectorCategories::DenseZeroOneVectorTag VectorCategory; 
};

template <class Endianness>
struct GF2VectorTraits<const BitVector<Endianness> >
{ 
	typedef BitVector<Endianness> VectorType;
	typedef VectorCategories::DenseZeroOneVectorTag VectorCategory; 
};

} // namespace LinBox

#include "linbox/vector/bit-vector.inl"

#endif // __LINBOX_bit_vector_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
