/* lela/vector/bit-vector.h
 * Copyright 2003 Bradford Hovinen
 *
 * Dense vectors over GF(2)
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_BIT_VECTOR_H
#define __LELA_BIT_VECTOR_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "lela/vector/traits.h"
#include "lela/vector/bit-iterator.h"

namespace LELA
{

// Forward-declarations
template <class Iterator, class ConstIterator = Iterator> class BitSubvector;
template <class Iterator, class ConstIterator, class _Endianness = DefaultEndianness<typename std::iterator_traits<Iterator>::value_type> > class BitSubvectorWordAligned;

/** A vector of boolean 0-1 values, stored compactly to save space. 
 *
 * BitVector provides an additional iterator, word_iterator, that gives
 * the bits in compact 32-bit words, so that vector operations may be done in
 * parallel. It is similar to the STL bit_vector except that it provides the
 * aforementioned additional iterator.
 \ingroup vector
 *
 * Notes about semantics:
 *  - word_iterator is guaranteed to iterate up to the *penultimate* word, not
 *    the last word
 *  - back_word gives a reference to the last word. It is not necessarily
 *    *(word_end () - 1)
 *  - Therefore, code working on BitVectors (as well as BitSubvectors) should
 *    iterate through the word_iterator and then do the same operation on 
 *    back_word ().
 */

template <class _Endianness = DefaultEndianness<uint64> >
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

	// Traits

	typedef VectorRepresentationTypes::Dense01 RepresentationType; 
	typedef VectorStorageTypes::Real StorageType;
	typedef BitVector ContainerType;
	typedef BitSubvector<iterator, const_iterator> SubvectorType;
	typedef BitSubvector<const_iterator, const_iterator> ConstSubvectorType;
	typedef BitSubvectorWordAligned<word_iterator, const_word_iterator> AlignedSubvectorType;
	typedef BitSubvectorWordAligned<const_word_iterator, const_word_iterator> ConstAlignedSubvectorType;
	static const int align = WordTraits<word_type>::bits;

	inline iterator                    begin      (void)
		{ return iterator (_v.begin (), 0UL); }
	inline const_iterator              begin      (void) const
		{ return const_iterator (_v.begin (), 0UL); }
	inline iterator                    end        (void);
	inline const_iterator              end        (void) const;

	inline reverse_iterator            rbegin     (void)
		{ return reverse_iterator (end () - 1UL); }
	inline const_reverse_iterator      rbegin     (void) const
		{ return const_reverse_iterator (end () - 1UL); }
	inline reverse_iterator            rend       (void)
		{ return reverse_iterator (begin () - 1UL); }
	inline const_reverse_iterator      rend       (void) const
		{ return const_reverse_iterator (begin () - 1UL); }

	inline word_iterator               word_begin  (void)       { return _v.begin (); }
	inline const_word_iterator         word_begin  (void) const { return _v.begin (); }
	inline word_iterator               word_end    (void)       { return _v.empty () ? _v.end () : _v.end () - 1; }
	inline const_word_iterator         word_end    (void) const { return _v.empty () ? _v.end () : _v.end () - 1; }

	inline reverse_word_iterator       word_rbegin (void)       { return _v.rbegin (); }
	inline const_reverse_word_iterator word_rbegin (void) const { return _v.rbegin (); }
	inline reverse_word_iterator       word_rend   (void)       { return _v.rend (); }
	inline const_reverse_word_iterator word_rend   (void) const { return _v.rend (); }

	// Element access

	inline reference       operator[] (size_type n)
		{ return *(begin () + n); }
	inline const_reference operator[] (size_type n) const
		{ return *(begin () + n); }

	inline reference at       (size_type n);
	inline const_reference at (size_type n) const;

	inline reference       front (void)
		{ return reference (_v.begin (), 0UL); }
	inline const_reference front (void) const
		{ return const_reference (_v.begin (), 0UL); }
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

	inline void resize (size_type new_size, bool val = false)
		{ _v.resize ((new_size >> WordTraits<word_type>::logof_size) + ((new_size & WordTraits<word_type>::pos_mask) ? 1UL : 0UL), val ? WordTraits<word_type>::all_ones : 0UL); _size = new_size; }

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

} // namespace LELA

#include "lela/vector/bit-vector.tcc"

#endif // __LELA_BIT_VECTOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
