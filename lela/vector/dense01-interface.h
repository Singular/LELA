/* lela/vector/dense01-interface.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Prototype for dense vectors. This file is provided for documentation only.
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 *
 */

#ifndef __LELA_VECTOR_DENSE01_INTERFACE_H
#define __LELA_VECTOR_DENSE01_INTERFACE_H

#include <iterator>

#include "linbox/vector/traits.h"
#include "linbox/vector/bit-iterator.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector.h"

namespace LELA
{

/** Dense 0-1 vector-interface
 *
 * This class is a prototype of a vector meeting the dense 0-1
 * interface in LELA. It can be used as a starting-point for those
 * wishing to implement their own dense 0-1 vector-types. It cannot be
 * used on its own.
 *
 * When feasible, suggested implementations of methods are provided.
 *
 * \ingroup vector
 */
class Dense01VectorInterface
{
public:
	/// @name Dense 0-1-specific types
	//@{

	/** Normal word-iterator
	 *
	 * This iterator iterates through the words of a dense 0-1
	 * vector, allowing fast parallel operations on entire words
	 * at a time.
	 *
	 * This iterator must stop at the *penultimate* word of the
	 * vector, not the last word! This convention makes
	 * bit-subvectors much faster and easier to implement.
	 *
	 * This must be a random access iterator, i.e. providing
	 * the operators ++, --, +, +=, - (both with an integral type
	 * and with another iterator), -=, ==, and !=.
	 *
	 * It must permit comparison with const_iterator as well as
	 * iterator.
	 *
	 * An ordinary C++ pointer satisfies all of these
	 * requirements and could in principle be used here.
	 */
	class word_iterator;

	/** Const iterator
	 *
	 * It must support the same operators as with the normal
	 * iterator. The same convention also applies.
	 *
	 * It must permit comparison with iterator as well as
	 * const_iterator and it must permit assignment from and
	 * construction from iterator.
	 *
	 * An ordinary C++ const pointer satisfies all of these
	 * requirements and could in principle be used here.
	 */
	class const_word_iterator;

	/// The type of a word, typically uint64
	typedef std::iterator_traits<word_iterator>::value_type word_type;

	/// A reference to the last word of the vector
	typedef word_type &back_word_reference;

	/// A const reference to the last word of the vector
	typedef const word_type &const_back_word_reference;

	/// Endianness: whether words are stored with the e_0 being
	/// the most significant bit (i.e. big-endian) or with the
	/// least significant bit (i.e. little-endian). It is
	/// essential that the endianness be applied consistently
	/// throughout the program: mixing different conventions will
	/// have unpredictable results!
	typedef LittleEndian<word_type> Endianness;

	//@}

	/// @name STL-defined types
	//@{

	/** Normal iterator
	 *
	 * This iterates through the entries of the vector one at a
	 * time.
	 *
	 * This must be a random access iterator, i.e. providing
	 * the operators ++, --, +, +=, - (both with an integral type
	 * and with another iterator), -=, ==, and !=.
	 *
	 * It must permit comparison with const_iterator as well as
	 * iterator.
	 *
	 * The class BitVectorIterator is provided for this purpose.
	 */
	typedef BitVectorIterator<word_iterator, const_word_iterator, Endianness> iterator;

	/** Const iterator
	 *
	 * This iterates through the entries of the vector one at a
	 * time.
	 *
	 * It must support the same operators as with the normal
	 * iterator.
	 *
	 * It must permit comparison with iterator as well as
	 * const_iterator and it must permit assignment from and
	 * construction from iterator.
	 *
	 * An ordinary C++ const pointer satisfies all of these
	 * requirements and could in principle be used here.
	 *
	 * The class BitVectorIterator is provided for this purpose.
	 */
	typedef BitVectorIterator<const_word_iterator, const_word_iterator, Endianness> iterator;

	/// Type of the elements of the vector
	typedef typename std::iterator_traits<iterator>::value_type value_type;
    
	/// Type of indices into the vector and of the size of the vector
	typedef size_t size_type;

	/// Type of the difference (i.e. distance) between two iterators
	typedef typename std::iterator_traits<iterator>::difference_type difference_type;

	/// Type of a pointer to an element in the vector
	typedef typename std::iterator_traits<iterator>::pointer pointer;

	/// Type of a reference to an element in the vector
	typedef typename std::iterator_traits<iterator>::reference reference;

	/// Type of a const reference
	typedef typename std::iterator_traits<const_iterator>::reference const_reference; 

	//@}

	/// @name LELA-defined types
	//@{

	/// Representation-type: always dense
	typedef VectorRepresentationTypes::Dense01 RepresentationType;

	/// Storage-type: see VectorStorageTypes
	/// Use Real if the vector represents a live array in memory, Generic if not
	typedef VectorStorageTypes::Generic StorageType;

	/// Container-type: this should be the type a user would use
	/// in generic code to construct another vector of the same
	/// type, e.g., if the user wishes to construct a copy.
	typedef BitVector<Endianness> ContainerType;

	/// Subvector-type: the type the user should use to construct
	/// a subvector of this vector.
	///
	/// The subvector-type must support a construct of the form
	/// Subvector (v, start_index, end_index)
	/// where v is the parent-vector and start_index and end_index
	/// are indices of the beginning and end of the subvector.
	///
	/// Normally it suffices to use Subvector for this type
	/// as given here.
	typedef BitSubvector<iterator, const_iterator> SubvectorType;

	/// Const subvector-type: to be used when the parent is const
	///
	/// Normally it suffices to use Subvector for this type
	/// as given here.
	typedef BitSubvector<const_iterator, const_iterator> ConstSubvectorType;

	/// Aligned subvector-type: similar to subvector, but may
	/// restrict the start and possibly end indices to be integer
	/// multiples of the constant align and should then provide
	/// better performance. In nearly all cases, it should be
	/// identical to Subvector.
	///
	/// Normally it suffices to use Subvector for this type
	/// as given here.
	typedef BitSubvectorWordAligned<iterator, const_iterator> AlignedSubvectorType;

	/// Const version of AlignedSubvectorType
	///
	/// Normally it suffices to use Subvector for this type
	/// as given here.
	typedef Subvector<ConstIterator, ConstIterator> ConstAlignedSubvectorType;

	/// Alignment-value: aligned subvectors must have start and
	/// possibly ending indices which are integer multiples of
	/// this value. In nearly all cases -- and particularly when
	/// SubvectorType and AlignedSubvectorType are the same, it
	/// should be one.
	static const int align = 1;

	//@}

	/// @name STL vector-interface
	//@{

	/// @name Iterator-access
	//@{

	/// Obtain an iterator for the beginning of the vector
	iterator begin ();

	/// Obtain a const iterator for the beginning of the vector
	const_iterator begin  () const;

	/// Obtain an iterator for one past the end of the vector
	iterator end ();

	/// Obtain a const iterator for one past the end of the vector
	const_iterator end () const;

	//@}

	/// @name Array-style element-access
	//@{

	/// Obtain a reference to the element at position n
	inline reference operator[] (size_type n) { return *(begin () + n); }

	/// Obtain a const reference to the element at position n
	inline const_reference operator[] (size_type n) const { return *(begin () + n); }

	//@}

	/// @name Information on the vector
	//@{

	/// Return the number of entries (i.e. bits) of the vector
	inline size_type size () const { return end () - begin (); }

	/// Determine whether the vector is empty; return true if so and false otherwise
	inline bool empty () const { return end () == begin (); }

	//@}

	//@}

	/// @name Word vector-interface
	///
	/// This mirrors the STL-interface above but provides direct
	/// access to the words of the vector
	//@{

	/// @name Iterator-access
	//@{

	/// Obtain a word-iterator for the beginning of the vector
	word_iterator word_begin ();

	/// Obtain a const word-iterator for the beginning of the vector
	const_word_iterator word_begin () const;

	/// Obtain a word-iterator for the last word of the vector
	word_iterator word_end ();

	/// Obtain a const word-iterator for last word of the vector
	const_word_iterator word_end () const;

	//@}

	/// @name Accessing the last word
	//@{

	/// Obtain a reference to the last word of the vector
	back_word_reference back_word ();

	/// Obtain a const reference to the last word of the vector
	const_back_word_reference back_word () const;

	//@}

	/// @name Information on the vector
	//@{

	/// Return the number of words of the vector
	inline size_type word_size () const { return word_end () - word_begin () + 1; }

	//@}

	//@}
};

} // namespace LELA

#endif // __LELA_VECTOR_SPARSE_INTERFACE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
