/* lela/vector/hybrid01-interface.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Prototype for hybrid 0-1 vectors. This file is provided for documentation only.
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_HYBRID01_INTERFACE_H
#define __LELA_VECTOR_HYBRID01_INTERFACE_H

#include <iterator>

#include "linbox/vector/traits.h"
#include "linbox/vector/sparse-subvector.h"
#include "linbox/vector/bit-iterator.h"

namespace LELA
{

/** Hybrid 0-1 vector interface
 *
 * This class is a prototype of a vector meeting the hybrid 0-1
 * interface in LELA. It can be used as a starting-point for those
 * wishing to implement their own hybrid 0-1 vector-types. It cannot
 * be used on its own.
 *
 * When feasible, suggested implementations of methods are provided.
 *
 * \ingroup vector
 */
class Hybrid01VectorInterface
{
public:
	/// @name Hybrid-specific types
	//@{

	/// The type of a block-index. Usually uint16 suffices
	typedef uint16 index_type;

	/// The type of a word, typically uint64
	typedef uint64 word_type;

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
	 * This must be at least an input-iterator, providing the
	 * increment-operator ++ as well as the comparison-operators
	 * == and !=. If it is a random access iterator,
	 * i.e. providing in addition the operators ++, --, +, +=, -
	 * (both with an integral type and with another iterator), -=,
	 * then certain operations, such as the construction of
	 * subvectors, may be much faster.
	 * 
	 * It must permit comparison with const_iterator as well as
	 * iterator.
	 * 
	 * An ordinary C++ pointer satisfies all of these
	 * requirements and could in principle be used here.
	 */
	class iterator;

	/** Const iterator
	 * 
	 * The same comments about operators apply as with the normal
	 * iterator.
	 * 
	 * It must permit comparison with iterator as well as
	 * const_iterator and it must permit assignment from and
	 * construction from iterator.
	 * 
	 * An ordinary C++ const pointer satisfies all of these
	 * requirements and could in principle be used here.
	 */
	class const_iterator;

	/// Type of the entries of the vector
	///
	/// This should in general be an STL-pair of an integral
	/// block-index-type and a word.
	typedef typename std::iterator_traits<iterator>::value_type value_type;
    
	/// Type of indices into the vector and of the size of the vector
	typedef size_t size_type;

	/// Type of the difference (i.e. distance) between two iterators
	typedef typename std::iterator_traits<iterator>::difference_type difference_type;

	/// Type of a pointer to an entry in the vector
	typedef typename std::iterator_traits<iterator>::pointer pointer;

	/// Type of a reference to an entry in the vector
	typedef typename std::iterator_traits<iterator>::reference reference;

	/// Type of a const reference
	typedef typename std::iterator_traits<const_iterator>::reference const_reference; 

	//@}

	/// @name LELA-defined types
	//@{

	/// Representation-type: always hybrid 0-1
	typedef VectorRepresentationTypes::Hybrid01 RepresentationType;

	/** Storage-type: see VectorStorageTypes
	 *
	 * Use Real if the vector represents a live array of
	 * index-word pairs in memory, Transformed if it transforms a
	 * pair of vectors into a vector of pairs (as with
	 * HybridVector), or Generic if none if these is true.
	 */
	typedef VectorStorageTypes::Generic StorageType;

	/// Container-type: this should be the type a user would use
	/// in generic code to construct another vector of the same
	/// type, e.g., if the user wishes to construct a copy.
	typedef HybridVector<typename Vector::word_type> ContainerType;

	/// Subvector-type: the type the user should use to construct
	/// a subvector of this vector.
	///
	/// The subvector-type must support a construct of the form
	/// Subvector (v, start_index, end_index)
	/// where v is the parent-vector and start_index and end_index
	/// are indices of the beginning and end of the subvector.
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	///
	/// Since SparseSubvector does not support non-const arbitrary
	/// (i.e. not word-aligned) subvectors of hybrid vectors, we
	/// set SubvectorType to be identical to
	/// ConstSubvectorType. It is then an error to construct an
	/// arbitrary subvector of a hybrid vector and then attempt to
	/// modify it (e.g. with push_back).
	typedef SparseSubvector<const Hybrid01VectorInterface, VectorRepresentationTypes::Hybrid01> SubvectorType;

	/// Const subvector-type: to be used when the parent is const
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	typedef SparseSubvector<const Hybrid01VectorInterface, VectorRepresentationTypes::Hybrid01> ConstSubvectorType;

	/// Aligned subvector-type: similar to subvector, but may
	/// restrict the start and possibly end indices to be integer
	/// multiples of the constant align and should then provide
	/// better performance. In the case of hybrid vectors, this
	/// usually means aligning the subvector along a
	/// word-boundary.
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	typedef SparseSubvector<Hybrid01VectorInterface, VectorRepresentationTypes::Hybrid01> AlignedSubvectorType;

	/// Const version of AlignedSubvectorType
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	typedef SparseSubvector<const Hybrid01VectorInterface, VectorRepresentationTypes::Hybrid01> ConstAlignedSubvectorType;

	/// Alignment-value: aligned subvectors must have start and
	/// possibly ending indices which are integer multiples of
	/// this value. For hybrid vectors this is normally the number
	/// of bits in a word.
	static const int align = WordTraits<word_type>::bits;

	//@}

	/// @name STL vector-interface
	//@{

	/// @name Iterator-access
	//@{

	/// Obtain an iterator for the beginning of the vector
	///
	/// Note: It must always be possible to derefence begin, even
	/// when the vector is empty. In that case it should yield the
	/// pair (0, 0).
	iterator begin ();

	/// Obtain a const iterator for the beginning of the vector
	///
	/// Note: It must always be possible to derefence begin, even
	/// when the vector is empty. In that case it should yield the
	/// pair (0, 0).
	const_iterator begin  () const;

	/// Obtain an iterator for one past the end of the vector
	///
	/// Important: Unlike the case with STL-vectors, is must
	/// always be possible to dereference end. The result should
	/// be a pair consisting of an index which is past the end of
	/// the vector (e.g. the last valid index plus one) and the
	/// zero word. This is required to avoid problems with the
	/// hybrid subvector.
	iterator end ();

	/// Obtain a const iterator for one past the end of the vector
	///
	/// Important: Unlike the case with STL-vectors, is must
	/// always be possible to dereference end. The result should
	/// be a pair consisting of an index which is past the end of
	/// the vector (e.g. the last valid index plus one) and the
	/// zero word. This is required to avoid problems with the
	/// hybrid subvector.
	const_iterator end () const;

	//@}

	/// @name Access to beginning and ending entries
	//@{

	/// Obtain a reference to the first entry
	inline reference front () { return *(begin ()); }

	/// Obtain a const reference to the first entry
	inline const_reference front () const { return *(begin ()); }

	/// Obtain a reference to the last entry
	inline reference back () { return *(end () - 1); }

	/// Obtain a const reference to the last entry
	inline const_reference back () const { return *(end () - 1); }

	//@}

	/// @name Manipulation of the vector
	///
	/// If the vector is constant then all of the methods of this
	/// section may be left out.
	//@{

	/** Append an entry to the vector
	 *
	 * This method is parametrised so that it may be used with
	 * STL-pairs as well as a derefenced iterator in case they are
	 * different (which is the case in, e.g., HybridVector). The
	 * type T should have the interface of an STL pair.
	 */
	template <class T>
	inline void push_back (T v) { insert (end (), v); }

	/// Clear all entries from the vector, yielding the zero-vector
	inline void clear () { erase (begin (), end ()); }

	/** Resize the vector a particular size
	 *
	 * This is not strictly required by the interface but may be
	 * useful in the implementation of assign below. It may leave
	 * the state of the newly allocated entries undefined.
	 */
	void resize (size_type s);

	/** Replace the vector with the contents of the iterators from first to last
	 *
	 * This should clear the vector and then copy the contents of
	 * the given iterators to the vector. It may be useful to
	 * specialise this on the iterator-tag so that it can proceed
	 * more quickly in the case of a forward-iterator.
	 */
	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
	{
		clear ();

		while (first != last)
			push_back (*first++);
	}

	/** Insert the given value into the given position
	 *
	 * The function is only required by SparseSubvector. If the
	 * user provides her own subvector-type which does not require
	 * this function, then it may be left out.
	 */
	template <class T>
	iterator insert (iterator pos, const T &x);

	/** Erase the entry at position pos
	 *
	 * This is in particular required by SparseMatrix::eraseEntry
	 * if this vector is used as a row-type.
	 */
	iterator erase (iterator pos);

	/** Erase all entries from first to last
	 *
	 * The function is only required by SparseSubvector. If the
	 * user provides her own subvector-type which does not require
	 * this function, then it may be left out.
	 */
	iterator erase (iterator first, iterator last);

	//@}

	/// @name Information on the vector
	//@{

	/// Return the number of entries of the vector
	inline size_type size () const { return end () - begin (); }

	/// Determine whether the vector is empty (i.e. the zero-vector); return true if so and false otherwise
	inline bool empty () const { return end () == begin (); }

	//@}

	//@}
};

} // namespace LELA

#endif // __LELA_VECTOR_HYBRID01_INTERFACE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
