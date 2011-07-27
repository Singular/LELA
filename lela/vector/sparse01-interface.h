/* lela/vector/sparse01-interface.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Prototype for sparse 0-1 vectors. This file is provided for documentation only.
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_SPARSE01_INTERFACE_H
#define __LELA_VECTOR_SPARSE01_INTERFACE_H

#include <iterator>

#include "linbox/vector/traits.h"
#include "linbox/vector/sparse-subvector.h"

namespace LELA
{

/** Sparse 0-1 vector interface
 *
 * This class is a prototype of a vector meeting the sparse 0-1
 * interface in LELA. It can be used as a starting-point for those
 * wishing to implement their own sparse 0-1 vector-types. It cannot
 * be used on its own.
 *
 * When feasible, suggested implementations of methods are provided.
 *
 * \ingroup vector
 */
class Sparse01VectorInterface
{
public:
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
	/// This should in general be an integral index-type
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

	/// Representation-type: always sparse 0-1
	typedef VectorRepresentationTypes::Sparse01 RepresentationType;

	/** Storage-type: see VectorStorageTypes
	 *
	 * Use Real if the vector represents a live array of indices
	 * in memory, Generic if not.
	 */
	typedef VectorStorageTypes::Generic StorageType;

	/// Container-type: this should be the type a user would use
	/// in generic code to construct another vector of the same
	/// type, e.g., if the user wishes to construct a copy.
	typedef std::vector<typename std::iterator_traits<Iterator>::value_type> ContainerType;

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
	typedef SparseSubvector<Sparse01VectorInterface, VectorRepresentationTypes::Sparse01> SubvectorType;

	/// Const subvector-type: to be used when the parent is const
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	typedef SparseSubvector<const Sparse01VectorInterface, VectorRepresentationTypes::Sparse01> ConstSubvectorType;

	/// Aligned subvector-type: similar to subvector, but may
	/// restrict the start and possibly end indices to be integer
	/// multiples of the constant align and should then provide
	/// better performance. In nearly all cases, it should be
	/// identical to Subvector.
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	typedef SparseSubvector<Sparse01VectorInterface, VectorRepresentationTypes::Sparse01> AlignedSubvectorType;

	/// Const version of AlignedSubvectorType
	///
	/// Normally it suffices to use SparseSubvector for this type
	/// as given here.
	typedef SparseSubvector<const Sparse01VectorInterface, VectorRepresentationTypes::Sparse01> ConstAlignedSubvectorType;

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

	/// Append an entry to the vector
	inline void push_back (value_type v) { insert (end (), v); }

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
	iterator insert (iterator pos, value_type x);

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

#endif // __LELA_VECTOR_SPARSE01_INTERFACE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
