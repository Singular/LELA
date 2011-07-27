/* lela/vector/dense-interface.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Prototype for dense vectors. This file is provided for documentation only.
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_DENSE_INTERFACE_H
#define __LELA_VECTOR_DENSE_INTERFACE_H

#include <iterator>

#include "linbox/vector/traits.h"
#include "linbox/vector/subvector.h"

namespace LELA
{

/** Dense vector-interface
 *
 * This class is a prototype of a vector meeting the dense interface
 * in LELA. It can be used as a starting-point for those wishing to
 * implement their own dense vector-types. It cannot be used on its
 * own.
 *
 * When feasible, suggested implementations of methods are provided.
 *
 * \ingroup vector
 */
class DenseVectorInterface
{
public:
	/// @name STL-defined types
	//@{

	/** Normal iterator
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
	class iterator;

	/** Const iterator
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
	 */
	class const_iterator;

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
	typedef VectorRepresentationTypes::Dense RepresentationType;

	/// Storage-type: see VectorStorageTypes
	/// Use Real if the vector represents a live array in memory, Generic if not
	typedef VectorStorageTypes::Generic StorageType;

	/// Container-type: this should be the type a user would use
	/// in generic code to construct another vector of the same
	/// type, e.g., if the user wishes to construct a copy.
	typedef std::vector<typename std::iterator_traits<iterator>::value_type> ContainerType;

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
	typedef Subvector SubvectorType;

	/// Const subvector-type: to be used when the parent is const
	///
	/// Normally it suffices to use Subvector for this type
	/// as given here.
	typedef Subvector<ConstIterator, ConstIterator> ConstSubvectorType;

	/// Aligned subvector-type: similar to subvector, but may
	/// restrict the start and possibly end indices to be integer
	/// multiples of the constant align and should then provide
	/// better performance. In nearly all cases, it should be
	/// identical to Subvector.
	///
	/// Normally it suffices to use Subvector for this type
	/// as given here.
	typedef Subvector AlignedSubvectorType;

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

	/// Return the number of entries of the vector
	inline size_type size () const { return end () - begin (); }

	/// Determine whether the vector is empty; return true if so and false otherwise
	inline bool empty () const { return end () == begin (); }

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
