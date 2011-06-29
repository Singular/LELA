/* linbox/vector/traits.h
 * Copyright 1999-2001 William J Turner,
 *           2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license details
 */

#ifndef __LINBOX_vector_traits_H
#define __LINBOX_vector_traits_H

#include <vector>	// STL vectors
#include <list>		// STL lists
#include <deque>	// STL deques
#include <utility>      // STL pairs
#include <functional>   // STL functions
#include <map>          // STL maps
#include <algorithm>    // STL algorithms

#include "linbox/ring/archetype.h"
#include "linbox/vector/bit-iterator.h"

namespace LinBox
{

/** @name Vector traits.
 * Vector traits are use to allow template specialization to choose different
 * code for dense and sparse vectors.
	\ingroup vector
 */
//@{

/** \brief Vector-representation-types.
 *
 * Over general rings there are two of vector-representations: dense
 * and sparse.
 *
 * A dense vector must support the interface of std::vector with the
 * exception of operations which invalidate iterators (e.g. insert,
 * erase, push_back, and so on). A dense vector exists in the free
 * module over its ring defined by its size. Thus it must always be
 * allocated to the correct size before any operation is performed on
 * it.
 *
 * A sparse vector must support an interface equivalent to
 * std::vector<std::pair<index_type, element_type> >, a vector of
 * index-entry pairs. Here index_type is an integral type giving the
 * index of a given entry and element_type is the type of the
 * underlying field-element. It must include appropriate typedefs for
 * index_type and element_type. Indices must be in strictly ascending
 * order and zero is not allowed as the element in a pair (such an
 * entry should be removed from the vector). The class @ref
 * SparseVector implements this interface.
 *
 * Note that the module in which a sparse vector exists is a priori
 * undefined, since arbitrary entries may be added.
 *
 * Over GF2 there are three representation-types: dense, sparse, and
 * hybrid.
 *
 * A dense 0-1 vector is similar to that over general rings with
 * the addition of a type word_iterator which iterates over entire
 * words (of type word_type -- usually uint64 -- which must be a
 * typedef in the vector) in the vector so that operations can be done
 * in parallel. It must correspondingly support the methods word_begin
 * and word_end. The class @ref BitVector implements this interface.
 *
 * A sparse 0-1 vector must support the interface of std::vector,
 * treated as a vector of indices, which must be sorted in strictly
 * ascending order. It must support the full interface, including
 * operations which invalidate iterators.
 *
 * A hybrid 0-1 vector must support the interface of
 * std::vector<std::pair<index_type, word_type> >, a vector of
 * index-word pairs. The vector-entries are divided into blocks, each
 * one the size of a word (of type word_type). The index is then the
 * index of the block and the word is the block itself. The word in a
 * pair must not be identically zero (such a pair should be removed
 * from the vector). Indices must be in strictly ascending order.
 *
 * In order that subvectors of hybrid vectors work properly, it must
 * always be possible to dereference an initialised iterator on a
 * hybrid vector, even if the iterator is initialised with end () and
 * even if the vector is empty. Dereferencing end () should result in
 * a pair consisting of the previous index plus one and the zero
 * word. This can be achieved by placing a sentry at the end of the
 * vector.
 *
 * Dense and hybrid 0-1 vectors must also define a type Endianness,
 * which indicates the order of entries in a word. The classes @ref
 * BigEndian<word_type> and @ref LittleEndian<word_type> define the
 * conventions respectively that the vector e_i corresponds to the
 * word with value 2^(N-i), where N is the number of bits in a word,
 * and that the vector e_i corresponds to the word with value 2^i.
 */
namespace VectorRepresentationTypes
{
	struct Generic {};
	struct Dense : public Generic {};
	struct Sparse : public Generic {};
	struct Dense01 : public Generic {};
	struct Sparse01 : public Generic {};
	struct Hybrid01 : public Generic {};
};

/** Vector storage-types
 *
 * The tags in this namespace indicate the underlying storage of a
 * vector. They are useful when interfacing with lower-level
 * libraries.
 *
 * There are three storage-types: generic, transformed, and real.
 *
 * Generic indicates that no assumptions at all should be made about
 * the underlying storage of the vector. This is useful if the code
 * using this library has created a virtual vector which does not
 * correspond directly to one in memory.
 *
 * Transformed (meaningful for sparse and hybrid 0-1 vectors only)
 * indicates that a vector has undergone a transformation from a pair
 * of vectors to a vector of pairs.
 *
 * Real means that the vector represents a true array in memory. A
 * pointer to the beginning of the vector may then be attained by
 * dereferencing and then taking the address of the iterator begin ().
 */
namespace VectorStorageTypes
{
	struct Generic {};
	struct Transformed : public Generic {};
	struct Real : public Generic {};
}

/** Vector traits template structure.
 *
 * This structure defines information about a vector on which methods
 * may specialise and from which they may generically derive
 * information.
 */
template <class Vector> struct DefaultVectorTraits
{
	/** Representation-type
	 *
	 * See @ref VectorRepresentationTypes for definitions. On this
	 * typedef methods which vary on the representation-type
	 * should specialise.
	 */
	typedef typename Vector::RepresentationType RepresentationType;

	/** Storage-type
	 *
	 * See @ref VectorStorageTypes for definitions.
	 */
	typedef typename Vector::StorageType StorageType;

	/** Container-type
	 *
	 * This defines a vector which can be declared in generic code
	 * with the same representation-type as Vector.
	 */
	typedef typename Vector::ContainerType ContainerType;

	/** Subvector
	 *
	 * This defines a subvector of the given vector which can be
	 * constructed in generic code. It must support a constructor
	 * SubvectorType (v, start_index, size)
	 */
	typedef typename Vector::SubvectorType SubvectorType;

	/** Const subvector
	 *
	 * Const version of @ref Subvector
	 */
	typedef typename Vector::ConstSubvectorType ConstSubvectorType;

	/** Aligned subvector
	 *
	 * This type may perform more quickly than @ref Subvector but
	 * the choice of start_index is restricted to an
	 * integer-multiple of @ref align. If align is one, then it is
	 * equivalent to Subvector.
	 *
	 * This is useful for recursive algorithms with some control
	 * over where a vector is to be split.
	 */
	typedef typename Vector::AlignedSubvectorType AlignedSubvectorType;

	/** Const aligned subvector
	 *
	 * Const version of @ref AlignedSubvector
	 */
	typedef typename Vector::ConstAlignedSubvectorType ConstAlignedSubvectorType;

	/** Alignment factor
	 *
	 * The starting-index and possibly size of a
	 * AlignedSubvectorType must be an integral multiple of this
	 * value.
	 */
	static const int align = 1;
};

/** Version of the above for vectors over GF2
 *
 * Identical to @ref DefaultVectorTraits, but used for vectors over
 * GF2.
 *
 * This is required because sparse vectors over GF2 and dense vectors
 * over a field with element-type unsigned int have the same type
 */
template <class Vector>
struct GF2VectorTraits
{
	typedef typename Vector::RepresentationType RepresentationType;
	typedef typename Vector::StorageType StorageType;
	typedef typename Vector::ContainerType ContainerType;
	typedef typename Vector::SubvectorType SubvectorType;
	typedef typename Vector::ConstSubvectorType ConstSubvectorType;
	typedef typename Vector::AlignedSubvectorType AlignedSubvectorType;
	typedef typename Vector::ConstAlignedSubvectorType ConstAlignedSubvectorType;
	static const int align = Vector::align;
};

/** Version which takes the field-element as a hint as to which
 * version of the above to use (useful when the field is not available)
 */
template <class Element, class Vector>
struct ElementVectorTraits
{
	typedef typename DefaultVectorTraits<Vector>::RepresentationType RepresentationType;
	typedef typename DefaultVectorTraits<Vector>::ContainerType ContainerType;
	typedef typename DefaultVectorTraits<Vector>::SubvectorType SubvectorType;
	typedef typename DefaultVectorTraits<Vector>::ConstSubvectorType ConstSubvectorType;
	typedef typename DefaultVectorTraits<Vector>::AlignedSubvectorType AlignedSubvectorType;
	typedef typename DefaultVectorTraits<Vector>::ConstAlignedSubvectorType ConstAlignedSubvectorType;
	static const int align = DefaultVectorTraits<Vector>::align;
};

/** Version using the field rather than its element-type
 */
template <class Ring, class Vector>
struct VectorTraits
{
	typedef typename ElementVectorTraits<typename Ring::Element, Vector>::RepresentationType RepresentationType;
	typedef typename ElementVectorTraits<typename Ring::Element, Vector>::ContainerType ContainerType;
	typedef typename ElementVectorTraits<typename Ring::Element, Vector>::SubvectorType SubvectorType;
	typedef typename ElementVectorTraits<typename Ring::Element, Vector>::ConstSubvectorType ConstSubvectorType;
	typedef typename ElementVectorTraits<typename Ring::Element, Vector>::AlignedSubvectorType AlignedSubvectorType;
	typedef typename ElementVectorTraits<typename Ring::Element, Vector>::ConstAlignedSubvectorType ConstAlignedSubvectorType;
	static const int align = ElementVectorTraits<typename Ring::Element, Vector>::align;
};

/** Specialisation of ElementVectorTraits for vectors over GF2 */
template <class Vector>
struct ElementVectorTraits<bool, Vector>
{
	typedef typename GF2VectorTraits<Vector>::RepresentationType RepresentationType;
	typedef typename GF2VectorTraits<Vector>::ContainerType ContainerType;
	typedef typename GF2VectorTraits<Vector>::SubvectorType SubvectorType;
	typedef typename GF2VectorTraits<Vector>::ConstSubvectorType ConstSubvectorType;
	typedef typename GF2VectorTraits<Vector>::AlignedSubvectorType AlignedSubvectorType;
	typedef typename GF2VectorTraits<Vector>::ConstAlignedSubvectorType ConstAlignedSubvectorType;
	static const int align = GF2VectorTraits<Vector>::align;
};

/** Utility-functions for vectors */

class VectorUtils
{
	template <class Element, class Vector>
	static inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorRepresentationTypes::Dense)
		{ a = v[i]; return true; }

	template <class Element, class Vector>
	static inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorRepresentationTypes::Sparse)
	{
		typename Vector::const_iterator j;

		if (v.size () == 0)
			return false;

		j = std::lower_bound (v.begin (), v.end (), i, CompareSparseEntries ());

		if (j == v.end () || j->first != i)
			return false;
		else {
			a = j->second;
			return true;
		}
	}

	template <class Element, class Vector>
	static inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorRepresentationTypes::Dense01)
		{ a = v[i]; return true; }

	template <class Element, class Vector>
	static inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorRepresentationTypes::Sparse01)
	{
		typename Vector::const_iterator j;

		if (v.size () == 0)
			return false;

		j = std::lower_bound (v.begin (), v.end (), i);

		if (j == v.end () || *j != i)
			return false;
		else {
			a = true;
			return true;
		}
	}

	template <class Element, class Vector>
	static inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorRepresentationTypes::Hybrid01)
	{
		typename Vector::const_iterator idx;

		idx = std::lower_bound (v.begin (), v.end (), i >> WordTraits<typename Vector::word_type>::logof_size, CompareSparseEntries ());

		if (idx != v.end () && idx->first == i >> WordTraits<typename Vector::word_type>::logof_size) {
			a = idx->second & Vector::Endianness::e_j (i & WordTraits<typename Vector::word_type>::pos_mask);
			return true;
		} else
			return false;
	}

	template <class Vector>
	static inline void ensureDimSpecialized (Vector &v, size_t n, VectorRepresentationTypes::Dense)
		{ v.resize (n); }

	template <class Vector>
	static inline void ensureDimSpecialized (Vector &v, size_t n, VectorRepresentationTypes::Sparse)
		{}

	template <class Vector>
	static inline void ensureDimSpecialized (Vector &v, size_t n, VectorRepresentationTypes::Dense01)
		{ v.resize (n); }

	template <class Vector>
	static inline void ensureDimSpecialized (Vector &v, size_t n, VectorRepresentationTypes::Sparse01)
		{}

	template <class Vector>
	static inline void ensureDimSpecialized (Vector &v, size_t n, VectorRepresentationTypes::Hybrid01)
		{}

	template <class Vector>
	static inline bool hasDimSpecialized (const Vector &v, size_t n, VectorRepresentationTypes::Dense)
		{ return v.size () == n; }

	template <class Vector>
	static inline bool hasDimSpecialized (const Vector &v, size_t n, VectorRepresentationTypes::Sparse)
		{ return v.empty () || v.back ().first < n; }

	template <class Vector>
	static inline bool hasDimSpecialized (const Vector &v, size_t n, VectorRepresentationTypes::Dense01)
		{ return v.size () == n; }

	template <class Vector>
	static inline bool hasDimSpecialized (const Vector &v, size_t n, VectorRepresentationTypes::Sparse01)
		{ return v.empty () || v.back () < n; }

	template <class Vector>
	static inline bool hasDimSpecialized (const Vector &v, size_t n, VectorRepresentationTypes::Hybrid01)
	{
		if (v.empty ())
			return true;
		else if ((v.back ().first << WordTraits<typename Vector::word_type>::logof_size) >= (long) n)
			return false;
		else if (v.back ().first == (n >> WordTraits<typename Vector::word_type>::logof_size))
			return (v.back ().second & Vector::Endianness::mask_right (n & WordTraits<typename Vector::word_type>::pos_mask)) == 0;
		else
			return true;
	}

	template <class Vector>
	static inline bool isValidSpecialized (const Vector &v, VectorRepresentationTypes::Dense)
		{ return true; }

	template <class Vector>
	static inline bool isValidSpecialized (const Vector &v, VectorRepresentationTypes::Sparse)
	{
		if (v.empty ())
			return true;

		typename Vector::const_iterator i = v.begin (), i_next = v.begin ();

		for (++i_next; i_next != v.end (); ++i_next, ++i)
			if (i->first >= i_next->first)
				return false;

		return true;
	}

	template <class Vector>
	static inline bool isValidSpecialized (const Vector &v, VectorRepresentationTypes::Dense01)
		{ return true; }

	template <class Vector>
	static inline bool isValidSpecialized (const Vector &v, VectorRepresentationTypes::Sparse01)
	{
		if (v.empty ())
			return true;

		typename Vector::const_iterator i = v.begin (), i_next = v.begin ();

		for (++i_next; i_next != v.end (); ++i_next, ++i)
			if (*i >= *i_next)
				return false;

		return true;
	}

	template <class Vector>
	static inline bool isValidSpecialized (const Vector &v, VectorRepresentationTypes::Hybrid01)
	{
		if (v.empty ())
			return true;

		typename Vector::const_iterator i = v.begin (), i_next = v.begin ();

		for (++i_next; i_next != v.end (); ++i_next, ++i)
			if (i->first >= i_next->first || i->second == 0)
				return false;

		return true;
	}

public:
	/** Closure to compare entries in a sparse vector
	 *
	 * This can be used either to sort entries by index
	 * (e.g. after performing a permutation) or to perform a
	 * binary search on entries.
	 */
	class CompareSparseEntries
	{
	public:
		template<typename PairType>
		inline bool operator () (const PairType &i, uint16 j) const
			{ return i.first < j; }

		template<typename PairType>
		inline bool operator () (const PairType &i, uint32 j) const
			{ return i.first < j; }

		template<typename PairType>
		inline bool operator () (const PairType &i, uint64 j) const
			{ return i.first < j; }

		template<typename PairType>
		inline bool operator () (size_t i, const PairType &j) const
			{ return i < j.first; }

		template<typename PairType1, typename PairType2>
		inline bool operator () (const PairType1 &i, const PairType2 &j) const
			{ return i.first < j.first; }
	};

	/** Retrieve the entry at index i from the vector and store it in the element a
	 *
	 * If there is no entry at index i, then the element a is left unchanged
	 *
	 * @returns true if the entry at index i exists in the vector, false otherwise
	 */
	template <class Element, class Vector>
	static inline bool getEntry (const Vector &v, Element &a, size_t i) 
		{ return getEntrySpecialised<Element, Vector> (v, a, i, typename ElementVectorTraits<Element, Vector>::RepresentationType()); }

	/** Ensure that the vector v is defined over the free module of rank n */
	template <class Ring, class Vector>
	static inline void ensureDim (Vector &v, size_t n) 
		{ ensureDimSpecialized (v, n, typename VectorTraits<Ring, Vector>::RepresentationType()); }

	/// Determines whether the vector v can represent a vector of dimension n.
	/// @returns true if v can represent a vector of dimension n and false otherwise
	template <class Ring, class Vector>
	static inline bool hasDim (const Vector &v, size_t n) 
		{ return hasDimSpecialized (v, n, typename VectorTraits<Ring, Vector>::RepresentationType()); }

	/// Determines whether v is a valid vector of its format.
	/// @returns true if v is valid or false if there is an error
	template <class Ring, class Vector>
	static inline bool isValid (const Vector &v) 
		{ return isValidSpecialized (v, typename VectorTraits<Ring, Vector>::RepresentationType()); }

	/// Compute the image of an index under a permutation
	template <class Iterator>
	static typename Iterator::value_type::first_type permutationImage (typename Iterator::value_type::first_type x, Iterator P_start, Iterator P_end)
	{
		for (Iterator i = P_start; i != P_end; ++i) {
			if (i->first == x)
				x = i->second;
			else if (i->second == x)
				x = i->first;
		}

		return x;
	}
}; // class VectorUtils

// Forward declarations of types we're about to use
template <typename Iterator, typename ConstIterator = Iterator> class Subvector;
template <class Element, class IndexVector = std::vector<uint32>, class ElementVector = std::vector<Element> > class SparseVector;
template <class Vector, class Trait = typename DefaultVectorTraits<Vector>::RepresentationType> class SparseSubvector;

template <class Element>
struct DefaultVectorTraits< std::vector<Element> >
{ 
	typedef VectorRepresentationTypes::Dense RepresentationType;
	typedef VectorStorageTypes::Real StorageType;
	typedef std::vector<Element> ContainerType;
	typedef Subvector<typename std::vector<Element>::iterator, typename std::vector<Element>::const_iterator> SubvectorType;
	typedef Subvector<typename std::vector<Element>::const_iterator, typename std::vector<Element>::const_iterator> ConstSubvectorType;
	typedef Subvector<typename std::vector<Element>::iterator, typename std::vector<Element>::const_iterator> AlignedSubvectorType;
	typedef Subvector<typename std::vector<Element>::const_iterator, typename std::vector<Element>::const_iterator> ConstAlignedSubvectorType;
	static const int align = 1;
};

template <class Element>
struct DefaultVectorTraits< const std::vector<Element> >
{ 
	typedef VectorRepresentationTypes::Dense RepresentationType;
	typedef VectorStorageTypes::Real StorageType;
	typedef std::vector<Element> ContainerType;
	typedef Subvector<typename std::vector<Element>::iterator, typename std::vector<Element>::const_iterator> SubvectorType;
	typedef Subvector<typename std::vector<Element>::const_iterator, typename std::vector<Element>::const_iterator> ConstSubvectorType;
	typedef Subvector<typename std::vector<Element>::iterator, typename std::vector<Element>::const_iterator> AlignedSubvectorType;
	typedef Subvector<typename std::vector<Element>::const_iterator, typename std::vector<Element>::const_iterator> ConstAlignedSubvectorType;
	static const int align = 1;
};

template <class index_type> 
struct GF2VectorTraits<std::vector<index_type> >
{ 
	typedef VectorRepresentationTypes::Sparse01 RepresentationType;
	typedef VectorStorageTypes::Real StorageType;
	typedef std::vector<index_type> ContainerType;
	typedef SparseSubvector<ContainerType> SubvectorType;
	typedef SparseSubvector<const ContainerType> ConstSubvectorType;
	typedef SparseSubvector<ContainerType> AlignedSubvectorType;
	typedef SparseSubvector<const ContainerType> ConstAlignedSubvectorType;
	static const int align = 1;
};

template <class index_type> 
struct GF2VectorTraits<const std::vector<index_type> >
{ 
	typedef VectorRepresentationTypes::Sparse01 RepresentationType;
	typedef VectorStorageTypes::Real StorageType;
	typedef const std::vector<index_type> ContainerType;
	typedef SparseSubvector<ContainerType> SubvectorType;
	typedef SparseSubvector<ContainerType> ConstSubvectorType;
	typedef SparseSubvector<ContainerType> AlignedSubvectorType;
	typedef SparseSubvector<ContainerType> ConstAlignedSubvectorType;
	static const int align = 1;
};

// Forward-declaration

/** Canonical vector types
 *
 * This class includes some typedefs that avoid the necessity to typedef
 * the vector type whenever it is used. In a typical case, one would say
 * Vector<Field>::Dense for a dense vector and Vector<Field>::Sparse for
 * a sparse vector.
 */

template <class Element>
struct RawVector 
{
	typedef std::vector<Element> Dense;
	typedef SparseVector<Element> Sparse;
};

template <class Ring>
struct Vector : public RawVector<typename Ring::Element>
{
	typedef typename Ring::Element Element;
};

//@} Vector traits

} // namespace LinBox

#endif // __LINBOX_vector_traits_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
