/* linbox/vector/vector-traits.h
 * Copyright 1999-2001 William J Turner,
 *           2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Implemented the Rootbeer meeting changes. Made VectorTag parametrized
 * and added typedef of those parameters to Traits. So now VectorCategory
 * for each vector has a "reference" back to the VectorTrait of each specific
 * vector (list of pairs, deque of pairs, etc.) through a typedef Trait. This
 * allows for generic manipulation of all vectors and placing the 
 * vector-implementation dependent code into VectorTraits only - as is done now
 * with the function sort.
 *
 * see COPYING for license details
 *
 * ------------------------------------ 
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

#include "linbox/field/archetype.h"
#include "linbox/field/rebind.h"
#include "linbox/vector/bit-iterator.h"

namespace LinBox
{

/** @name Vector traits.
 * Vector traits are use to allow template specialization to choose different
 * code for dense and sparse vectors.
	\ingroup vector
 */
//@{

/** \brief List of vector categories.
 *
 * This structure contains three structures: one relating to dense vectors,
 * one relating to sparse vectors implemented as sequences of pairs, and 
 * one relating to sparse vectors implemented as associative containers.
 * These types allow us to use template specialization to use different 
 * code for different types of vectors.
 */
struct VectorCategories
{
	struct GenericVectorTag
	{ 
                friend std::ostream &operator << (std::ostream &o, 
                                                 const GenericVectorTag &)
			{ return o << "GenericVectorTag"; } 
	};
            
	// These are valid for GF2 only

	struct SparseZeroOneVectorTag : public GenericVectorTag
	{ 
                friend std::ostream &operator << (std::ostream &o, 
						  const SparseZeroOneVectorTag &)
			{ return o << "SparseZeroOneVectorTag"; } 
	};

	// These are valid for all fields

	struct DenseVectorTag : public SparseZeroOneVectorTag
	{ 
		// Inherits from SparseZeroOneVectorTag:
		// This simplifies gf2 management of vectors
		// and makes Dense vector also a generic vector
                friend std::ostream &operator << (std::ostream &o, 
						  const DenseVectorTag &)
			{ return o << "DenseVectorTag"; } 
	};
            
	struct SparseVectorTag : public GenericVectorTag
	{ 
                friend std::ostream &operator << (std::ostream &o, 
						  const SparseVectorTag &)
			{ return o << "SparseVectorTag"; } 
	};
            
	struct DenseZeroOneVectorTag : public DenseVectorTag
	{
                friend std::ostream &operator << (std::ostream &o, 
						  const DenseZeroOneVectorTag &)
			{ return o << "DenseZeroOneVectorTag"; } 
	};

	struct HybridZeroOneVectorTag : public GenericVectorTag
	{
                friend std::ostream &operator << (std::ostream &o, 
						  const HybridZeroOneVectorTag &)
			{ return o << "HybridZeroOneVectorTag"; } 
	};

	struct HybridZeroOneSequenceVectorTag : public GenericVectorTag
	{ 
                friend std::ostream &operator << (std::ostream &o, 
						  const HybridZeroOneSequenceVectorTag &)
			{ return o << "HybridZeroOneSequenceVectorTag"; } 
	};
};

/** Vector traits template structure.
 * By default, it tries to take all information from the vector class,
 * but it cannot usually do this.  For example, the vector_category
 * type is not defined in STL types, so this must be done through 
 * template specialization.
 * @param Vector \Ref{LinBox} dense or sparse vector.
 */
template <class Vector> struct VectorTraits
{
	typedef typename Vector::VectorCategory VectorCategory;
	typedef Vector VectorType;

	// These are defined for all STL vectors and sequence containers.
};

// Namespace containing some useful generic functions

namespace VectorWrapper 
{
	class CompareSparseEntries
	{
	public:
		template<typename PairType>
		inline bool operator () (const PairType &i, const size_t j) const
			{ return i.first < j; }

		template<typename PairType1, typename PairType2>
		inline bool operator () (const PairType1 &i, const PairType2 &j) const
			{ return i.first < j.first; }
	};

	template <class Element, class Vector>
	inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorCategories::DenseVectorTag)
		{ a = v[i]; return true; }

	template <class Element, class Vector>
	inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorCategories::SparseVectorTag)
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
	inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorCategories::SparseZeroOneVectorTag)
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
	inline bool getEntrySpecialised (const Vector &v, Element &a, size_t i, VectorCategories::HybridZeroOneVectorTag)
	{
		typedef typename Vector::second_type::Endianness Endianness;
		typedef typename std::iterator_traits<typename Vector::second_type::const_word_iterator>::value_type word_type;

		typename Vector::first_type::const_iterator idx;

		idx = std::lower_bound (v.first.begin (), v.first.end (), i >> WordTraits<word_type>::logof_size);

		if (idx != v.first.end () && *idx == i >> WordTraits<word_type>::logof_size) {
			a = *(v.second.wordBegin () + (idx - v.first.begin ())) & Endianness::e_j (i & WordTraits<word_type>::pos_mask);
			return a;
		} else
			return false;
	}

	template <class Element, class Vector>
	inline bool getEntry (const Vector &v, Element &a, size_t i) 
		{ return getEntrySpecialised<Element, Vector> (v, a, i, typename VectorTraits<Vector>::VectorCategory()); }

	template <class Vector>
	inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::DenseVectorTag)
		{ v.resize (n); }

	template <class Vector>
	inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseVectorTag)
		{}

	template <class Vector>
	inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::DenseZeroOneVectorTag)
		{ v.resize (n); }

	template <class Vector>
	inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseZeroOneVectorTag)
		{}

	template <class Vector>
	inline void ensureDim (Vector &v, size_t n) 
		{ ensureDimSpecialized (v, n, typename VectorTraits<Vector>::VectorCategory()); }
} // Namespace VectorWrapper

// Specialization for STL vectors
template <class Element>
struct VectorTraits< std::vector<Element> >
{ 
	typedef std::vector<Element> VectorType;
	typedef typename VectorCategories::DenseVectorTag VectorCategory; 
};

// Specialization for STL vectors of pairs of size_t and elements
template <class Element> 
struct VectorTraits< std::vector< std::pair<uint32, Element> > >
{ 
	typedef std::vector< std::pair<uint32, Element> > VectorType;
	typedef typename VectorCategories::SparseVectorTag VectorCategory; 

	static void sort (VectorType& v) { std::stable_sort (v.begin (), v.end (), VectorWrapper::CompareSparseEntries ()); }
};

// Specialization for STL lists of pairs of size_t and elements
template <class Element> 
struct VectorTraits< std::list< std::pair<size_t, Element> > >
{ 
	typedef std::list< std::pair<size_t, Element> > VectorType;
	typedef typename VectorCategories::SparseVectorTag VectorCategory; 

	static void sort (VectorType& v) { v.sort (VectorWrapper::CompareSparseEntries ()); }
};

// Specialization for STL singly linked lists of pairs of size_t and elements
template <class Element> 
struct VectorTraits< std::deque< std::pair<size_t, Element> > >
{ 
	typedef std::deque< std::pair<size_t, Element> > VectorType;
	typedef typename VectorCategories::SparseVectorTag VectorCategory; 

	static void sort (VectorType& v) { std::stable_sort (v.begin (), v.end (), VectorWrapper::CompareSparseEntries ()); }
};
  
// Specialization for a const STL vector of size_t's
template <> 
struct VectorTraits< const std::vector<size_t> >
{ 
	typedef std::vector<size_t> VectorType;
	typedef VectorCategories::SparseZeroOneVectorTag VectorCategory; 
};

// Now we create some "canonical" vector types, so that users don't 
// always have to typedef everything

// Forward-declaration
template <class Element, class IndexVector = std::vector<uint32>, class ElementVector = std::vector<Element> >
class SparseVector;

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

template <class Field>
struct Vector : public RawVector<typename Field::Element>
{
	typedef typename Field::Element Element;
        
	template<class U>
	struct rebind
	{
                typedef Vector<U> other;
	};
};

template<class T, class U>
struct Rebind< std::vector<T>, U >
{
	typedef typename Vector<U>::Dense other;
};

template<class T, class U>
struct Rebind< std::pair<std::vector<size_t>, std::vector<T> >,U >
{
	typedef typename Vector<U>::Sparse other;
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
