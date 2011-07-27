/* lela/vector/hybrid.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * 0-1 vector in the hybrid format
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_HYBRID_H
#define __LELA_VECTOR_HYBRID_H

#include "lela/vector/sparse.h"

namespace LELA
{

struct HybridSubvectorWordAlignedTag;

/** Hybrid vector
 *
 * This is a canonical implementation of the hybrid 0-1
 * vector-representation-type.
 *
 * \ingroup vector
 */
template <class _Endianness, class IndexType, class WordType>
class HybridVector : public SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> >
{
public:
	typedef SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> > parent_type;

	typedef VectorRepresentationTypes::Hybrid01 RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef HybridVector ContainerType;
	typedef SparseSubvector<HybridVector, VectorRepresentationTypes::Hybrid01> ConstSubvectorType;
	typedef SparseSubvector<HybridVector, HybridSubvectorWordAlignedTag> AlignedSubvectorType;
	typedef SparseSubvector<const HybridVector, HybridSubvectorWordAlignedTag> ConstAlignedSubvectorType;
	static const int align = WordTraits<WordType>::bits;

	typedef _Endianness Endianness;
	typedef IndexType index_type;
	typedef WordType word_type;

	typedef typename parent_type::iterator iterator;
	typedef typename parent_type::const_iterator const_iterator;
	typedef typename parent_type::reverse_iterator reverse_iterator;
	typedef typename parent_type::const_reverse_iterator const_reverse_iterator;
	typedef typename parent_type::reference reference;
	typedef typename parent_type::const_reference const_reference;
	typedef typename parent_type::size_type size_type;

	HybridVector ()
		{ init_vectors (); }
		
	template <class IV, class EV>
	HybridVector (IV &iv, EV &ev)
		: SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> > (iv, ev)
	{
		parent_type::_idx.insert (parent_type::_idx.begin (), iv.front () - 1);
		parent_type::_elt.insert (parent_type::_elt.begin (), 0ULL);
		parent_type::_idx.insert (parent_type::_idx.end (), iv.back () + 1);
		parent_type::_elt.insert (parent_type::_elt.end (), 0ULL);
	}

	template <class IIt, class EIt>
	HybridVector (IIt idx_begin, IIt idx_end, EIt elt_begin)
		: SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> > (idx_begin, idx_end, elt_begin)
	{
		parent_type::_idx.insert (parent_type::_idx.begin (), *idx_begin - 1);
		parent_type::_elt.insert (parent_type::_elt.begin (), 0ULL);
		parent_type::_idx.insert (parent_type::_idx.end (), (idx_begin == idx_end) ? 0 : (*(idx_end - 1) + 1));
		parent_type::_elt.insert (parent_type::_elt.end (), 0ULL);
	}

	inline iterator               begin  ()       { return iterator (parent_type::_idx.begin () + 1, parent_type::_elt.begin () + 1); }
	inline const_iterator         begin  () const { return const_iterator (parent_type::_idx.begin () + 1, parent_type::_elt.begin () + 1); }
	inline iterator               end    ()       { return iterator (parent_type::_idx.end () - 1, parent_type::_elt.end () - 1); }
	inline const_iterator         end    () const { return const_iterator (parent_type::_idx.end () - 1, parent_type::_elt.end () - 1); }

	inline reverse_iterator       rbegin ()       { return reverse_iterator (end ()); }
	inline const_reverse_iterator rbegin () const { return const_reverse_iterator (end ()); }
	inline reverse_iterator       rend   ()       { return reverse_iterator (begin ()); }
	inline const_reverse_iterator rend   () const { return const_reverse_iterator (begin ()); }

	inline reference       operator[] (size_type n)       { return *(begin () + n); }
	inline const_reference operator[] (size_type n) const { return *(begin () + n); }

	inline reference       at (size_type n)
	{
		if (n >= size ())
			throw std::out_of_range ("n");
		else
			return (*this)[n];
	}

	inline const_reference at (size_type n) const
	{
		if (n >= size ())
			throw std::out_of_range ("n");
		else
			return (*this)[n];
	}

	inline reference       front     ()            { return *(begin ()); }
	inline const_reference front     () const      { return *(begin ()); }
	inline reference       back      ()            { return *(end () - 1); }
	inline const_reference back      () const      { return *(end () - 1); }

	template <class T>
	inline void            push_back (T v)
	{
		parent_type::_idx.back () = v.first;
		parent_type::_elt.back () = v.second;
		parent_type::_idx.push_back (v.first + 1);
		parent_type::_elt.push_back (0ULL);
	}

	inline void            clear     ()            { parent_type::_idx.clear (); parent_type::_elt.clear (); init_vectors (); }
	inline void            resize    (size_type s) { parent_type::_idx.resize (s + 2); parent_type::_elt.resize (s + 2); }

	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
	{
		clear ();

		while (first != last)
			push_back (*first++);

		fix ();
	}

	void assign (const_iterator first, const_iterator last)
		{ parent_type::assign (first - 1, last + 1); }

	inline size_type       size      () const      { return parent_type::_idx.size () - 2;  }
	inline bool            empty     () const      { return parent_type::_idx.size () == 2; }
	inline size_type       max_size  () const      { return std::min (parent_type::_idx.max_size (), parent_type::_elt.max_size ()) - 2;  }

	// Call this after building a hybrid-vector, or at least after pushing the first entry onto it
	inline void fix ()
		{ parent_type::_idx[0] = parent_type::_idx[1] - 1; }

private:
	void init_vectors ()
	{
		parent_type::_idx.push_back ((index_type) -1);
		parent_type::_idx.push_back (0);
		parent_type::_elt.push_back (0ULL);
		parent_type::_elt.push_back (0ULL);
	}
};

} // namespace LELA

namespace std
{

// Specialisation of std::swap to sparse vectors
template <class Endianness, class IndexType, class WordType>
void swap (LELA::HybridVector<Endianness, IndexType, WordType> &v1, LELA::HybridVector<Endianness, IndexType, WordType> &v2)
	{ v1.swap (v2); }

} // namespace std

#endif // __LELA_VECTOR_HYBRID_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
