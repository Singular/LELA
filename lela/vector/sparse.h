/* lela/vector/sparse.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Sparse vector which makes a pair of vectors appear to be a vector of pairs
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_SPARSE_H
#define __LELA_VECTOR_SPARSE_H

#include <vector>
#include <algorithm>
#include <cstddef>

#include "lela/util/property.h"
#include "lela/vector/traits.h"

namespace LELA
{

/// Reference to an (index, entry)-pair
template <class IndexIterator, class ElementIterator>
class SparseVectorReference {
public:
	typedef typename std::iterator_traits<IndexIterator>::value_type first_type;
	typedef typename std::iterator_traits<ElementIterator>::value_type second_type;

	Property<IndexIterator> first;
	Property<ElementIterator> second;

	SparseVectorReference () {}
	SparseVectorReference (IndexIterator idx, ElementIterator elt) : first (idx), second (elt) {}

	template <class IIt2, class EIt2>
	SparseVectorReference (const SparseVectorReference<IIt2, EIt2> &ref) : first (ref.first), second (ref.second) {}

	SparseVectorReference &operator = (const std::pair<first_type, second_type> &r)
		{ *first._i = r.first; *second._i = r.second; return *this; }

	SparseVectorReference &operator = (const SparseVectorReference &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }

	template <class II, class EI>
	SparseVectorReference &operator = (const SparseVectorReference<II, EI> &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }

	operator std::pair<first_type, second_type> () const
		{ return std::pair<first_type, second_type> (*first._i, *second._i); }

	template <class II, class EI>
	bool operator == (const SparseVectorReference<II, EI> &r) const
		{ return (first == r.first) && (second == r.second); }

	template <class II, class EI>
	bool operator != (const SparseVectorReference<II, EI> &r) const
		{ return (first != r.first) || (second != r.second); }
};

/// Iterator for sparse vectors
template <class IndexIterator, class ElementIterator, class ConstIndexIterator = IndexIterator, class ConstElementIterator = ElementIterator>
class SparseVectorIterator
{
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef SparseVectorReference<IndexIterator, ElementIterator> reference;
	typedef std::pair<typename IndexIterator::value_type, typename ElementIterator::value_type> value_type;
	typedef const SparseVectorReference<ConstIndexIterator, ConstElementIterator> const_reference;
	typedef reference *pointer;
	typedef const_reference *const_pointer;
	typedef ptrdiff_t difference_type;
	typedef size_t size_type;

	SparseVectorIterator () {}
	SparseVectorIterator (IndexIterator idx, ElementIterator elt) : _ref (idx, elt) {}
	SparseVectorIterator (const SparseVectorIterator &i) : _ref (i._ref) {}

	template <class IIt2, class EIt2, class CIIt2, class CEIt2>
	SparseVectorIterator (const SparseVectorIterator<IIt2, EIt2, CIIt2, CEIt2> &i) : _ref (i._ref) {}

	SparseVectorIterator &operator = (const SparseVectorIterator &i)
	{
		_ref.first._i = i._ref.first._i;
		_ref.second._i = i._ref.second._i;
		return *this;
	}

	template <class IIt, class EIt, class CIIt, class CEIt>
	SparseVectorIterator &operator = (const SparseVectorIterator<IIt, EIt, CIIt, CEIt> &i)
	{
		_ref.first._i = i._ref.first._i;
		_ref.second._i = i._ref.second._i;
		return *this;
	}

	SparseVectorIterator &operator ++ () 
	{
		++_ref.first._i;
		++_ref.second._i;
		return *this;
	}

	SparseVectorIterator operator ++ (int) 
	{
		SparseVectorIterator tmp (*this);
		++*this;
		return tmp;
	}

	SparseVectorIterator operator + (difference_type i) const
		{ return SparseVectorIterator (_ref.first._i + i, _ref.second._i + i); }

	SparseVectorIterator &operator += (difference_type i) 
	{
		_ref.first._i += i;
		_ref.second._i += i;
		return *this;
	}

	SparseVectorIterator &operator -- () 
	{
		--_ref.first._i;
		--_ref.second._i;
		return *this;
	}

	SparseVectorIterator operator -- (int) 
	{
		SparseVectorIterator tmp (*this);
		--*this;
		return tmp;
	}

	SparseVectorIterator operator - (difference_type i) const
		{ return *this + -i; }

	SparseVectorIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SparseVectorIterator i) const 
		{ return _ref.first._i - i._ref.first._i; }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ return _ref; }

	const_reference operator * () const 
		{ return _ref; }

	pointer operator -> () 
		{ return &_ref; }

	const_pointer operator -> () const
		{ return reinterpret_cast<const_pointer> (&_ref); }

	bool operator == (const SparseVectorIterator &c) const 
		{ return (_ref.first._i == c._ref.first._i); }

	bool operator < (const SparseVectorIterator &c) const 
		{ return (_ref.first._i < c._ref.first._i); }

	bool operator != (const SparseVectorIterator &c) const 
		{ return (_ref.first._i != c._ref.first._i); }

	template <class IIt, class EIt, class CIIt, class CEIt>
	bool operator == (const SparseVectorIterator<IIt, EIt, CIIt, CEIt> &c) const 
		{ return (_ref.first._i == c._ref.first._i); }

	template <class IIt, class EIt, class CIIt, class CEIt>
	bool operator < (const SparseVectorIterator<IIt, EIt, CIIt, CEIt> &c) const 
		{ return (_ref.first._i < c._ref.first._i); }

	template <class IIt, class EIt, class CIIt, class CEIt>
	bool operator != (const SparseVectorIterator<IIt, EIt, CIIt, CEIt> &c) const 
		{ return (_ref.first._i != c._ref.first._i); }

private:
	template <class IIt, class EIt, class CIIt, class CEIt>
	friend class SparseVectorIterator;

	template <class Element, class IV, class EV>
	friend class SparseVector;

	reference _ref;
};

/** Sparse vector wrapper
 *
 * This class acts as a wrapper around a pair of vectors -- the first
 * for column-indices, the second for entries -- making it appear to
 * be a vector of (index, entry)-pairs.
 *
 * \ingroup vector
 */
template <class IndexIterator, class ElementIterator, class ConstIndexIterator = IndexIterator, class ConstElementIterator = ElementIterator>
class ConstSparseVector
{
public:
	typedef VectorRepresentationTypes::Sparse RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef ConstSparseVector ContainerType;
	typedef SparseSubvector<ConstSparseVector, VectorRepresentationTypes::Sparse> ConstSubvectorType;
	typedef SparseSubvector<ConstSparseVector, VectorRepresentationTypes::Sparse> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef SparseVectorIterator<IndexIterator, ElementIterator, ConstIndexIterator, ConstElementIterator> iterator;
	typedef SparseVectorIterator<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> const_iterator;

	typedef IndexIterator index_iterator;
	typedef ElementIterator element_iterator;
	typedef ConstIndexIterator const_index_iterator;
	typedef ConstElementIterator const_element_iterator;

	typedef typename iterator::reference reference;
	typedef typename iterator::const_reference const_reference;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef typename iterator::value_type value_type;
	typedef size_t size_type;

	template <class IndexVector, class ElementVector>
	ConstSparseVector (IndexVector &iv, ElementVector &ev)
		: _idx_begin (iv.begin ()), _idx_end (iv.end ()), _elt_begin (ev.begin ()) {}

	ConstSparseVector (IndexIterator &idx_begin, IndexIterator &idx_end, ElementIterator &elt_begin)
		: _idx_begin (idx_begin), _idx_end (idx_end), _elt_begin (elt_begin) {}

	ConstSparseVector () {}

	template <class IIt, class EIt, class CIIt, class CEIt>
	inline ConstSparseVector &operator = (const ConstSparseVector<IIt, EIt, CIIt, CEIt> &v)
	{
		_idx_begin = v._idx_begin;
		_idx_end = v._idx_end;
		_elt_begin = v._elt_begin;
		return *this;
	}

	inline iterator               begin  (void)       { return iterator (_idx_begin, _elt_begin); }
	inline const_iterator         begin  (void) const { return const_iterator (_idx_begin, _elt_begin); }
	inline iterator               end    (void)       { return iterator (_idx_end, _elt_begin + size ()); }
	inline const_iterator         end    (void) const { return const_iterator (_idx_end, _elt_begin + size ()); }

	inline reverse_iterator       rbegin (void)       { return reverse_iterator (_idx_end, _elt_begin + size ()); }
	inline const_reverse_iterator rbegin (void) const { return const_reverse_iterator (_idx_end, _elt_begin + size ()); }
	inline reverse_iterator       rend   (void)       { return reverse_iterator (_idx_begin, _elt_begin); }
	inline const_reverse_iterator rend   (void) const { return const_reverse_iterator (_idx_begin, _elt_begin); }

	inline reference       operator[] (size_type n)       { return *(begin () + n); }
	inline const_reference operator[] (size_type n) const { return *(begin () + n); }

	inline reference at (size_type n)
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

	inline reference       front     ()       { return *(begin ()); }
	inline const_reference front     () const { return *(begin ()); }
	inline reference       back      ()       { return *(end () - 1); }
	inline const_reference back      () const { return *(end () - 1); }

	inline size_type       size      () const { return _idx_end - _idx_begin;  }
	inline bool            empty     () const { return _idx_end == _idx_begin; }
	inline size_type       max_size  () const { return _idx_end - _idx_begin;  }

	inline bool operator == (const ConstSparseVector &v) const
		{ return (_idx_begin == v._idx_begin) && (_idx_end == v._idx_end) && (_elt_begin == v._elt_begin); }

protected:
	template <class V, class T>
	friend class SparseSubvector;

	template <class Element, class IndexVector, class ElementVector>
	friend class SparseVector;

	ConstIndexIterator index_begin () const { return _idx_begin; }
	ConstIndexIterator index_end () const { return _idx_end; }
	ConstElementIterator element_begin () const { return _elt_begin; }

	IndexIterator _idx_begin, _idx_end;
	ElementIterator _elt_begin;
};

/** Generic sparse vector
 *
 * This class represents a sparse vector (stored as a pair of vectors)
 * as a vector of (column-index, entry)-pairs.
 *
 * \ingroup vector
 */
template <class Element, class IndexVector, class ElementVector> // N.B. default argument in forward-declaration in traits.h
class SparseVector
{
public:
	typedef VectorRepresentationTypes::Sparse RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef SparseVector ContainerType;
	typedef SparseSubvector<SparseVector, VectorRepresentationTypes::Sparse> SubvectorType;
	typedef SparseSubvector<const SparseVector, VectorRepresentationTypes::Sparse> ConstSubvectorType;
	typedef SparseSubvector<SparseVector, VectorRepresentationTypes::Sparse> AlignedSubvectorType;
	typedef SparseSubvector<const SparseVector, VectorRepresentationTypes::Sparse> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef typename IndexVector::iterator IndexIterator;
	typedef typename IndexVector::const_iterator ConstIndexIterator;
	typedef typename ElementVector::iterator ElementIterator;
	typedef typename ElementVector::const_iterator ConstElementIterator;

	typedef SparseVectorIterator<IndexIterator, ElementIterator, ConstIndexIterator, ConstElementIterator> iterator;
	typedef SparseVectorIterator<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> const_iterator;

	typedef typename iterator::reference reference;
	typedef typename iterator::const_reference const_reference;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef typename iterator::value_type value_type;
	typedef size_t size_type;

	SparseVector () {}
		
	template <class IV, class EV>
	SparseVector (IV &iv, EV &ev)
		: _idx (iv.size ()), _elt (ev.size ())
	{
		std::copy (iv.begin (), iv.end (), _idx.begin ());
		std::copy (ev.begin (), ev.end (), _elt.begin ());
	}

	template <class IIt, class EIt>
	SparseVector (IIt idx_begin, IIt idx_end, EIt elt_begin)
		: _idx (idx_end - idx_begin), _elt (idx_end - idx_begin)
	{
		std::copy (idx_begin, idx_end, _idx.begin ());
		std::copy (elt_begin, elt_begin + (idx_end - idx_begin), _elt.begin ());
	}

	inline SparseVector &operator = (const SparseVector &v)
	{
		_idx = v._idx;
		_elt = v._elt;
		return *this;
	}

	template <class IIt, class EIt, class CIIt, class CEIt>
	inline SparseVector &operator = (const ConstSparseVector<IIt, EIt, CIIt, CEIt> &v)
	{
		_idx.resize (v.size ());
		_elt.resize (v.size ());
		std::copy (v._idx_begin, v._idx_end, _idx.begin ());
		std::copy (v._elt_begin, v._elt_begin + v.size (), _elt.begin ());
		return *this;
	}

	inline iterator               begin  ()       { return iterator (_idx.begin (), _elt.begin ()); }
	inline const_iterator         begin  () const { return const_iterator (_idx.begin (), _elt.begin ()); }
	inline iterator               end    ()       { return iterator (_idx.end (), _elt.end ()); }
	inline const_iterator         end    () const { return const_iterator (_idx.end (), _elt.end ()); }

	inline reverse_iterator       rbegin ()       { return reverse_iterator (_idx.end (), _elt.end ()); }
	inline const_reverse_iterator rbegin () const { return const_reverse_iterator (_idx.end (), _elt.end ()); }
	inline reverse_iterator       rend   ()       { return reverse_iterator (_idx.begin (), _elt.begin ()); }
	inline const_reverse_iterator rend   () const { return const_reverse_iterator (_idx.begin (), _elt.begin ()); }

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
	inline void            push_back (T v)         { _idx.push_back (v.first); _elt.push_back (v.second); }
	inline void            clear     ()            { _idx.clear (); _elt.clear (); }
	inline void            resize    (size_type s) { _idx.resize (s); _elt.resize (s); }

	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
	{
		clear ();

		while (first != last)
			push_back (*first++);
	}

	void assign (const_iterator first, const_iterator last)
	{
		_idx.assign (first._ref.first._i, last._ref.first._i);
		_elt.assign (first._ref.second._i, last._ref.second._i);
	}

	template <class T>
	inline iterator insert (iterator pos, const T &x)
	{
		typename IndexVector::iterator i_idx;
		typename ElementVector::iterator i_elt;

		i_idx = _idx.insert (pos._ref.first._i, x.first); 
		i_elt = _elt.insert (pos._ref.second._i, x.second);

		return iterator (i_idx, i_elt);
	}

	template <class It>
	void insert (iterator pos, It begin, It end)
	{
		typename IndexVector::value_type idx;
		typename ElementVector::value_type elt;

		typename IndexVector::difference_type p = pos._ref.first._i - _idx.begin ();

		_idx.insert (pos._ref.first._i, end - begin, idx);
		_elt.insert (pos._ref.second._i, end - begin, elt);

		std::copy (begin, end, iterator (_idx.begin () + p, _elt.begin () + p));
	}

	inline iterator erase (iterator pos)
	{
		typename IndexVector::iterator i_idx;
		typename ElementVector::iterator i_elt;

		i_idx = _idx.erase (pos._ref.first._i); 
		i_elt = _elt.erase (pos._ref.second._i);

		return iterator (i_idx, i_elt);
	}

	inline iterator erase (iterator first, iterator last)
	{
		typename IndexVector::iterator i_idx;
		typename ElementVector::iterator i_elt;

		i_idx = _idx.erase (first._ref.first._i, last._ref.first._i); 
		i_elt = _elt.erase (first._ref.second._i, last._ref.second._i);

		return iterator (i_idx, i_elt);
	}

	inline size_type       size      () const      { return _idx.size ();  }
	inline bool            empty     () const      { return _idx.empty (); }
	inline size_type       max_size  () const      { return std::min (_idx.max_size (), _elt.max_size ());  }

	inline bool operator == (const SparseVector &v) const
		{ return (_idx == v._idx) && (_elt == v._elt); }

	void swap (SparseVector &v)
		{ std::swap (_idx, v._idx); std::swap (_elt, v._elt); }

protected:
	IndexVector _idx;
	ElementVector _elt;
};

} // namespace LELA

namespace std
{

// Specialisation of std::swap to sparse vectors
template <class Element, class IndexVector, class ElementVector>
void swap (LELA::SparseVector<Element, IndexVector, ElementVector> &v1, LELA::SparseVector<Element, IndexVector, ElementVector> &v2)
	{ v1.swap (v2); }

} // namespace std

#endif // __LELA_VECTOR_SPARSE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
