/* linbox/vector/sparse.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#ifndef __LINBOX_VECTOR_SPARSE_H
#define __LINBOX_VECTOR_SPARSE_H

#include <vector>
#include <algorithm>

#include "linbox/vector/traits.h"

namespace LinBox
{

/// Class to mimic a C#-style property
template <class Iterator>
class SparseVectorProperty {
	Iterator _i;

	template <class IIt, class EIt, class CIIt, class CEIt>
	friend class SparseVectorIterator;

	template <class IIt, class EIt, class CIIt, class CEIt>
	friend class SparseVectorReference;

	template <class Element, class IV, class EV>
	friend class SparseVector;

public:
	typedef typename std::iterator_traits<Iterator>::value_type value_type;

	SparseVectorProperty () {}
	SparseVectorProperty (Iterator i) : _i (i) {}

	SparseVectorProperty &operator = (const value_type &v)
		{ *_i = v; return *this; }

	operator typename std::iterator_traits<Iterator>::reference ()
		{ return *_i; }

	operator const typename std::iterator_traits<Iterator>::reference () const
		{ return *_i; }
};

/// Reference to an (index, entry)-pair
template <class IndexIterator, class ElementIterator, class ConstIndexIterator, class ConstElementIterator>
class SparseVectorReference {
public:
	typedef typename std::iterator_traits<IndexIterator>::value_type first_type;
	typedef typename std::iterator_traits<ElementIterator>::value_type second_type;

	SparseVectorProperty<IndexIterator> first;
	SparseVectorProperty<ElementIterator> second;

	SparseVectorReference () {}
	SparseVectorReference (IndexIterator idx, ElementIterator elt) : first (idx), second (elt) {}

	SparseVectorReference &operator = (const std::pair<first_type, second_type> &r)
		{ *first._i = r.first; *second._i = r.second; return *this; }

	SparseVectorReference &operator = (const SparseVectorReference &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }

	SparseVectorReference &operator = (const SparseVectorReference<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }

	operator std::pair<first_type, second_type> () const
		{ return std::pair<first_type, second_type> (*first._i, *second._i); }
};

/// Specialisation of above for const references
template <class ConstIndexIterator, class ConstElementIterator>
class SparseVectorReference<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> {
public:
	typedef typename std::iterator_traits<ConstIndexIterator>::value_type first_type;
	typedef typename std::iterator_traits<ConstElementIterator>::value_type second_type;

	SparseVectorProperty<ConstIndexIterator> first;
	SparseVectorProperty<ConstElementIterator> second;

	SparseVectorReference () {}
	SparseVectorReference (ConstIndexIterator idx, ConstElementIterator elt) : first (idx), second (elt) {}

	SparseVectorReference &operator = (const SparseVectorReference &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }

	bool operator == (const SparseVectorReference &r) const
		{ return (first == r.first) && (second == r.second); }

	bool operator != (const SparseVectorReference &r) const
		{ return (first != r.first) || (second != r.second); }

private:
	template <class Element, class IV, class EV>
	friend class SparseVector;
};

/// Iterator for sparse vectors
template <class IndexIterator, class ElementIterator, class ConstIndexIterator = IndexIterator, class ConstElementIterator = ElementIterator>
class SparseVectorIterator
{
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef SparseVectorReference<IndexIterator, ElementIterator, ConstIndexIterator, ConstElementIterator> reference;
	typedef std::pair<typename IndexIterator::value_type, typename ElementIterator::value_type> value_type;
	typedef const SparseVectorReference<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> const_reference;
	typedef reference *pointer;
	typedef const_reference *const_pointer;
	typedef ptrdiff_t difference_type;
	typedef size_t size_type;

	SparseVectorIterator () {}
	SparseVectorIterator (IndexIterator idx, ElementIterator elt) : _ref (idx, elt) {}
	SparseVectorIterator (const SparseVectorIterator &i) : _ref (i._ref) {}

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
		{ return &_ref; }

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

/// Forward declarations
template <class Vector, class Trait>
class SparseSubvector;

/** Sparse vector wrapper
 *
 * This class acts as a wrapper around a pair of vectors -- the first
 * for column-indices, the second for entries -- making it appear to
 * be a vector of (index, entry)-pairs.
 */
template <class IndexIterator, class ElementIterator, class ConstIndexIterator = IndexIterator, class ConstElementIterator = ElementIterator>
class ConstSparseVector
{
public:
	typedef typename VectorCategories::SparseVectorTag VectorCategory; 

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

	template <class Vector, class Trait>
	class SparseSubvector;

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
 */
template <class Element, class IndexVector, class ElementVector> // N.B. default argument in forward-declaration in traits.h
class SparseVector
{
public:
	typedef typename VectorCategories::SparseVectorTag VectorCategory; 

	typedef typename IndexVector::iterator IndexIterator;
	typedef typename IndexVector::const_iterator ConstIndexIterator;
	typedef typename ElementVector::iterator ElementIterator;
	typedef typename ElementVector::const_iterator ConstElementIterator;

	typedef IndexIterator index_iterator;
	typedef ElementIterator element_iterator;
	typedef ConstIndexIterator const_index_iterator;
	typedef ConstElementIterator const_element_iterator;

	typedef SparseVectorIterator<IndexIterator, ElementIterator, ConstIndexIterator, ConstElementIterator> iterator;
	typedef SparseVectorIterator<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> const_iterator;

	typedef typename iterator::reference reference;
	typedef typename iterator::const_reference const_reference;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	typedef IndexVector index_vector;
	typedef typename IndexVector::value_type index_type;

	typedef ElementVector element_vector;
	typedef Element element_type;

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

	inline size_type       size      () const      { return _idx.size ();  }
	inline bool            empty     () const      { return _idx.empty (); }
	inline size_type       max_size  () const      { return std::min (_idx.max_size (), _elt.max_size ());  }

	inline bool operator == (const SparseVector &v) const
		{ return (_idx == v._idx) && (_elt == v._elt); }

	void swap (SparseVector &v)
		{ std::swap (_idx, v._idx); std::swap (_elt, v._elt); }

protected:
	template <class V, class T>
	friend class SparseSubvector;

	IndexIterator        index_begin   ()       { return _idx.begin (); }
	IndexIterator        index_end     ()       { return _idx.end (); }
	ConstIndexIterator   index_begin   () const { return _idx.begin (); }
	ConstIndexIterator   index_end     () const { return _idx.end (); }
	ElementIterator      element_begin ()       { return _elt.begin (); }
	ConstElementIterator element_begin () const { return _elt.begin (); }

	IndexVector _idx;
	ElementVector _elt;
};

} // namespace LinBox

namespace std
{

// Specialisation of std::swap to sparse vectors
template <class Element, class IndexVector, class ElementVector>
void swap (LinBox::SparseVector<Element, IndexVector, ElementVector> &v1, LinBox::SparseVector<Element, IndexVector, ElementVector> &v2)
	{ v1.swap (v2); }

} // namespace std

#endif // __LINBOX_VECTOR_SPARSE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
