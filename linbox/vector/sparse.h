/* linbox/vector/sparse.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#ifndef __LINBOX_VECTOR_SPARSE_H
#define __LINBOX_VECTOR_SPARSE_H

#include <vector>

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

public:
	typedef typename std::iterator_traits<Iterator>::value_type T;

	SparseVectorProperty () {}
	SparseVectorProperty (Iterator i) : _i (i) {}

	T &operator = (const T &v)
		{ return *_i = v; }

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

	SparseVectorReference &operator = (const SparseVectorReference &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }

	SparseVectorReference &operator = (const SparseVectorReference<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> &r)
		{ *first._i = *r.first._i; *second._i = *r.second._i; return *this; }
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
};

/// Iterator for sparse vectors
template <class IndexIterator, class ElementIterator, class ConstIndexIterator = IndexIterator, class ConstElementIterator = ElementIterator>
class SparseVectorIterator
{
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef SparseVectorReference<IndexIterator, ElementIterator, ConstIndexIterator, ConstElementIterator> reference;
	typedef reference value_type;
	typedef const SparseVectorReference<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> const_reference;
	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef long difference_type;
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

	value_type _ref;
};

/// Forward declaration
template <class Element, class IndexVector = std::vector<size_t>, class ElementVector = std::vector<Element> >
class SparseVector;

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
	typedef SparseVectorIterator<IndexIterator, ElementIterator, ConstIndexIterator, ConstElementIterator> iterator;
	typedef SparseVectorIterator<ConstIndexIterator, ConstElementIterator, ConstIndexIterator, ConstElementIterator> const_iterator;

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

	inline ConstSparseVector &operator = (const ConstSparseVector &v)
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

	inline reference       at (size_type n)
	{
		if (n >= size ())
			throw std::out_of_range ("LinBox::SparseVector");
		else
			return (*this)[n];
	}

	inline const_reference at (size_type n) const
	{
		if (n >= size ())
			throw std::out_of_range ("LinBox::SparseVector");
		else
			return (*this)[n];
	}

	inline reference       front     (void)       { return *(begin ()); }
	inline const_reference front     (void) const { return *(begin ()); }
	inline reference       back      (void)       { return *(end () - 1); }
	inline const_reference back      (void) const { return *(end () - 1); }

	inline size_type       size      (void) const { return _idx_end - _idx_begin;  }
	inline bool            empty     (void) const { return _idx_end == _idx_begin; }
	inline size_type       max_size  (void) const { return _idx_end - _idx_begin;  }

	inline bool operator == (const ConstSparseVector &v) const
		{ return (_idx_begin == v._idx_begin) && (_idx_end == v._idx_end) && (_elt_begin == v._elt_begin); }

private:
	template <class Element, class IndexVector, class ElementVector>
	friend class SparseVector;

	IndexIterator _idx_begin, _idx_end;
	ElementIterator _elt_begin;
};

/** Generic sparse vector
 *
 * This class represents a sparse vector (stored as a pair of vectors)
 * as a vector of (column-index, entry)-pairs.
 */
template <class Element, class IndexVector, class ElementVector>
class SparseVector
{
public:
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

	SparseVector (IndexIterator &idx_begin, IndexIterator &idx_end, ElementIterator &elt_begin)
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
			throw std::out_of_range ("LinBox::SparseVector");
		else
			return (*this)[n];
	}

	inline const_reference at (size_type n) const
	{
		if (n >= size ())
			throw std::out_of_range ("LinBox::SparseVector");
		else
			return (*this)[n];
	}

	inline reference       front     ()            { return *(begin ()); }
	inline const_reference front     () const      { return *(begin ()); }
	inline reference       back      ()            { return *(end () - 1); }
	inline const_reference back      () const      { return *(end () - 1); }

	template <class T>
	inline void            push_back (const T &v)  { _idx.push_back (v.first); _elt.push_back (v.second); }
	inline void            clear     ()            { _idx.clear (); _elt.clear (); }
	inline void            resize    (size_type s) { _idx.resize (s); _elt.resize (s); }

	inline size_type       size      () const      { return _idx.size ();  }
	inline bool            empty     () const      { return _idx.empty (); }
	inline size_type       max_size  () const      { return std::min (_idx.max_size (), _elt.max_size ());  }

	inline bool operator == (const SparseVector &v) const
		{ return (_idx == v._idx) && (_elt == v._elt); }

private:
	IndexVector _idx;
	ElementVector _elt;
};

} // namespace LinBox

#endif // __LINBOX_VECTOR_SPARSE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
