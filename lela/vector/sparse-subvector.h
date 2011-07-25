/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* lela/vector/sparse-subvector.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from sparse-vector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SPARSE_SUBVECTOR_H
#define __VECTOR_SPARSE_SUBVECTOR_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "lela/vector/traits.h"
#include "lela/vector/shifted-vector.h"
#include "lela/vector/sparse.h"
#include "lela/vector/subvector.h"

namespace LELA
{

/** Empty struct to be used when constructing a sparse subvector of a
 * hybrid vector with word-aligned boundaries. This version is much
 * faster.
 *
 * \ingroup vector
 */
struct HybridSubvectorWordAlignedTag {};

/** A subvector of a sparse vector. It provides an interface which
 * mimics that of the containing vector, but does not allow any
 * operations which modify the vector's contents.
 *
 * \ingroup vector
 */

template <class Vector, class Trait>
class SparseSubvector
{
    public:
	typedef Trait RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef Vector ContainerType;
	typedef SparseSubvector SubvectorType;
	typedef SparseSubvector<const Vector, Trait> ConstSubvectorType;
	typedef SparseSubvector AlignedSubvectorType;
	typedef SparseSubvector<const Vector, Trait> ConstAlignedSubvectorType;
	static const int align = 1;

	SparseSubvector ();
	SparseSubvector (Vector &v, size_t start, size_t end);
	SparseSubvector (SparseSubvector &v, size_t start, size_t end);
	SparseSubvector (const SparseSubvector &v);

	~SparseSubvector () {}
}; // template <class Vector, Trait> class SparseSubvector

// Prototype for reference into sparse subvector
template <class Iterator>
struct SparseSubvectorReferencePT
{
	typedef typename ShiftedPropertyFirstEntry<Iterator>::value_type first_type;
	typedef typename Property<Iterator, SecondEntryAccessor<Iterator> >::value_type second_type;

	// This is idiotic: we must maintain two copies of the same
	// iterator -- one for each property. Ordinarily we would use
	// a union to avoid this, but this is impossible in C++
	// because we cannot guarantee that Iterator does not have a
	// nontrivial copy-constructor, and there is no way around
	// this.
	ShiftedPropertyFirstEntry<Iterator> first;
	Property<Iterator, SecondEntryAccessor<Iterator> > second;

	SparseSubvectorReferencePT () {}

	SparseSubvectorReferencePT (Iterator i, typename std::iterator_traits<Iterator>::value_type::first_type shift)
		: first (i, shift), second (i)
		{}

	SparseSubvectorReferencePT (const SparseSubvectorReferencePT &r)
		{ first._i = r.first._i; first._shift = r.first._shift; second._i = r.second._i; }
};

// Prototype for iterator over sparse subvector
template <class Iterator, class ConstIterator>
class SparseSubvectorIteratorPT
{
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef SparseSubvectorReferencePT<Iterator> reference;
	typedef SparseSubvectorReferencePT<ConstIterator> const_reference;
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	typedef reference *pointer;
	typedef const_reference *const_pointer;
	typedef ptrdiff_t difference_type;
	typedef size_t size_type;

	SparseSubvectorIteratorPT () {}
	SparseSubvectorIteratorPT (Iterator i, typename std::iterator_traits<Iterator>::value_type::first_type shift)
		: _ref (i, shift) {}
	SparseSubvectorIteratorPT (const SparseSubvectorIteratorPT &i)
		: _ref (i._ref) {}

	SparseSubvectorIteratorPT &operator = (const SparseSubvectorIteratorPT &i)
	{
		_ref.first._i = i._ref.first._i;
		_ref.first._shift = i._ref.first._shift;
		_ref.second._i = i._ref.second._i;
		return *this;
	}

	template <class It, class CIt>
	SparseSubvectorIteratorPT &operator = (const SparseSubvectorIteratorPT<It, CIt> &i)
	{
		_ref.first._i = i._ref.first._i;
		_ref.first._shift = i._ref.first._shift;
		_ref.second._i = i._ref.second._i;
		return *this;
	}

	SparseSubvectorIteratorPT &operator ++ () 
		{ ++_ref.first._i; ++_ref.second._i; return *this; }

	SparseSubvectorIteratorPT operator ++ (int) 
	{
		SparseSubvectorIteratorPT tmp (*this);
		++*this;
		return tmp;
	}

	SparseSubvectorIteratorPT operator + (difference_type i) const
		{ return SparseSubvectorIteratorPT (_ref.first._i + i, _ref.first._shift); }

	SparseSubvectorIteratorPT &operator += (difference_type i) 
		{ _ref.first._i += i; _ref.second._i += i; return *this; }

	SparseSubvectorIteratorPT &operator -- () 
		{ --_ref.first._i; --_ref.second._i; return *this; }

	SparseSubvectorIteratorPT operator -- (int) 
	{
		SparseSubvectorIteratorPT tmp (*this);
		--*this;
		return tmp;
	}

	SparseSubvectorIteratorPT operator - (difference_type i) const
		{ return *this + -i; }

	SparseSubvectorIteratorPT &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SparseSubvectorIteratorPT i) const 
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

	bool operator == (const SparseSubvectorIteratorPT &c) const 
		{ return (_ref.first._i == c._ref.first._i); }

	bool operator < (const SparseSubvectorIteratorPT &c) const 
		{ return (_ref.first._i < c._ref.first._i); }

	bool operator != (const SparseSubvectorIteratorPT &c) const 
		{ return (_ref.first._i != c._ref.first._i); }

	template <class It, class CIt>
	bool operator == (const SparseSubvectorIteratorPT<It, CIt> &c) const 
		{ return (_ref.first._i == c._ref.first._i); }

	template <class It, class CIt>
	bool operator < (const SparseSubvectorIteratorPT<It, CIt> &c) const 
		{ return (_ref.first._i < c._ref.first._i); }

	template <class It, class CIt>
	bool operator != (const SparseSubvectorIteratorPT<It, CIt> &c) const 
		{ return (_ref.first._i != c._ref.first._i); }

private:
	template <class It, class CIt>
	friend class SparseSubvectorIteratorPT;

	template <class Vector, class Trait>
	friend class SparseSubvector;

	reference _ref;
};

// Specialisation of SparseSubvector to const vector in sparse format

template <class Vector>
class SparseSubvector<const Vector, VectorRepresentationTypes::Sparse>
{
    public:
	typedef VectorRepresentationTypes::Sparse RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef const Vector ContainerType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> SubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> ConstSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> AlignedSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef typename Vector::value_type value_type;
	typedef typename Vector::size_type size_type;

	typedef SparseSubvectorReferencePT<typename Vector::const_iterator> const_reference;
	typedef SparseSubvectorIteratorPT<typename Vector::const_iterator, typename Vector::const_iterator> const_iterator;

	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		_begin = std::lower_bound (v.begin (), v.end (), start, VectorUtils::FindSparseEntryLB ());
		_end = std::lower_bound (v.begin (), v.end (), finish, VectorUtils::FindSparseEntryLB ());
		_shift = start;
	}

	SparseSubvector (const SparseSubvector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		_begin = std::lower_bound (v._begin, v._end, start + v._shift, VectorUtils::FindSparseEntryLB ());
		_end = std::lower_bound (v._begin, v._end, finish + v._shift, VectorUtils::FindSparseEntryLB ());
		_shift = start + v._shift;
	}

	~SparseSubvector () {}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ _begin = v._begin; _end = v._end; _shift = v._shift; return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<const V, VectorRepresentationTypes::Sparse> &v)
		{ _begin = v._begin; _end = v._end; _shift = v._shift; return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<V, VectorRepresentationTypes::Sparse> &v)
		{ _begin = v._v->begin () + v._begin_idx; _end = v._v->begin () + v._end_idx; _shift = v._shift; return *this; }

	inline const_iterator         begin  () const { return const_iterator (_begin, _shift); }
	inline const_iterator         end    () const { return const_iterator (_end, _shift); }

	inline const_reverse_iterator rbegin () const { return const_reverse_iterator (end ()); }
	inline const_reverse_iterator rend   () const { return const_reverse_iterator (begin ()); }

	inline const_reference front     () const      { return *(begin ()); }
	inline const_reference back      () const      { return *(end () - 1); }

	inline size_type       size      () const { return _end - _begin;  }
	inline bool            empty     () const { return _end == _begin; }

	inline bool operator == (const SparseSubvector &v) const
		{ return (_begin == v._begin) && (_end == v._end) && (_shift == v._shift); }

    protected:

	typename Vector::const_iterator _begin, _end;
	typename Vector::value_type::first_type _shift;

	SparseSubvector (const Vector &v, typename Vector::const_iterator begin, typename Vector::const_iterator end,
			 typename std::iterator_traits<typename Vector::const_iterator>::value_type::first_type start)
		: _begin (begin), _end (end), _shift (start)
	{}

}; // template <class Vector> class SparseSubvector<const Vector, Sparse>

// Specialisation of SparseSubvector to vector in sparse format
//
// Note: only const subvectors of a mutable sparse subvector may be created!

template <class Vector>
class SparseSubvector<Vector, VectorRepresentationTypes::Sparse>
{
    public:
	typedef VectorRepresentationTypes::Sparse RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef const Vector ContainerType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> SubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> ConstSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> AlignedSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef typename Vector::value_type value_type;
	typedef typename Vector::size_type size_type;

	typedef SparseSubvectorReferencePT<typename Vector::iterator> reference;
	typedef SparseSubvectorReferencePT<typename Vector::const_iterator> const_reference;
	typedef SparseSubvectorIteratorPT<typename Vector::iterator, typename Vector::const_iterator> iterator;
	typedef SparseSubvectorIteratorPT<typename Vector::const_iterator, typename Vector::const_iterator> const_iterator;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::value_type::first_type start, typename Vector::value_type::first_type finish)
	{
		_v = &v;
		_begin_idx = std::lower_bound (v.begin (), v.end (), start, VectorUtils::FindSparseEntryLB ()) - v.begin ();
		_end_idx = std::lower_bound (v.begin (), v.end (), finish, VectorUtils::FindSparseEntryLB ()) - v.begin ();
		_shift = start;
	}

	SparseSubvector (const SparseSubvector &v)
		: _v (v._v), _begin_idx (v._begin_idx), _end_idx (v._end_idx), _shift (v._shift) {}

	~SparseSubvector () {}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ _v = v._v; _begin_idx = v._begin_idx; _end_idx = v._end_idx; _shift = v._shift; return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<V, VectorRepresentationTypes::Sparse> &v)
		{ _v = v._v; _begin_idx = v._begin_idx; _end_idx = v._end_idx; _shift = v._shift; return *this; }

	inline iterator               begin  ()       { return iterator (_v->begin () + _begin_idx, _shift); }
	inline const_iterator         begin  () const { return const_iterator (_v->begin () + _begin_idx, _shift); }
	inline iterator               end    ()       { return iterator (_v->begin () + _end_idx, _shift); }
	inline const_iterator         end    () const { return const_iterator (_v->begin () + _end_idx, _shift); }

	inline reverse_iterator       rbegin ()       { return reverse_iterator (end ()); }
	inline const_reverse_iterator rbegin () const { return const_reverse_iterator (end ()); }
	inline reverse_iterator       rend   ()       { return reverse_iterator (begin ()); }
	inline const_reverse_iterator rend   () const { return const_reverse_iterator (begin ()); }

	inline reference       front     ()            { return *(begin ()); }
	inline const_reference front     () const      { return *(begin ()); }
	inline reference       back      ()            { return *(end () - 1); }
	inline const_reference back      () const      { return *(end () - 1); }

	template <class T>
	inline void            push_back (T v)         { insert (end (), v); }
	inline void            clear     ()            { _v->erase (_v->begin () + _begin_idx, _v->begin () + _end_idx); _end_idx = _begin_idx; }

	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
	{
		clear ();

		while (first != last)
			push_back (*first++);
	}

	template <class T>
	inline iterator insert (iterator pos, const T &x)
		{ ++_end_idx; return iterator (_v->insert (pos._ref.first._i, value_type (x.first + _shift, x.second)), _shift); }

	inline iterator erase (iterator pos)
		{ --_end_idx; return iterator (_v->erase (pos._ref._i)); }

	inline size_type       size      () const { return _end_idx - _begin_idx;  }
	inline bool            empty     () const { return _end_idx == _begin_idx; }

	inline bool operator == (const SparseSubvector &v) const
		{ return (_v == v._v) && (_begin_idx == v._begin_idx) && (_end_idx == v._end_idx) && (_shift == v._shift); }

    protected:

	friend class SparseSubvector<const Vector, VectorRepresentationTypes::Sparse>;

	Vector *_v;
	typename Vector::size_type _begin_idx, _end_idx;
	typename Vector::value_type::first_type _shift;

	SparseSubvector (Vector &v, typename Vector::const_iterator begin, typename Vector::const_iterator end,
			 typename std::iterator_traits<typename Vector::const_iterator>::value_type::first_type start)
		: _v (&v), _begin_idx (begin - v.begin ()), _end_idx (end - v.begin ()), _shift (start)
	{}

}; // template <class Vector> class SparseSubvector<Vector, Sparse>

// Specialisation of SparseSubvector to const vector in sparse zero-one format

template <class Vector>
class SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> : public ConstShiftedVector<typename Vector::const_iterator>
{
    public:
	typedef VectorRepresentationTypes::Sparse01 RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef Vector ContainerType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> SubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> ConstSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> AlignedSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef ConstShiftedVector<typename Vector::const_iterator> parent_type;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: parent_type (std::lower_bound (v.begin (), v.end (), start), std::lower_bound (v.begin (), v.end (), finish), start)
		{}
	SparseSubvector (const SparseSubvector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: parent_type (std::lower_bound (v._v.begin (), v._v.end (), start + v._shift), std::lower_bound (v._v.begin (), v._v.end (), finish + v._shift), start + v._shift)
		{}
	SparseSubvector (const SparseSubvector &v)
		: parent_type (v)
		{}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ parent_type::operator = (v); return *this; }

	~SparseSubvector () {}
}; // template <class Vector> class SparseSubvector<const Vector, Sparse01>

// Specialisation of SparseSubvector to vector in sparse zero-one format

template <class Vector>
class SparseSubvector<Vector, VectorRepresentationTypes::Sparse01> : public ShiftedVector<MutableSubvector<Vector> >
{
	MutableSubvector<Vector> _v;

    public:
	typedef VectorRepresentationTypes::Sparse01 RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef Vector ContainerType;
	typedef SparseSubvector SubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> ConstSubvectorType;
	typedef SparseSubvector AlignedSubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Sparse01> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef ShiftedVector<MutableSubvector<Vector> > parent_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: _v (v, std::lower_bound (v.begin (), v.end (), start), std::lower_bound (v.begin (), v.end (), finish))
		{ parent_type::set_parent (_v); parent_type::set_shift (start); }
	SparseSubvector (SparseSubvector &v, typename Vector::value_type start, typename Vector::value_type finish)
		: _v (v._v.parent (), std::lower_bound (v._v.begin (), v._v.end (), start + v._shift), std::lower_bound (v._v.begin (), v._v.end (), finish + v._shift))
		{ parent_type::set_parent (_v); parent_type::set_shift (start + v._shift); }

	SparseSubvector (const SparseSubvector &v)
		: _v (v._v)
		{ parent_type::set_parent (_v); parent_type::set_shift (v._shift); }

	SparseSubvector &operator = (const SparseSubvector &v)
		{ _v = v._v; parent_type::set_shift (v._shift); parent_type::set_parent (_v); return *this; }

	~SparseSubvector () {}
}; // template <class Vector> class SparseSubvector<Vector, Sparse01>

// Specialisation of SparseSubvector to hybrid 0-1 vector; must word aligned

template <class Vector>
class SparseSubvector<Vector, HybridSubvectorWordAlignedTag>
	: public SparseSubvector<Vector, VectorRepresentationTypes::Sparse>
{
    public:
	typedef VectorRepresentationTypes::Hybrid01 RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef Vector ContainerType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Hybrid01> SubvectorType;
	typedef SparseSubvector<const Vector, VectorRepresentationTypes::Hybrid01> ConstSubvectorType;
	typedef SparseSubvector<Vector, HybridSubvectorWordAlignedTag> AlignedSubvectorType;
	typedef SparseSubvector<const Vector, HybridSubvectorWordAlignedTag> ConstAlignedSubvectorType;
	static const int align = WordTraits<typename Vector::word_type>::bits;

	typedef SparseSubvector<Vector, VectorRepresentationTypes::Sparse> parent_type;
	typedef typename Vector::Endianness Endianness;
	typedef typename Vector::index_type index_type;
	typedef typename Vector::word_type word_type;

	SparseSubvector () {}
	SparseSubvector (Vector &v, size_t start, size_t finish)
		: parent_type (v, std::lower_bound (v.begin (), v.end (), start >> WordTraits<word_type>::logof_size, VectorUtils::FindSparseEntryLB ()),
			       std::lower_bound (v.begin (), v.end (), (finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size,
						 VectorUtils::FindSparseEntryLB ()),
			       start >> WordTraits<word_type>::logof_size)
	{
		lela_check ((start & WordTraits<word_type>::pos_mask) == 0);
		lela_check ((finish & WordTraits<word_type>::pos_mask) == 0 || v.empty () || finish > (size_t) (v.back ().first << WordTraits<word_type>::logof_size));
	}

	SparseSubvector (SparseSubvector &v, size_t start, size_t finish)
		: parent_type (v, start >> WordTraits<word_type>::logof_size, (finish + WordTraits<word_type>::bits - 1) >> WordTraits<word_type>::logof_size)
	{
		lela_check ((start & WordTraits<word_type>::pos_mask) == 0);
		lela_check ((finish & WordTraits<word_type>::pos_mask) == 0 || v.empty () || finish > (size_t) (v.back ().first << WordTraits<word_type>::logof_size));
	}

	SparseSubvector &operator = (const SparseSubvector &v)
		{ this->SparseSubvector<Vector, VectorRepresentationTypes::Sparse>::operator = (v); return *this; }

	template <class V>
	SparseSubvector &operator = (const SparseSubvector<V, VectorRepresentationTypes::Sparse> &v)
		{ this->SparseSubvector<Vector, VectorRepresentationTypes::Sparse>::operator = (v); return *this; }

}; // template <class Vector> class SparseSubvector<Vector, HybridSubvectorWordAlignedTag>

} // namespace LELA

#include "lela/vector/sparse-subvector-hybrid.h"

#endif // __VECTOR_SPARSE_SUBVECTOR_H
