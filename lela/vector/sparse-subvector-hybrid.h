/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* lela/vector/sparse-subvector-hybrid.h
 * Copyright 2011 Bradford Hovinen
 *
 * Evolved from sparse-subvector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_SPARSE_SUBVECTOR_HYBRID_H
#define __VECTOR_SPARSE_SUBVECTOR_HYBRID_H

#include <vector>

#include "lela/vector/traits.h"
#include "lela/vector/sparse-subvector.h"
#include "lela/vector/hybrid.h"

namespace LELA
{

// Specialisation of SparseSubvector to const vector of hybrid zero-one format

template <class Vector>
class SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01>
{
    public:
	typedef VectorRepresentationTypes::Hybrid01 RepresentationType; 
	typedef VectorStorageTypes::Transformed StorageType;
	typedef HybridVector<typename Vector::Endianness, typename Vector::index_type, typename Vector::word_type> ContainerType;
	typedef SparseSubvector<ContainerType, VectorRepresentationTypes::Hybrid01> SubvectorType;
	typedef SparseSubvector<ContainerType, VectorRepresentationTypes::Hybrid01> ConstSubvectorType;
	typedef SparseSubvector<ContainerType, VectorRepresentationTypes::Hybrid01> AlignedSubvectorType;
	typedef SparseSubvector<ContainerType, VectorRepresentationTypes::Hybrid01> ConstAlignedSubvectorType;
	static const int align = 1;

	typedef typename Vector::index_type index_type;
	typedef typename Vector::word_type word_type;
	typedef std::pair<index_type, word_type> value_type;
	typedef typename Vector::Endianness Endianness;

	class const_iterator;
	friend class const_iterator;

	SparseSubvector () {}
	SparseSubvector (const Vector &v, size_t start, size_t finish)
		: _start (start), _finish (finish)
		{ set_start_end (v.begin (), v.end ()); }
	SparseSubvector (const SparseSubvector &v, size_t start, size_t finish)
		: _start (v._start + start), _finish (v._start + finish)
		{ set_start_end (v._begin, v._end); }
	SparseSubvector (const Vector &v, const_iterator &iter, size_t start, size_t finish)
		: _start (start), _finish (finish), _begin (iter._pos), _end (v.end ())
	{
		if (_begin != v.begin ()) {
			typename Vector::const_iterator t = _begin;

			--t;

			if (start >> WordTraits<word_type>::logof_size == t->first)
			    _begin = t;
		}
	}

	~SparseSubvector () {}

	inline const_iterator begin () const
		{ return const_iterator (_start, _finish, _begin, false); }
	inline const_iterator end   () const
		{ return const_iterator (_start, _finish, _end, true); }

	//	inline size_t         size  () const { return _end - _begin; }
	inline bool           empty () const { return begin () == end (); }

	inline value_type     front () const { return *(begin ()); }

    private:
	size_t _start, _finish;
	typename Vector::const_iterator _begin, _end;

	void set_start_end (typename Vector::const_iterator begin, typename Vector::const_iterator end);

    public:
	class const_iterator
	{
	public:
		typedef std::forward_iterator_tag iterator_category;
		typedef typename Vector::index_type index_type;
		typedef typename Vector::word_type word_type;
		typedef typename std::pair<index_type, word_type> value_type;
		typedef typename std::pair<index_type, word_type> *pointer;
		typedef typename std::iterator_traits<typename Vector::const_iterator>::difference_type difference_type;
		typedef typename Vector::Endianness Endianness;

		typedef SparseSubvector<Vector, VectorRepresentationTypes::Hybrid01> container_type;

		class const_reference
		{
		public:
			const_reference () {}
			const_reference (const const_reference &r)
				: first (r.first), second (r.second) {}

			index_type first;
			word_type second;

			operator value_type () const { return value_type (first, second); }
		};

		typedef const_reference reference;

		friend class SparseSubvector;

		const_iterator () {}
		const_iterator (index_type start_idx, index_type end_idx, const typename Vector::const_iterator &pos, bool end)
			: _pos (pos), _start_idx (start_idx), _end_idx (end_idx), _ref_first_valid (false), _ref_second_valid (false), _start (true), _end (end)
			{}

		const_iterator (const const_iterator &i)
			: _pos (i._pos), _start_idx (i._start_idx), _end_idx (i._end_idx),
			  _ref (i._ref), _ref_first_valid (i._ref_first_valid), _ref_second_valid (i._ref_second_valid), _start (i._start), _end (i._end) {}

		const_iterator &operator = (const const_iterator &i)
		{
			_pos = i._pos;
			_start_idx = i._start_idx;
			_end_idx = i._end_idx;
			_ref.first = i._ref.first;
			_ref.second = i._ref.second;
			_ref_first_valid = i._ref_first_valid;
			_ref_second_valid = i._ref_second_valid;
			_start = i._start;
			_end = i._end;
			return *this;
		}

		const_iterator &operator ++ ();

		const_iterator operator ++ (int) 
		{
			const_iterator tmp (*this);
			++*this;
			return tmp;
		}

		difference_type operator - (const_iterator &i) const 
			{ return _pos - i._pos; }

		const const_reference &operator * ()
			{ update_ref (); return _ref; }

		const const_reference *operator -> ()
			{ update_ref (); return &_ref; }

		bool operator == (const const_iterator &c);

		bool operator != (const const_iterator &c)
			{ return !(*this == c); }

	private:
		typename Vector::const_iterator _pos;
		index_type _start_idx;
		index_type _end_idx;
		const_reference _ref;
		bool _ref_first_valid;
		bool _ref_second_valid;
		bool _start;
		bool _end;

		void update_ref ();
	};
}; // template <class Vector> class SparseSubvector<Vector, Hybrid01>

} // namespace LELA

#include "lela/vector/sparse-subvector-hybrid.tcc"

#endif // __VECTOR_SPARSE_SUBVECTOR_HYBRID_H
