/* lela/blas/level1-gf2.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for GF2
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_GF2_TCC
#define __BLAS_LEVEL1_GF2_TCC

#include <algorithm>
#include <iostream>

#include "lela/blas/level1-gf2.h"
#include "lela/blas/level1-ll.h"
#include "lela/blas/level1-generic.tcc"

namespace LELA
{

namespace BLAS1 
{

template <class Modules, class reference, class Vector1, class Vector2>
reference &_dot<GF2, GenericModule<GF2>::Tag>::dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
							 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01)
{
	lela_check (x.size () == y.size ());

	if (x.empty ())
		return res = false;

	typename Vector1::word_type t = 0;
	typename Vector1::const_word_iterator i = x.word_begin ();
	typename Vector2::const_word_iterator j = y.word_begin ();

	if (x.empty ())
		return res = false;
	else if (i == x.word_end ())
		t = x.back_word () & y.back_word ();
	else {
		while (i != x.word_end ())
			t ^= *i++ & *j++;
        
		t ^= x.back_word () & y.back_word ();
	}

        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class Modules, class reference, class Vector1, class Vector2>
reference &_dot<GF2, GenericModule<GF2>::Tag>::dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
							 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01)
{
	lela_check (VectorUtils::hasDim<GF2> (y, x.size ()));

	typename Vector2::const_iterator i = y.begin ();

	res = false;

	for (; i != y.end (); ++i)
		if (x[*i]) res = !res;

	return res;
}

template <class Modules, class reference, class Vector1, class Vector2>
reference &_dot<GF2, GenericModule<GF2>::Tag>::dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
							 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01)
{
	typename Vector1::word_type t = 0;
	typename Vector1::const_word_iterator i = x.word_begin ();
	typename Vector2::const_iterator j;

	for (j = y.begin (); j != y.end (); ++j) {
		if (j->first == x.word_size () - 1)
			t ^= x.back_word () & j->second;
		else
			t ^= *(i + j->first) & j->second;
	}
        
        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class Modules, class reference, class Vector1, class Vector2>
reference &_dot<GF2, GenericModule<GF2>::Tag>::dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
							 VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01)
{
	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	res = false;

	while (i != x.end () || j != y.end ()) {
		while (i != x.end () && (j == y.end () || *i < *j)) ++i;
		while (j != y.end () && (i == x.end () || *j < *i)) ++j;
		if (i != x.end () && j != y.end () && *i == *j) { res = !res; ++i; ++j; }
	}

	return res;
}

template <class Modules, class reference, class Vector1, class Vector2>
reference &_dot<GF2, GenericModule<GF2>::Tag>::dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
							 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	res = false;

	for (i = x.begin (), j = y.begin (); i != x.end () && j != y.end (); ++i) {
		for (; j != y.end () && *j >> WordTraits<typename Vector1::word_type>::logof_size < i->first; ++j);

		for (; j != y.end () && *j >> WordTraits<typename Vector1::word_type>::logof_size == i->first; ++j)
			if (i->second & Vector1::Endianness::e_j (*j & WordTraits<typename Vector1::word_type>::pos_mask))
				res = !res;
	}

	return res;
}

template <class Modules, class reference, class Vector1, class Vector2>
reference &_dot<GF2, GenericModule<GF2>::Tag>::dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
							 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01)
{
	if (x.empty () || y.empty ())
		return res = false;

	typename Vector1::word_type t = 0;
	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	while (i != x.end () && j != y.end ()) {
		while (i != x.end () && i->first < j->first)
			++i;

		if (i == x.end ())
			break;

		while (j != y.end () && j->first < i->first)
			++j;

		t ^= i->second & j->second;
		++i; ++j;
	}

        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<GF2, GenericModule<GF2>::Tag>::copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01)
{
	typedef typename std::iterator_traits<typename Vector2::const_iterator>::value_type index_type;

	typename Vector1::const_iterator i;
	index_type idx = 0;

	y.clear ();

	for (i = x.begin (); i != x.end (); ++i, ++idx)
		if (*i) y.push_back (idx);

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<GF2, GenericModule<GF2>::Tag>::copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01)
{
	typename Vector1::const_word_iterator i;
	typename Vector2::index_type idx = 0;

	y.clear ();

	for (i = x.word_begin (); i != x.word_end (); ++i, ++idx)
		if (*i)
			y.push_back (typename Vector2::value_type (idx, *i));

	if (x.back_word ())
		y.push_back (typename Vector2::value_type (idx, x.back_word ()));

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<GF2, GenericModule<GF2>::Tag>::copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01)
{
	lela_check (VectorUtils::hasDim<GF2> (x, y.size ()));

	typename Vector1::const_iterator i;

	std::fill (y.word_begin (), y.word_end (), 0);
	y.back_word () = 0;

	for (i = x.begin (); i != x.end (); ++i)
        	y[*i] = true;

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<GF2, GenericModule<GF2>::Tag>::copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01)
{
	typedef typename std::iterator_traits<typename Vector1::const_iterator>::value_type index_type;

	typename Vector1::const_iterator i;

	y.clear ();

	for (i = x.begin (); i != x.end (); ++i) {
		if (y.empty () || y.back ().first != *i >> WordTraits<typename Vector2::word_type>::logof_size)
			y.push_back (typename Vector2::value_type (*i >> WordTraits<typename Vector2::word_type>::logof_size, 0));

		y.back ().second |= Vector2::Endianness::e_j (*i & WordTraits<typename Vector2::word_type>::pos_mask);
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<GF2, GenericModule<GF2>::Tag>::copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01)
{
	typename Vector1::const_iterator i;

	std::fill (y.word_begin (), y.word_end (), 0);
	y.back_word () = 0;

	for (i = x.begin (); i != x.end (); ++i) {
		if (i->first < y.size () >> WordTraits<typename Vector2::word_type>::logof_size)
			*(y.word_begin () + i->first) = i->second;
		else
			y.back_word () = i->second;
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<GF2, GenericModule<GF2>::Tag>::copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01)
{
	typedef typename std::iterator_traits<typename Vector2::iterator>::value_type index_type;

	typename Vector1::const_iterator i;
	index_type idx = 0;
	typename Vector1::word_type t;

	y.clear ();

	for (i = x.begin (); i != x.end (); ++i)
		for (t = Vector1::Endianness::e_0, idx = i->first << WordTraits<typename Vector1::word_type>::logof_size; t != 0; t = Vector1::Endianness::shift_right (t, 1), ++idx)
			if (i->second & t) y.push_back (idx);

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, typename GenericModule<GF2>::Tag>::axpy_impl
	(const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Generic)
{
	if (a) {
		typename Vector<GF2>::Dense tmp (x.size ());

		_copy<GF2, typename Modules::Tag>::op (F, M, x, tmp);
		_axpy<GF2, typename Modules::Tag>::op (F, M, F.one (), y, tmp);
		_copy<GF2, typename Modules::Tag>::op (F, M, tmp, y);
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, GenericModule<GF2>::Tag>::axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01)
{
	lela_check (y.size () == x.size ());

	if (a) {
		typename Vector2::word_iterator i = y.word_begin ();
		typename Vector1::const_word_iterator j = x.word_begin ();

		for (; i != y.word_end (); ++i, ++j)
			*i ^= *j;

		y.back_word () ^= x.back_word ();
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, GenericModule<GF2>::Tag>::axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01)
{
	lela_check (VectorUtils::hasDim<GF2> (x, y.size ()));

	if (a) {
		typename Vector1::const_iterator i;

		for (i = x.begin (); i != x.end (); ++i)
			y[*i] = !y[*i];
	}

	return y;
}

template <class Modules, class Vector>
void fast_copy (const GF2 &F, Modules &M, std::vector<typename Vector::value_type> &v, Vector &w)
	{ _copy<GF2, typename Modules::Tag>::op (F, M, v, w); }

template <class Modules, class Vector>
void fast_copy (const GF2 &F, Modules &M, HybridVector<typename Vector::Endianness, typename Vector::index_type, typename Vector::word_type> &v, Vector &w)
	{ _copy<GF2, typename Modules::Tag>::op (F, M, v, w); }

template <class Modules, class index_type>
void fast_copy (const GF2 &F, Modules &M, std::vector<index_type> &v, std::vector<index_type> &w)
	{ std::swap (v, w); }

template <class Modules, class Endianness, class index_type, class word_type>
void fast_copy (const GF2 &F, Modules &M, HybridVector<Endianness, index_type, word_type> &v, HybridVector<Endianness, index_type, word_type> &w)
	{ std::swap (v, w); }

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, GenericModule<GF2>::Tag>::axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01)
{
	if (a) {
		std::vector<typename Vector2::value_type> res;

		typename Vector2::const_iterator i = y.begin ();
		typename Vector1::const_iterator j = x.begin ();

		while (i != y.end () || j != x.end ()) {
			while (i != y.end () && (j == x.end () || *i < *j)) { res.push_back (*i); ++i; }
			while (j != x.end () && (i == y.end () || *j < *i)) { res.push_back (*j); ++j; }
			if (i != y.end () && j != x.end () && *i == *j) { ++i; ++j; }
		}

		fast_copy (F, M, res, y);
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, typename GenericModule<GF2>::Tag>::axpy_impl
	(const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01)
{
	if (a) {
		typename Vector<GF2>::Hybrid tmp;

		_copy<GF2, typename Modules::Tag>::op (F, M, x, tmp);
		_axpy<GF2, typename Modules::Tag>::op (F, M, F.one (), y, tmp);
		_copy<GF2, typename Modules::Tag>::op (F, M, tmp, y);
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, GenericModule<GF2>::Tag>::axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01)
{
	if (a) {
		typename Vector1::const_iterator i;
		typename Vector2::word_iterator j = y.word_begin ();

		for (i = x.begin (); i != x.end (); ++i) {
			if (i->first < y.size () >> WordTraits<typename Vector2::word_type>::logof_size)
				*(j + i->first) ^= i->second;
			else
				y.back_word () ^= i->second;
		}
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, typename GenericModule<GF2>::Tag>::axpy_impl
	(const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01)
{
	if (a) {
		typename Vector<GF2>::Sparse tmp;

		_copy<GF2, typename Modules::Tag>::op (F, M, x, tmp);
		_axpy<GF2, typename Modules::Tag>::op (F, M, F.one (), y, tmp);
		_copy<GF2, typename Modules::Tag>::op (F, M, tmp, y);
	}

	return y;
}

template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<GF2, GenericModule<GF2>::Tag>::axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
							 VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01)
{
	if (a) {
		HybridVector<typename Vector2::Endianness, typename Vector2::index_type, typename Vector2::word_type> res;

		typename Vector2::iterator i = y.begin ();
		typename Vector1::const_iterator j = x.begin ();

		while (i != y.end () || j != x.end ()) {
			while (i != y.end () && (j == x.end () || i->first < j->first)) {
				res.push_back (typename Vector1::value_type (i->first, i->second));
				++i;
			}
			while (j != x.end () && (i == y.end () || j->first < i->first)) {
				res.push_back (typename Vector1::value_type (j->first, j->second));
				++j;
			}
			if (i != y.end () && j != x.end () && i->first == j->first) {
				if (i->second ^ j->second)
					res.push_back (typename Vector1::value_type (i->first, i->second ^ j->second));

				++i;
				++j;
			}
		}

		fast_copy (F, M, res, y);
	}

	return y;
}

template <class Modules, class Iterator, class Vector>
Vector &_permute<GF2, GenericModule<GF2>::Tag>::permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v,
							      VectorRepresentationTypes::Dense01)
{
	for (; P_begin != P_end; ++P_begin)
		std::swap (v[P_begin->first], v[P_begin->second]);

	return v;
}

template <class Modules, class Iterator, class Vector>
Vector &_permute<GF2, GenericModule<GF2>::Tag>::permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v,
							      VectorRepresentationTypes::Sparse01)
{
	typename Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		*j = image (*j, P_begin, P_end);

	std::sort (v.begin (), v.end ());

	return v;
}

template <class Modules, class Iterator, class Vector>
Vector &_permute<GF2, GenericModule<GF2>::Tag>::permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v,
							      VectorRepresentationTypes::Hybrid01)
{
	::LELA::Vector<GF2>::Sparse w;

	_copy<GF2, typename Modules::Tag>::op (F, M, v, w);
	_permute<GF2, typename Modules::Tag>::op (F, M, P_begin, P_end, w);
	_copy<GF2, typename Modules::Tag>::op (F, M, w, v);

	return v;
}

template <class Modules, class Vector1, class Vector2>
bool _equal<GF2, GenericModule<GF2>::Tag>::equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
						       VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01)
{
	if (x.size () != y.size ())
		return false;

	typename Vector1::const_word_iterator i = x.word_begin ();
	typename Vector2::const_word_iterator j = y.word_begin ();

	for (; j != y.word_end (); ++j, ++i)
		if (*i != *j) return false;

	if (x.back_word () != y.back_word ())
		return false;

	return true;
}

template <class Modules, class Vector1, class Vector2>
bool _equal<GF2, GenericModule<GF2>::Tag>::equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
						       VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01)
{
	if (x.empty () && !y.empty ())
		return false;

	typedef typename std::iterator_traits<typename Vector2::const_iterator>::value_type index_type;

	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();
	index_type idx = 0;

	for (; j != y.end (); ++j, ++i, ++idx) {
		while (idx < *j) {
			if (*i) return false;
			++idx;
			++i;

			if (i == x.end ())
				return false;
		}

		if (!*i) return false;
	}

	for (; i != x.end (); ++i)
		if (*i) return false;

	return true;
}

template <class Modules, class Vector1, class Vector2>
bool _equal<GF2, GenericModule<GF2>::Tag>::equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
						       VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01)
{
	typename Vector1::const_iterator i_1;
	typename Vector2::const_iterator i_2;

	if (x.size () != y.size ())
		return false;

	for (i_1 = x.begin (), i_2 = y.begin (); i_1 != x.end (); ++i_1, ++i_2)
		if (*i_1 != *i_2)
			return false;

	return true;
}

template <class Modules, class Vector1, class Vector2>
bool _equal<GF2, GenericModule<GF2>::Tag>::equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
						       VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01)
{
	if (x.empty () && !y.empty ())
		return false;

	typename Vector1::const_word_iterator i = x.word_begin ();
	typename Vector2::const_iterator j = y.begin ();
	typename Vector2::index_type idx = 0;

	for (; j != y.end (); ++j) {
		while (idx < j->first) {
			if (*i) return false;
			++idx;
			++i;

			if (i == x.word_end () && (j->first != idx || j->second != x.back_word ()))
				return false;
		}

		if (i == x.word_end ())
			return j->second == x.back_word ();

		if (*i++ != j->second)
			return false;

		++idx;
	}

	for (; i != x.word_end (); ++i)
		if (*i) return false;

	return x.back_word () == 0;
}

template <class Modules, class Vector1, class Vector2>
bool _equal<GF2, GenericModule<GF2>::Tag>::equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
						       VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j = y.begin ();
	typename Vector1::index_type idx;
	typename Vector1::word_type t;

	for (i = x.begin (); i != x.end (); ++i) {
		if (j == y.end () || *j < i->first)
			return false;

		idx = i->first << WordTraits<typename Vector1::word_type>::logof_size;
		t = Vector1::Endianness::e_0;

		for (; t != 0; t = Vector1::Endianness::shift_right (t, 1), ++idx) {
			if (i->second & t) {
				if (j == y.end () || *j != idx)
					return false;
				else
					++j;
			}
		}
	}

	return j == y.end ();
}

template <class Modules, class Vector1, class Vector2>
bool _equal<GF2, GenericModule<GF2>::Tag>::equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
						       VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01)
{
	typename Vector1::const_iterator i_1;
	typename Vector2::const_iterator i_2;

	if (x.size () != y.size ())
		return false;

	for (i_1 = x.begin (), i_2 = y.begin (); i_1 != x.end (); ++i_1, ++i_2)
		if (*i_1 != *i_2)
			return false;

	return true;
}

template <class Modules, class Vector>
bool _is_zero<GF2, GenericModule<GF2>::Tag>::is_zero_impl (const GF2 &F, Modules &M, const Vector &x, VectorRepresentationTypes::Dense01)
{
	typename Vector::const_word_iterator i;

	for (i = x.word_begin (); i != x.word_end (); ++i)
		if (*i) return false;

	if (x.back_word ())
		return false;

	return true;
}

template <class word, class Endianness>
inline int head_in_word (word v)
{
	// FIXME: This can be made faster...
	word mask = Endianness::e_0;
	int idx = 0;

	for (; mask != 0; mask = Endianness::shift_right (mask, 1), ++idx)
		if (v & mask)
			return idx;

	return -1;
}

template <class Modules, class reference, class Vector>
int _head<GF2, GenericModule<GF2>::Tag>::head_impl (const GF2 &F, Modules &M, reference &a, const Vector &x, VectorRepresentationTypes::Dense01)
{
	typename Vector::const_word_iterator i;
	size_t idx;

	for (i = x.word_begin (), idx = 0; i != x.word_end (); ++i, ++idx) {
		if (*i) {
			a = true;
			return head_in_word<typename Vector::word_type, typename Vector::Endianness> (*i) +
				(idx << WordTraits<typename Vector::word_type>::logof_size);
		}
	}

	if (x.back_word ()) {
		a = true;
		return head_in_word<typename Vector::word_type, typename Vector::Endianness> (x.back_word ()) +
			(idx << WordTraits<typename Vector::word_type>::logof_size);
	}

	return -1;
}

template <class Modules, class reference, class Vector>
int _head<GF2, GenericModule<GF2>::Tag>::head_impl (const GF2 &F, Modules &M, reference &a, const Vector &x, VectorRepresentationTypes::Hybrid01)
{
	if (x.empty ())
		return -1;
	else {
		a = true;
		return head_in_word<typename Vector::word_type, typename Vector::Endianness> (x.front ().second)
			+ (x.front ().first << WordTraits<typename Vector::word_type>::logof_size);
	}
}

template <class Modules, class Vector>
std::istream &_read<GF2, GenericModule<GF2>::Tag>::read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &x,
							      VectorRepresentationTypes::Dense01)
{
	typename Vector::iterator i;
	char c;

	do { is >> c ; } while (!std::isdigit (c));

	is.unget ();

	for (i = x.begin (); i != x.end (); ++i)
		is >> *i;

	return is;
}

template <class Modules, class Vector>
std::istream &_read<GF2, GenericModule<GF2>::Tag>::read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &x,
							      VectorRepresentationTypes::Sparse01)
{
	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type index_type;

	char c;
	index_type idx;

	do { is >> c ; } while (!std::isdigit (c));

	is.unget ();
	x.clear ();

	while (1) {
		is >> c;

		if (!std::isdigit (c) && c != ' ') break;
		is.unget ();
		is >> idx;
		x.push_back (idx);
	}

	return is;
}

template <class Modules, class Vector>
std::ostream &_write<GF2, GenericModule<GF2>::Tag>::write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &x,
								VectorRepresentationTypes::Dense01)
{
	os << "[ ";

	for (typename Vector::const_iterator i = x.begin (); i != x.end (); ++i)
		os << *i << ' ';

	os << ']';

	return os;
}

template <class Modules, class Vector>
std::ostream &_write<GF2, GenericModule<GF2>::Tag>::write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &x,
								VectorRepresentationTypes::Sparse01)
{
	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type index_type;

	typename Vector::const_iterator i;
	index_type idx = 0;

	os << "[ ";

	for (i = x.begin (); i != x.end (); ++i) {
		while (++idx <= *i)
			os << 0 << ' ';

		os << 1 << ' ';
	}
	os << ']';

	return os;
}

template <class Modules, class Vector>
std::ostream &_write<GF2, GenericModule<GF2>::Tag>::write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &x,
								VectorRepresentationTypes::Hybrid01)
{
	typename Vector::const_iterator i;
	typename Vector::index_type idx = 0;
	typename Vector::word_type mask;

	os << "[ ";

	for (i = x.begin (); i != x.end (); ++i) {
		while (++idx <= i->first << WordTraits<typename Vector::word_type>::logof_size)
			os << "0 ";

		for (mask = Vector::Endianness::e_0; mask != 0; mask = Vector::Endianness::shift_right (mask, 1)) {
			if (i->second & mask)
				os << "1 ";
			else
				os << "0 ";
			++idx;
		}
	}
	os << ']';

	return os;
}

} // namespace BLAS1

} // namespace LELA

#endif // __BLAS_LEVEL1_GF2_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
