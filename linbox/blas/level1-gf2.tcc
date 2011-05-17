/* linbox/blas/level1-gf2.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for GF2
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_GF2_TCC
#define __BLAS_LEVEL1_GF2_TCC

#include <algorithm>
#include <iostream>

#include "linbox/blas/level1-gf2.h"
#include "linbox/blas/level1-generic.h"

namespace LinBox
{

namespace BLAS1 
{

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	linbox_check (x.size () == y.size ());
	linbox_check (start_idx <= end_idx);
	linbox_check (start_idx <= x.size ());

	if (x.empty ())
		return res = false;

	typename Vector1::word_type t = 0;
	typename Vector1::const_word_iterator i = x.wordBegin () + (start_idx >> WordTraits<typename Vector1::word_type>::logof_size);
	typename Vector2::const_word_iterator j = y.wordBegin () + (start_idx >> WordTraits<typename Vector1::word_type>::logof_size);

	typename Vector1::const_word_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		x.wordEnd () : x.wordBegin () + ((end_idx + WordTraits<typename Vector1::word_type>::bits - 1) >> WordTraits<typename Vector1::word_type>::logof_size);

	if (i == i_end)
		return res = false;
	else if (i == i_end - 1)
		t = *i & *j &
			Vector1::Endianness::mask_right (start_idx & WordTraits<typename Vector1::word_type>::pos_mask) &
			Vector1::Endianness::mask_left (end_idx & WordTraits<typename Vector1::word_type>::pos_mask);
	else {
		if ((start_idx & WordTraits<typename Vector1::word_type>::pos_mask) != 0)
			t = *i++ & *j++ & Vector1::Endianness::mask_right (start_idx & WordTraits<typename Vector1::word_type>::pos_mask);

		while (i != i_end - 1)
			t ^= *i++ & *j++;
        
		t ^= *i & *j & Vector1::Endianness::mask_left (end_idx & WordTraits<typename Vector1::word_type>::pos_mask);
	}

        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
{
	linbox_check (VectorWrapper::hasDim<GF2> (y, x.size ()));
	linbox_check (start_idx <= end_idx);
	linbox_check (start_idx <= x.size ());

	typename Vector2::const_iterator i = (start_idx == 0) ? y.begin () : std::lower_bound (y.begin (), y.end (), start_idx);

	typename Vector2::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		y.end () : std::lower_bound (y.begin (), y.end (), end_idx);

	res = false;

	for (; i != i_end; ++i)
		if (x[*i]) res = !res;

	return res;
}

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
{
	linbox_check (VectorWrapper::hasDim<GF2> (y, x.size ()));
	linbox_check (start_idx <= end_idx);
	linbox_check (start_idx <= x.size ());

	typename Vector1::word_type t = 0;
	typename Vector1::const_word_iterator i = x.wordBegin ();
	typename Vector2::const_iterator j = (start_idx == 0) ?
		y.begin () : std::lower_bound (y.begin (), y.end (), start_idx >> WordTraits<typename Vector2::word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

	typename Vector2::index_type search_idx = (end_idx + WordTraits<typename Vector1::word_type>::bits - 1) >> WordTraits<typename Vector2::word_type>::logof_size;

	typename Vector2::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		y.end () : std::lower_bound (y.begin (), y.end (), search_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector2::const_iterator j_stop;

	if (j == j_end)
		return res = false;
	else if (j == j_end - 1) {
		t = *(i + j->first) & j->second;

		if (j->first == start_idx >> WordTraits<typename Vector1::word_type>::logof_size)
			t &= Vector1::Endianness::mask_right (start_idx & WordTraits<typename Vector1::word_type>::pos_mask);

		if (j->first == end_idx >> WordTraits<typename Vector1::word_type>::logof_size)
			t &= Vector1::Endianness::mask_left (end_idx & WordTraits<typename Vector1::word_type>::pos_mask);
	} else {
		if (end_idx == static_cast<size_t> (-1))
			j_stop = j_end;
		else
			j_stop = j_end - 1;

		if (j->first == start_idx >> WordTraits<typename Vector1::word_type>::logof_size) {
			t = *(i + j->first) & j->second & Vector1::Endianness::mask_right (start_idx & WordTraits<typename Vector1::word_type>::pos_mask);
			++j;
		}

		for (; j != j_stop; ++j)
			t ^= *(i + j->first) & j->second;

		if (end_idx != static_cast<size_t> (-1) && j->first == end_idx >> WordTraits<typename Vector1::word_type>::logof_size)
			t ^= *(i + j->first) & j->second & Vector1::Endianness::mask_left (end_idx & WordTraits<typename Vector1::word_type>::pos_mask);
	}
        
        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class reference, class Vector1, class Vector2>
reference &dot_impl (const GF2 &F, GenericModule &M, reference &res, const Vector1 &x, const Vector2 &y,
		     size_t start_idx, size_t end_idx,
		     VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
{
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx);
	typename Vector2::const_iterator j = (start_idx == 0) ? y.begin () : std::lower_bound (y.begin (), y.end (), start_idx);

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		x.end () : std::lower_bound (x.begin (), x.end (), end_idx);
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		y.end () : std::lower_bound (y.begin (), y.end (), end_idx);

	res = false;

	while (i != i_end || j != j_end) {
		while (i != i_end && (j == j_end || *i < *j)) ++i;
		while (j != j_end && (i == i_end || *j < *i)) ++j;
		if (i != i_end && j != j_end && *i == *j) { res = !res; ++i; ++j; }
	}

	return res;
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
{
	typedef typename std::iterator_traits<typename Vector2::const_iterator>::value_type index_type;

	typename Vector1::const_iterator i;
	index_type idx = 0;

	y.clear ();

	for (i = x.begin (); i != x.end (); ++i, ++idx)
		if (*i) y.push_back (idx);

	return y;
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
{
	typename Vector1::const_word_iterator i;
	typename Vector2::index_type idx = 0;

	y.clear ();

	for (i = x.wordBegin (); i != x.wordEnd (); ++i, ++idx)
		if (*i)
			y.push_back (typename Vector1::value_type (idx, *i));

	return y;
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	linbox_check (VectorWrapper::hasDim<GF2> (x, y.size ()));

	typename Vector1::const_iterator i;

	std::fill (y.wordBegin (), y.wordEnd (), 0);

	for (i = x.begin (); i != x.end (); ++i)
        	y[*i] = true;

	return y;
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
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

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	typename Vector1::const_iterator i;

	std::fill (y.wordBegin (), y.wordEnd (), 0);

	for (i = x.begin (); i != x.end (); ++i)
        	*(y.wordBegin () + i->first) = i->second;

	return y;
}

template <class Vector1, class Vector2>
Vector2 &copy_impl (const GF2 &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
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

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	linbox_check (y.size () == x.size ());

	if (a) {
		typename Vector2::word_iterator i = y.wordBegin ();
		typename Vector1::const_word_iterator j = x.wordBegin ();

		for (; i != y.wordEnd (); ++i, ++j)
			*i ^= *j;
	}

	return y;
}

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	linbox_check (VectorWrapper::hasDim<GF2> (x, y.size ()));

	if (a) {
		typename Vector1::const_iterator i;

		for (i = x.begin (); i != x.end (); ++i)
			y[*i] = !y[*i];
	}

	return y;
}

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
{
	if (a) {
		Vector2 res;

		typename Vector2::const_iterator i = y.begin ();
		typename Vector1::const_iterator j = x.begin ();

		while (i != y.end () || j != x.end ()) {
			while (i != y.end () && (j == x.end () || *i < *j)) { res.push_back (*i); ++i; }
			while (j != x.end () && (i == y.end () || *j < *i)) { res.push_back (*j); ++j; }
			if (i != y.end () && j != x.end () && *i == *j) { ++i; ++j; }
		}

		std::swap (y, res);
	}

	return y;
}

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	linbox_check (VectorWrapper::hasDim<GF2> (x, y.size ()));

	if (a) {
		typename Vector1::const_iterator i;
		typename Vector2::word_iterator j = y.wordBegin ();

		for (i = x.begin (); i != x.end (); ++i)
			*(j + i->first) ^= i->second;
	}

	return y;
}

template <class Vector1, class Vector2>
Vector2 &axpy_impl (const GF2 &F, GenericModule &M, bool a, const Vector1 &x, Vector2 &y,
		    VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
{
	if (a) {
		Vector2 res;

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

		std::swap (y, res);
	}

	return y;
}

template <class Modules, class Iterator, class Vector>
Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::SparseZeroOneVectorTag)
{
	typename Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		*j = image (*j, P_begin, P_end);

	std::sort (v.begin (), v.end ());

	return v;
}

template <class Modules, class Iterator, class Vector>
Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::HybridZeroOneVectorTag)
{
	::LinBox::Vector<GF2>::Sparse w;

	_copy (F, M, v, w);
	_permute (F, M, P_begin, P_end, w);
	_copy (F, M, w, v);

	return v;
}

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseZeroOneVectorTag, VectorCategories::DenseZeroOneVectorTag)
{
	typename Vector1::const_word_iterator i = x.wordBegin ();
	typename Vector2::const_word_iterator j = y.wordBegin ();

	for (; j != y.wordEnd (); ++j, ++i)
		if (*i != *j) return false;

	return true;
}

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
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

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
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

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
{
	if (x.empty () && !y.empty ())
		return false;

	typename Vector1::const_word_iterator i = x.wordBegin ();
	typename Vector2::const_iterator j = y.begin ();
	typename Vector2::index_type idx = 0;

	for (; j != y.end (); ++j) {
		while (idx < j->first) {
			if (*i) return false;
			++idx;
			++i;

			if (i == x.wordEnd ())
				return false;
		}

		if (*i++ != j->second)
			return false;

		++idx;
	}

	for (; i != x.wordEnd (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::HybridZeroOneVectorTag, VectorCategories::SparseZeroOneVectorTag)
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

template <class Vector1, class Vector2>
bool equal_impl (const GF2 &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::HybridZeroOneVectorTag, VectorCategories::HybridZeroOneVectorTag)
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

template <class Vector>
bool is_zero_impl (const GF2 &F, GenericModule &M, const Vector &x, VectorCategories::DenseZeroOneVectorTag)
{
	typename Vector::const_word_iterator i;

	for (i = x.wordBegin (); i != x.wordEnd (); ++i)
		if (*i) return false;

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

template <class reference, class Vector>
int head_impl (const GF2 &F, GenericModule &M, reference &a, const Vector &x, VectorCategories::DenseZeroOneVectorTag)
{
	// FIXME: This can be made faster...
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); ++i) {
		if (*i) {
			a = true;
			return i - x.begin ();
		}
	}

	return -1;
}

template <class reference, class Vector>
int head_impl (const GF2 &F, GenericModule &M, reference &a, const Vector &x, VectorCategories::HybridZeroOneVectorTag)
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
std::istream &read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &x, VectorCategories::DenseZeroOneVectorTag)
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
std::istream &read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &x, VectorCategories::SparseZeroOneVectorTag)
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
std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &x, VectorCategories::DenseZeroOneVectorTag)
{
	os << "[ ";

	for (typename Vector::const_iterator i = x.begin (); i != x.end (); ++i)
		os << *i << ' ';

	os << ']';

	return os;
}

template <class Modules, class Vector>
std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &x, VectorCategories::SparseZeroOneVectorTag)
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
std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &x, VectorCategories::HybridZeroOneVectorTag)
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

} // namespace LinBox

#endif // __BLAS_LEVEL1_GF2_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
