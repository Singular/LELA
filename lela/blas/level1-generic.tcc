/* lela/blas/level1-generic.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 1 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_GENERIC_TCC
#define __BLAS_LEVEL1_GENERIC_TCC

#include "lela/blas/level1-generic.h"
#include "lela/util/debug.h"
#include "lela/blas/level1-ll.h"

namespace LELA
{

namespace BLAS1 
{

template <class Ring>
template <class Modules, class T, class Vector1, class Vector2>
T &_dot<Ring, typename GenericModule<Ring>::Tag>::dot_impl
	(const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	lela_check (x.size () == y.size ());

	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	F.copy (res, F.zero ());

	for (i = x.begin (), j = y.begin (); i != x.end (); ++i, ++j)
		F.axpyin (res, *i, *j);

	return res;
}

template <class Ring>
template <class Modules, class T, class Vector1, class Vector2>
T &_dot<Ring, typename GenericModule<Ring>::Tag>::dot_impl
	(const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Ring> (x, y.size ()));

	typename Vector1::const_iterator i = x.begin ();

	F.copy (res, F.zero ());

	for (; i != x.end (); ++i)
		F.axpyin (res, i->second, y[i->first]);

	return res;
}

template <class Ring>
template <class Modules, class T, class Vector1, class Vector2>
T &_dot<Ring, typename GenericModule<Ring>::Tag>::dot_impl
	(const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	F.copy (res, F.zero ());

	for (; i != x.end () && j != y.end (); ++i) {
		while (j != y.end () && j->first < i->first) ++j;

		if (j != y.end () && j->first == i->first)
			F.axpyin (res, i->second, j->second);
	}

	return res;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<Ring, typename GenericModule<Ring>::Tag>::copy_impl
	(const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	lela_check (x.size () == y.size ());

	typename Vector1::const_iterator i;
	typename Vector2::iterator j;

	for (i = x.begin (), j = y.begin (); i != x.end (); ++i, ++j)
		F.copy (*j, *i);

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<Ring, typename GenericModule<Ring>::Tag>::copy_impl
	(const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i;
	size_t idx;

	y.clear ();

	for (i = x.begin (), idx = 0; i != x.end (); ++i, ++idx) {
		if (!F.isZero (*i)) {
			y.push_back (typename Vector2::value_type (idx, typename Ring::Element ()));
			F.copy (y.back ().second, *i);
		}
	}

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<Ring, typename GenericModule<Ring>::Tag>::copy_impl
	(const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i;

	y.clear ();

	for (i = x.begin (); i != x.end (); ++i) {
		y.push_back (typename Vector2::value_type (i->first, typename Ring::Element ()));
		F.copy (y.back ().second, i->second);
	}

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<Ring, typename GenericModule<Ring>::Tag>::copy_impl
	(const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Ring> (x, y.size ()));

	typename Vector1::const_iterator i;

	_scal<Ring, typename Modules::Tag>::op (F, M, F.zero (), y);

	for (i = x.begin (); i != x.end (); i++)
		F.copy (y[i->first], i->second);

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, typename GenericModule<Ring>::Tag>::axpy_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	lela_check (y.size () == x.size ());

	typename Vector2::iterator i;
	typename Vector1::const_iterator j;

	if (F.isOne (a))
		for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
			F.addin (*i, *j);
	else if (F.areEqual (a, F.minusOne ()))
		for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
			F.subin (*i, *j);
	else
		for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
			F.axpyin (*i, a, *j);

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, typename GenericModule<Ring>::Tag>::axpy_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
{
	typename Vector<Ring>::Dense tmp (x.size ());
	_copy<Ring, typename Modules::Tag>::op (F, M, x, tmp);
	_scal<Ring, typename Modules::Tag>::op (F, M, a, tmp);
	_axpy<Ring, typename Modules::Tag>::op (F, M, F.one (), y, tmp);
	_copy<Ring, typename Modules::Tag>::op (F, M, tmp, y);
	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, typename GenericModule<Ring>::Tag>::axpy_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Ring> (x, y.size ()));

	typename Vector1::const_iterator j;

	for (j = x.begin (); j != x.end (); ++j)
		F.axpyin (y[j->first], a, j->second);

	return y;
}

template <class Ring, class Modules, class Vector>
void fast_copy (const Ring &F, Modules &M, SparseVector<typename Ring::Element, std::vector<typename Vector::value_type::first_type>, std::vector<typename Vector::value_type::second_type> > &v, Vector &w)
	{ w.assign (v.begin (), v.end ()); }

template <class Ring, class Modules, class index_type>
void fast_copy (const Ring &F, Modules &M,
		SparseVector<typename Ring::Element, std::vector<index_type>, std::vector<typename Ring::Element> > &v,
		SparseVector<typename Ring::Element, std::vector<index_type>, std::vector<typename Ring::Element> > &w)
	{ std::swap (v, w); }

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, typename GenericModule<Ring>::Tag>::axpy_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	SparseVector<typename Ring::Element, std::vector<typename Vector2::value_type::first_type>, std::vector<typename Vector2::value_type::second_type> > tmp;

	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;
	typename Ring::Element c;

	for (i = x.begin (), j = y.begin (); i != x.end (); i++) {
		while (j != y.end () && j->first < i->first) {
			tmp.push_back (*j);
			j++;
		}

		if (j != y.end () && i->first == j->first) {
			F.axpy (c, a, i->second, j->second);
			j++;
		}
		else
			F.mul (c, a, i->second);

		if (!F.isZero (c)) {
			tmp.push_back (typename Vector2::value_type (i->first, typename Ring::Element ()));
			F.copy (tmp.back ().second, c);
		}
	}

	while (j != y.end ()) {
		tmp.push_back (*j);
		j++;
	}

	fast_copy (F, M, tmp, y);

	return y;
}

template <class Ring>
template <class Modules, class Vector>
Vector &_scal<Ring, typename GenericModule<Ring>::Tag>::scal_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x, VectorRepresentationTypes::Dense)
{
	typename Vector::iterator i;

	if (F.isZero (a))
		for (i = x.begin (); i != x.end (); i++)
			F.copy (*i, F.zero ());
	else if (F.areEqual (a, F.minusOne ()))
		for (i = x.begin (); i != x.end (); i++)
			F.negin (*i);
	else
		for (i = x.begin (); i != x.end (); i++)
			F.mulin (*i, a);

	return x;
}

template <class Ring>
template <class Modules, class Vector>
Vector &_scal<Ring, typename GenericModule<Ring>::Tag>::scal_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x, VectorRepresentationTypes::Sparse)
{
	typename Vector::iterator i;

	if (F.isZero (a)) {
		x.clear ();
		return x;
	}

	for (i = x.begin (); i != x.end (); i++)
		F.mulin (i->second, a);

	return x;
}

template <class Ring>
template <class Modules, class Iterator, class Vector>
Vector &_permute<Ring, typename GenericModule<Ring>::Tag>::permute_impl
	(const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Dense)
{
	for (; P_begin != P_end; ++P_begin)
		std::swap (v[P_begin->first], v[P_begin->second]);

	return v;
}

template <class Iterator>
typename std::iterator_traits<Iterator>::value_type::first_type image
	(typename std::iterator_traits<Iterator>::value_type::first_type x, Iterator P_begin, Iterator P_end)
{
	for (; P_begin != P_end; ++P_begin) {
		if (P_begin->first == x)
			x = P_begin->second;
		else if (P_begin->second == x)
			x = P_begin->first;
	}

	return x;
}

template <class Ring>
template <class Modules, class Iterator, class Vector>
Vector &_permute<Ring, typename GenericModule<Ring>::Tag>::permute_impl
	(const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Sparse)
{
	typename Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		j->first = image (j->first, P_begin, P_end);

	std::sort (v.begin (), v.end (), VectorUtils::CompareSparseEntries ());

	return v;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
bool _equal<Ring, typename GenericModule<Ring>::Tag>::equal_impl
	(const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	if (x.size () != y.size ())
		return false;

	for (i = x.begin (), j = y.begin (); i != x.end (); ++i, ++j)
		if (!F.areEqual (*i, *j))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
bool _equal<Ring, typename GenericModule<Ring>::Tag>::equal_impl
	(const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;
	size_t idx;

	for (i = x.begin (), j = y.begin (), idx = 0; j != y.end (); ++j, ++idx) {
		if (i != x.end () && i->first == idx) {
			if (!F.areEqual (i->second, *j))
				return false;
			i++;
		}
		else if (!F.isZero (*j))
			return false;
	}

	// If i did not reach the end of x, then the input was invalid
	lela_check (i == x.end ());

	return true;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
bool _equal<Ring, typename GenericModule<Ring>::Tag>::equal_impl
	(const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	if (x.size () != y.size ())
		return false;

	for (i = x.begin (), j = y.begin (); i != x.end (); ++i, ++j) {
		if (i->first != j->first || !F.areEqual (i->second, j->second))
			return false;
	}

	return true;
}

template <class Ring>
template <class Modules, class Vector>
bool _is_zero<Ring, typename GenericModule<Ring>::Tag>::is_zero_impl
	(const Ring &F, Modules &M, const Vector &x, VectorRepresentationTypes::Dense)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); i++)
		if (!F.isZero (*i))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Vector>
bool _is_zero<Ring, typename GenericModule<Ring>::Tag>::is_zero_impl
	(const Ring &F, Modules &M, const Vector &x, VectorRepresentationTypes::Sparse)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); i++)
		if (!F.isZero (i->second))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Vector>
int _head<Ring, typename GenericModule<Ring>::Tag>::head_impl
	(const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x, VectorRepresentationTypes::Dense)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end () && F.isZero (*i); ++i);

	if (i == x.end ())
		return -1;
	else {
		F.copy (a, *i);
		return i - x.begin ();
	}
}

template <class Ring>
template <class Modules, class Vector>
int _head<Ring, typename GenericModule<Ring>::Tag>::head_impl
	(const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x, VectorRepresentationTypes::Sparse)
{
	if (x.empty ())
		return -1;
	else {
		F.copy (a, x.front ().second);
		return x.front ().first;
	}
}

template <class Ring>
template <class Modules, class Vector>
std::istream &_read<Ring, typename GenericModule<Ring>::Tag>::read_impl
	(const Ring &F, Modules &M, std::istream &is, Vector &v, VectorRepresentationTypes::Dense)
{
	typename Vector::iterator i;
	char c;
	bool seekrightbracket = false;

	i = v.begin ();
	do is >> c; while (is && isspace (c));

	if (c == '[')
		seekrightbracket = true;
	else
		is.unget ();

	while (i != v.end() && is) {
		do is >> c; while (!isdigit (c) && c != '-');
		is.unget ();
		F.read (is, *i++);
	}

	if (seekrightbracket)
		do is >> c; while (is && c != ']');

	return is;
}

template <class Ring>
template <class Modules, class Vector>
std::istream &_read<Ring, typename GenericModule<Ring>::Tag>::read_impl
	(const Ring &F, Modules &M, std::istream &is, Vector &v, VectorRepresentationTypes::Sparse)
{
	typename Ring::Element tmp;
	char c;
	int idx;

	do is >> c; while (is && isspace (c));

	if (isdigit (c))
		is.unget ();

	c = ','; v.clear (); idx = 0;

	while (is && c == ',') {
		do is >> c; while (is && isspace (c));
		is.unget ();
		F.read (is, tmp);

		if (!F.isZero (tmp))
			v.push_back (std::pair <size_t, typename Ring::Element> (idx, tmp));

		is >> c;
		idx++;
	}

	return is;
}

template <class Ring>
template <class Modules, class Vector>
std::ostream &_write<Ring, typename GenericModule<Ring>::Tag>::write_impl
	(const Ring &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Dense)
{
	typename Vector::const_iterator i;

	os << '[';

	for (i = v.begin (); i != v.end ();) {
		F.write (os, *i);

		if (++i != v.end ())
			os << ", ";
	}

	os << ']';

	return os;
}

template <class Ring>
template <class Modules, class Vector>
std::ostream &_write<Ring, typename GenericModule<Ring>::Tag>::write_impl
	(const Ring &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Sparse)
{
	typename Vector::const_iterator i;
	size_t idx;

	os << '[';

	for (i = v.begin (), idx = 0; i != v.end ();) {
		while (idx++ < i->first)
			os << "0, ";

		F.write (os, i->second);

		if (++i != v.end ())
			os << ", ";
	}

	os << ']';

	return os;
}

template <class Iterator>
std::ostream &write_permutation (std::ostream &os, Iterator P_begin, Iterator P_end)
{
	for (; P_begin != P_end; ++P_begin)
		os << "(" << P_begin->first << " " << P_begin->second << ")";

	return os;
}

} // namespace BLAS1

} // namespace LELA

#endif // __BLAS_LEVEL1_GENERIC_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
