/* linbox/blas/level1-generic.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 1 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_GENERIC_TCC
#define __BLAS_LEVEL1_GENERIC_TCC

#include "linbox/blas/level1-generic.h"
#include "linbox/util/debug.h"
#include "linbox/blas/level1-ll.h"

namespace LinBox
{

namespace BLAS1 
{

template <class Ring>
template <class Modules, class Vector1, class Vector2>
typename Ring::Element &_dot<Ring, GenericModule::Tag>::dot_impl (const Ring &F, Modules &M, typename Ring::Element &res, const Vector1 &x, const Vector2 &y,
								  size_t start_idx, size_t end_idx,
								  VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	linbox_check (x.size () == y.size ());
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i, i_end = x.begin () + std::min (x.size (), end_idx);
	typename Vector2::const_iterator j;

	F.assign (res, F.zero ());

	for (i = x.begin () + start_idx, j = y.begin () + start_idx; i != i_end; ++i, ++j)
		F.axpyin (res, *i, *j);

	return res;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
typename Ring::Element &_dot<Ring, GenericModule::Tag>::dot_impl (const Ring &F, Modules &M, typename Ring::Element &res, const Vector1 &x, const Vector2 &y,
								  size_t start_idx, size_t end_idx,
								  VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, y.size ()));
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorUtils::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		x.end () : std::lower_bound (x.begin (), x.end (), end_idx, VectorUtils::CompareSparseEntries ());
		
	F.assign (res, F.zero ());

	for (; i != i_end; ++i)
		F.axpyin (res, i->second, y[i->first]);

	return res;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
typename Ring::Element &_dot<Ring, GenericModule::Tag>::dot_impl (const Ring &F, Modules &M, typename Ring::Element &res, const Vector1 &x, const Vector2 &y,
								  size_t start_idx, size_t end_idx,
								  VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorUtils::CompareSparseEntries ());
	typename Vector2::const_iterator j = (start_idx == 0) ? y.begin () : std::lower_bound (y.begin (), y.end (), start_idx, VectorUtils::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		x.end () : std::lower_bound (x.begin (), x.end (), end_idx, VectorUtils::CompareSparseEntries ());
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		y.end () : std::lower_bound (y.begin (), y.end (), end_idx, VectorUtils::CompareSparseEntries ());

	F.assign (res, F.zero ());

	for (; i != i_end && j != j_end; ++i) {
		while (j != j_end && j->first < i->first) ++j;

		if (j != j_end && j->first == i->first)
			F.axpyin (res, i->second, j->second);
	}

	return res;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<Ring, GenericModule::Tag>::copy_impl (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
						     VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i;
	int idx;

	y.clear ();

	for (i = x.begin (), idx = 0; i != x.end (); i++, idx++)
		if (!F.isZero (*i))
			y.push_back (typename Vector2::value_type (idx, *i));

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_copy<Ring, GenericModule::Tag>::copy_impl (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
						     VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, y.size ()));

	typename Vector1::const_iterator j;
	typename Vector2::iterator i;
	size_t idx;

	for (i = y.begin (), j = x.begin (), idx = 0; j != x.end (); i++, j++, idx++) {
		while (idx < j->first) {
			F.assign (*i, F.zero ());
			i++; idx++;
		}

		*i = j->second;
	}

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, GenericModule::Tag>::axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
						     VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	linbox_check (y.size () == x.size ());

	typename Vector2::iterator i;
	typename Vector1::const_iterator j;

	for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
		F.axpyin (*i, a, *j);

	return y;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, GenericModule::Tag>::axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
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
Vector2 &_axpy<Ring, GenericModule::Tag>::axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
						     VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, y.size ()));

	typename Vector1::const_iterator j;

	for (j = x.begin (); j != x.end (); ++j)
		F.axpyin (y[j->first], a, j->second);

	return y;
}

template <class Ring, class Modules, class Vector>
void fast_copy (const Ring &F, Modules &M, SparseVector<typename Ring::Element, std::vector<typename Vector::index_type>, std::vector<typename Vector::element_type> > &v, Vector &w)
	{ _copy<Ring, typename Modules::Tag>::op (F, M, v, w); }

template <class Ring, class Modules, class index_type>
void fast_copy (const Ring &F, Modules &M,
		SparseVector<typename Ring::Element, std::vector<index_type>, std::vector<typename Ring::Element> > &v,
		SparseVector<typename Ring::Element, std::vector<index_type>, std::vector<typename Ring::Element> > &w)
	{ std::swap (v, w); }

template <class Ring>
template <class Modules, class Vector1, class Vector2>
Vector2 &_axpy<Ring, GenericModule::Tag>::axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
						     VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	SparseVector<typename Ring::Element, std::vector<typename Vector2::index_type>, std::vector<typename Vector2::element_type> > tmp;

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

		if (!F.isZero (c))
			tmp.push_back (typename Vector2::value_type (i->first, c));
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
Vector &_scal<Ring, GenericModule::Tag>::scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x, VectorRepresentationTypes::Dense)
{
	typename Vector::iterator i;

	for (i = x.begin (); i != x.end (); i++)
		F.mulin (*i, a);

	return x;
}

template <class Ring>
template <class Modules, class Vector>
Vector &_scal<Ring, GenericModule::Tag>::scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x, VectorRepresentationTypes::Sparse)
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
Vector &_permute<Ring, GenericModule::Tag>::permute_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Dense)
{
	for (; P_begin != P_end; ++P_begin)
		std::swap (v[P_begin->first], v[P_begin->second]);

	return v;
}

template <class Iterator>
typename Iterator::value_type::first_type image (typename Iterator::value_type::first_type x, Iterator P_begin, Iterator P_end)
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
Vector &_permute<Ring, GenericModule::Tag>::permute_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Sparse)
{
	typename Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		j->first = image (j->first, P_begin, P_end);

	std::sort (v.begin (), v.end (), VectorUtils::CompareSparseEntries ());

	return v;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
bool _equal<Ring, GenericModule::Tag>::equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
						   VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	if (x.size () != y.size ())
		return false;

	for (i = x.begin (), j = y.begin (); i != x.end (); i++, j++)
		if (!F.areEqual (*i, *j))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
bool _equal<Ring, GenericModule::Tag>::equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
						   VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;
	size_t idx;

	for (i = x.begin (), j = y.begin (), idx = 0; i != x.end () && j != y.end (); j++, idx++) {
		if (i->first == idx) {
			if (!F.areEqual (i->second, *j))
				return false;
			i++;
		}
		else if (!F.isZero (*j))
			return false;
	}

	return true;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2>
bool _equal<Ring, GenericModule::Tag>::equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
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
bool _is_zero<Ring, GenericModule::Tag>::is_zero_impl (const Ring &F, Modules &M, const Vector &x, VectorRepresentationTypes::Dense)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); i++)
		if (!F.isZero (*i))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Vector>
bool _is_zero<Ring, GenericModule::Tag>::is_zero_impl (const Ring &F, Modules &M, const Vector &x, VectorRepresentationTypes::Sparse)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); i++)
		if (!F.isZero (i->second))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Vector>
int _head<Ring, GenericModule::Tag>::head_impl (const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x, VectorRepresentationTypes::Dense)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end () && F.isZero (*i); ++i);

	if (i == x.end ())
		return -1;
	else {
		F.assign (a, *i);
		return i - x.begin ();
	}
}

template <class Ring>
template <class Modules, class Vector>
int _head<Ring, GenericModule::Tag>::head_impl (const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x, VectorRepresentationTypes::Sparse)
{
	if (x.empty ())
		return -1;
	else {
		F.assign (a, x.front ().second);
		return x.front ().first;
	}
}

template <class Ring>
template <class Modules, class Vector>
std::istream &_read<Ring, GenericModule::Tag>::read_impl (const Ring &F, Modules &M, std::istream &is, Vector &v, VectorRepresentationTypes::Dense)
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
std::istream &_read<Ring, GenericModule::Tag>::read_impl (const Ring &F, Modules &M, std::istream &is, Vector &v, VectorRepresentationTypes::Sparse)
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
std::ostream &_write<Ring, GenericModule::Tag>::write_impl (const Ring &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Dense)
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
std::ostream &_write<Ring, GenericModule::Tag>::write_impl (const Ring &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Sparse)
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

} // namespace LinBox

#endif // __BLAS_LEVEL1_GENERIC_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
