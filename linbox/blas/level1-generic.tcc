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

namespace LinBox
{

namespace BLAS1 
{

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
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

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, y.size ()));
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		x.end () : std::lower_bound (x.begin (), x.end (), end_idx, VectorWrapper::CompareSparseEntries ());
		
	F.assign (res, F.zero ());

	for (; i != i_end; ++i)
		F.axpyin (res, i->second, y[i->first]);

	return res;
}

template <class Field, class Vector1, class Vector2>
typename Field::Element &dot_impl (const Field &F, GenericModule &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag)
{
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator j = (start_idx == 0) ? y.begin () : std::lower_bound (y.begin (), y.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		x.end () : std::lower_bound (x.begin (), x.end (), end_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		y.end () : std::lower_bound (y.begin (), y.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	F.assign (res, F.zero ());

	for (; i != i_end && j != j_end; ++i) {
		while (j != j_end && j->first < i->first) ++j;

		if (j != j_end && j->first == i->first)
			F.axpyin (res, i->second, j->second);
	}

	return res;
}

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag)
{
	typename Vector1::const_iterator i;
	int idx;

	y.clear ();

	for (i = x.begin (), idx = 0; i != x.end (); i++, idx++)
		if (!F.isZero (*i))
			y.push_back (typename Vector2::value_type (idx, *i));

	return y;
}

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, y.size ()));

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

template <class Field, class Vector1, class Vector2>
Vector2 &copy_impl (const Field &F, GenericModule &M, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag)
{
	y.clear ();
	y.insert (y.begin (), x.begin (), x.end ());
	return y;
}

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (y.size () == x.size ());

	typename Vector2::iterator i;
	typename Vector1::const_iterator j;

	for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
		F.axpyin (*i, a, *j);

	return y;
}

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::DenseVectorTag, VectorCategories::SparseVectorTag)
{
	typename Vector<Field>::Dense tmp (x.size ());
	_copy (F, M, x, tmp);
	_scal (F, M, a, tmp);
	_axpy (F, M, F.one (), y, tmp);
	_copy (F, M, tmp, y);
	return y;
}

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, y.size ()));

	typename Vector1::const_iterator j;

	for (j = x.begin (); j != x.end (); ++j)
		F.axpyin (y[j->first], a, j->second);

	return y;
}

template <class Field, class Vector1, class Vector2>
Vector2 &axpy_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y,
		    VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag)
{
	Vector2 tmp;

	typename Vector2::const_iterator i;
	typename Vector1::const_iterator j;
	typename Field::Element c;

	for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
		while (i != y.end () && i->first < j->first) {
			tmp.push_back (*i);
			i++;
		}

		if (i != y.end () && i->first == j->first) {
			F.axpy (c, a, j->second, i->second);
			i++;
		}
		else
			F.mul (c, a, j->second);

		if (!F.isZero (c))
			tmp.push_back (typename SparseVector<typename Field::Element>::value_type (j->first, c));
	}

	while (i != y.end ()) {
		tmp.push_back (*i);
		i++;
	}

	std::swap (y, tmp);

	return y;
}

template <class Field, class Vector>
Vector &scal_impl (const Field &F, GenericModule &M, const typename Field::Element &a, Vector &x, VectorCategories::DenseVectorTag)
{
	typename Vector::iterator i;

	for (i = x.begin (); i != x.end (); i++)
		F.mulin (*i, a);

	return x;
}

template <class Field, class Vector>
Vector &scal_impl (const Field &F, GenericModule &M, const typename Field::Element &a, Vector &x, VectorCategories::SparseVectorTag)
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

template <class Field, class Modules, class Iterator, class Vector>
Vector &permute_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::DenseVectorTag)
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

template <class Field, class Modules, class Iterator, class Vector>
Vector &permute_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorCategories::SparseVectorTag)
{
	typename Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		j->first = image (j->first, P_begin, P_end);

	std::sort (v.begin (), v.end (), VectorWrapper::CompareSparseEntries ());

	return v;
}

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag)
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

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseVectorTag, VectorCategories::DenseVectorTag)
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

template <class Field, class Vector1, class Vector2>
bool equal_impl (const Field &F, GenericModule &M, const Vector1 &x, const Vector2 &y,
		 VectorCategories::SparseVectorTag, VectorCategories::SparseVectorTag)
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

template <class Field, class Vector>
bool is_zero_impl (const Field &F, GenericModule &M, const Vector &x, VectorCategories::DenseVectorTag)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); i++)
		if (!F.isZero (*i))
			return false;

	return true;
}

template <class Field, class Vector>
bool is_zero_impl (const Field &F, GenericModule &M, const Vector &x, VectorCategories::SparseVectorTag)
{
	typename Vector::const_iterator i;

	for (i = x.begin (); i != x.end (); i++)
		if (!F.isZero (i->second))
			return false;

	return true;
}

template <class Field, class Vector>
int head_impl (const Field &F, GenericModule &M, typename Field::Element &a, const Vector &x, VectorCategories::DenseVectorTag)
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

template <class Field, class Vector>
int head_impl (const Field &F, GenericModule &M, typename Field::Element &a, const Vector &x, VectorCategories::SparseVectorTag)
{
	if (x.empty ())
		return -1;
	else {
		F.assign (a, x.front ().second);
		return x.front ().first;
	}
}

template <class Field, class Modules, class Vector>
std::istream &read_impl (const Field &F, Modules &M, std::istream &is, Vector &v, VectorCategories::DenseVectorTag)
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

template <class Field, class Modules, class Vector>
std::istream &read_impl (const Field &F, Modules &M, std::istream &is, Vector &v, VectorCategories::SparseVectorTag)
{
	typename Field::Element tmp;
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
			v.push_back (std::pair <size_t, typename Field::Element> (idx, tmp));

		is >> c;
		idx++;
	}

	return is;
}

template <class Field, class Modules, class Vector>
std::ostream &write_impl (const Field &F, Modules &M, std::ostream &os, const Vector &v, VectorCategories::DenseVectorTag)
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

template <class Field, class Modules, class Vector>
std::ostream &write_impl (const Field &F, Modules &M, std::ostream &os, const Vector &v, VectorCategories::SparseVectorTag)
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
