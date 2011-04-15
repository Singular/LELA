/* linbox/vector/vector-domain.inl
 * Copyright 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-07-24 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added support for the new SparseParallel vector type; this involves quite a
 * few new specializations.
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>
 *
 * Added the modifications for categories and vector traits that were designed 
 * at the Rootbeer meeting. Added parametrization of VectorTags by VectorTraits.
 * 
 * ------------------------------------
 * 2002-06-04 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Updated function definitions according to the new policy set in
 * vector-domain.h  This means the functions each take a parameter tag (or, for
 * dot product, tag1 and tag2) that allows specialization by vector type.
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_vector_domain_INL
#define __LINBOX_field_vector_domain_INL

#include "linbox/linbox-config.h"

#include <iostream>
#include <cctype>

#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"

namespace LinBox
{

template <class Field>
template <class Vector>
std::ostream &VectorDomain<Field>::writeSpecialized (std::ostream &os, const Vector &x,
						     VectorCategories::DenseVectorTag) const
{
	typename Vector::const_iterator i;

	os << '[';

	for (i = x.begin (); i != x.end ();) {
		VectorDomainBase<Field>::_F.write (os, *i);

		if (++i != x.end ())
			os << ", ";
	}

	os << ']';

	return os;
}

template <class Field>
template <class Vector>
std::ostream &VectorDomain<Field>::writeSpecialized (std::ostream &os, const Vector &x,
						     VectorCategories::SparseVectorTag) const
{
	typename Vector::const_iterator i;
	size_t idx;

	os << '[';

	for (i = x.begin (), idx = 0; i != x.end ();) {
		while (idx++ < i->first)
			os << "0, ";

		VectorDomainBase<Field>::_F.write (os, i->second);

		if (++i != x.end ())
			os << ", ";
	}

	os << ']';

	return os;
}

template <class Field>
template <class Vector>
std::istream &VectorDomain<Field>::readSpecialized (std::istream &is, Vector &x,
						    VectorCategories::DenseVectorTag) const
{
	typename Vector::iterator i;
	char c;
	bool seekrightbracket = false;

	i = x.begin ();
	do is >> c; while (is && isspace (c));
	if (c == '[') seekrightbracket = true;
	else is.unget ();

	while (i != x.end() && is) {
		do is >> c; while (!isdigit(c) && c != '-');
		is.unget ();
		VectorDomainBase<Field>::_F.read (is, *i++);
		//std::cerr << std::endl << "just read this: ";
		//VectorDomainBase<Field>::_F.write(cerr, *(i-1)) << " at index " << (i-x.begin());
	}
	if (seekrightbracket) do is >> c; while (is && c != ']');

	return is;
}

template <class Field>
template <class Vector>
std::istream &VectorDomain<Field>::readSpecialized (std::istream &is, Vector &x,
						    VectorCategories::SparseVectorTag) const
{
	typename Field::Element tmp;
	char c;
	int idx;

	do is >> c; while (is && isspace (c));

	if (isdigit (c))
		is.unget ();

	c = ','; x.clear (); idx = 0;

	while (is && c == ',') {
		do is >> c; while (is && isspace (c));
		is.unget ();
		VectorDomainBase<Field>::_F.read (is, tmp);
		if (!VectorDomainBase<Field>::_F.isZero (tmp))
			x.push_back (std::pair <size_t, typename Field::Element> (idx, tmp));
		is >> c;
		idx++;
	}

	return is;
}

template <class Field>
template <class Vector1, class Vector2>
bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					       VectorCategories::DenseVectorTag,
					       VectorCategories::DenseVectorTag) const
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	if (v1.size () != v2.size ()) return false;

	for (i = v1.begin (), j = v2.begin (); i != v1.end (); i++, j++)
		if (!VectorDomainBase<Field>::_F.areEqual (*i, *j))
			return false;

	return true;
}

template <class Field>
template <class Vector1, class Vector2>
bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					       VectorCategories::SparseVectorTag,
					       VectorCategories::DenseVectorTag) const
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;
	size_t idx;

	for (i = v1.begin (), j = v2.begin (), idx = 0; i != v1.end () && j != v2.end (); j++, idx++) {
		if (i->first == idx) {
			if (!VectorDomainBase<Field>::_F.areEqual (i->second, *j))
				return false;
			i++;
		}
		else if (!VectorDomainBase<Field>::_F.isZero (*j))
			return false;
	}

	return true;
}

template <class Field>
template <class Vector1, class Vector2>
bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					       VectorCategories::SparseVectorTag,
					       VectorCategories::SparseVectorTag) const
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	for (i = v1.begin (), j = v2.begin (); i != v1.end () || j != v2.end ();) {
		while (i != v1.end () && (j == v2.end () || i->first < j->first)) {
			if (!VectorDomainBase<Field>::_F.isZero (i->second))
				return false;
			i++;
		}

		while (j != v2.end () && (i == v1.end () || j->first < i->first)) {
			if (!VectorDomainBase<Field>::_F.isZero (j->second))
				return false;
			j++;
		}

		if (i != v1.end () && j != v2.end () && i->first == j->first) {
			if (!VectorDomainBase<Field>::_F.areEqual (i->second, j->second))
				return false;

			i++; j++;
		}
	}

	return true;
}

template <class Field>
template <class Vector>
bool VectorDomain<Field>::isZeroSpecialized (const Vector &v, VectorCategories::DenseVectorTag) const
{
	typename Vector::const_iterator i;

	for (i = v.begin (); i != v.end (); i++)
		if (!VectorDomainBase<Field>::_F.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Vector>
bool VectorDomain<Field>::isZeroSpecialized (const Vector &v, VectorCategories::SparseVectorTag) const
{
	typename Vector::const_iterator i;

	for (i = v.begin (); i != v.end (); i++)
		if (!VectorDomainBase<Field>::_F.isZero (i->second))
			return false;

	return true;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
					       VectorCategories::SparseVectorTag,
					       VectorCategories::DenseVectorTag) const
{
	typename Vector2::const_iterator i;
	int idx;

	res.clear ();

	for (i = v.begin (), idx = 0; i != v.end (); i++, idx++)
		if (!VectorDomainBase<Field>::_F.isZero (*i))
			res.push_back (std::pair <size_t, typename Field::Element> (idx, *i));

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
					       VectorCategories::DenseVectorTag,
					       VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (v, res.size ()));

	typename Vector1::iterator i;
	typename Vector2::const_iterator j;
	size_t idx;

	for (i = res.begin (), j = v.begin (), idx = 0; j != v.end (); i++, j++, idx++) {
		while (idx < j->first) {
			VectorDomainBase<Field>::_F.init (*i, 0);
			i++; idx++;
		}

		*i = j->second;
	}

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
					       VectorCategories::DenseVectorTag) const
{
	if (i == 0 && len == 0) {
		copy (res, v);
	} else {
		Vector1 res_part;

		copy (res_part, v);

		std::copy (res_part.begin (), (len == 0) ? res_part.end () : res_part.begin () + len, res.begin () + i);
	}

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
					       VectorCategories::SparseVectorTag) const
{
	if (i == 0 && len == 0)
		return copy (res, v);

	typename Vector1::iterator r_begin, r_end, part_end, iter;
	Vector1 res_part;

	copy (res_part, v);

	if (len == 0)
		part_end = res_part.end ();
	else
		part_end = std::lower_bound (res_part.begin (), res_part.end (), len,
					     VectorWrapper::CompareSparseEntries ());

	for (iter = res_part.begin (); iter != part_end; iter++)
		iter->first += i;

	r_begin = std::lower_bound (res.begin (), res.end (), i, VectorWrapper::CompareSparseEntries ());
	r_end = (len == 0) ? res.end () : std::lower_bound (r_begin, res.end (), i + len,
							    VectorWrapper::CompareSparseEntries ());
	r_begin = res.erase (r_begin, r_end);
	res.insert (r_begin, res_part.begin (), part_end);

	return res;
}

template <class Field>
template <class Vector>
Vector &VectorDomain<Field>::copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
					      VectorCategories::DenseVectorTag) const
{
	std::copy (v.begin (), (len == 0) ? v.end () : v.begin () + len, res.begin () + i);
	return res;
}

template <class Field>
template <class Vector>
Vector &VectorDomain<Field>::copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
					      VectorCategories::SparseVectorTag) const
{
	typename Vector::const_iterator v_end;
	typename Vector::iterator r_begin, r_end, iter;
	typename Vector::difference_type offset;

	if (len == 0)
		v_end = v.end ();
	else
		v_end = std::lower_bound (v.begin (), v.end (), len,
					  VectorWrapper::CompareSparseEntries ());

	r_begin = std::lower_bound (res.begin (), res.end (), i, VectorWrapper::CompareSparseEntries ());
	r_end = (len == 0) ? res.end () : std::lower_bound (r_begin, res.end (), i + len,
							    VectorWrapper::CompareSparseEntries ());
	r_begin = res.erase (r_begin, r_end);
	offset = r_begin - res.begin ();
	res.insert (r_begin, v.begin (), v_end);
	r_begin = res.begin () + offset;
	r_end = r_begin + (v_end - v.begin ());

	for (iter = r_begin; iter != r_end; iter++)
		iter->first += i;

	return res;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					      VectorCategories::DenseVectorTag,
					      VectorCategories::DenseVectorTag,
					      VectorCategories::DenseVectorTag) const
{
	linbox_check (y.size () == x.size ());
	linbox_check (res.size () == x.size ());

	typename Vector2::const_iterator i;
	typename Vector3::const_iterator j;
	typename Vector1::iterator k;

	for (i = y.begin (), j = x.begin (), k = res.begin (); i != y.end (); i++, j++, k++)
		VectorDomainBase<Field>::_F.add (*k, *i, *j);

	return res;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					      VectorCategories::SparseVectorTag,
					      VectorCategories::SparseVectorTag,
					      VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator i;
	typename Vector3::const_iterator j;
	Element tmp;

	res.clear ();

	for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
		while (i != y.end () && i->first < j->first) {
			res.push_back (*i);
			i++;
		}

		if (i != y.end () && i->first == j->first) {
			VectorDomainBase<Field>::_F.add (tmp, i->second, j->second);
			if (!VectorDomainBase<Field>::_F.isZero (tmp))
				res.push_back (std::pair <size_t, Element> (j->first, tmp));
			i++;
		} else {
			res.push_back (*j);
		}
	}

	while (i != y.end ()) {
		res.push_back (*i);
		i++;
	}

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::addinSpecialized (Vector1 &y, const Vector2 &x,
						VectorCategories::DenseVectorTag,
						VectorCategories::DenseVectorTag) const
{
	linbox_check (y.size () == x.size ());

	typename Vector1::iterator i;
	typename Vector2::const_iterator j;

	for (i = y.begin (), j = x.begin (); i != y.end (); i++, j++)
		VectorDomainBase<Field>::_F.addin (*i, *j);

	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::addinSpecialized (Vector1 &y, const Vector2 &x,
						VectorCategories::SparseVectorTag,
						VectorCategories::SparseVectorTag) const
{
	Vector1 res;

	add (res, y, x);
	copy (y, res);
	return y;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					      VectorCategories::DenseVectorTag,
					      VectorCategories::DenseVectorTag,
					      VectorCategories::DenseVectorTag) const
{
	linbox_check (y.size () == x.size ());
	linbox_check (res.size () == x.size ());

	typename Vector2::const_iterator i;
	typename Vector3::const_iterator j;
	typename Vector1::iterator k;

	for (i = y.begin (), j = x.begin (), k = res.begin (); i != y.end (); i++, j++, k++)
		VectorDomainBase<Field>::_F.sub (*k, *i, *j);

	return res;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					      VectorCategories::SparseVectorTag,
					      VectorCategories::SparseVectorTag,
					      VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator i;
	typename Vector3::const_iterator j;
	Element tmp;

	res.clear ();

	for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
		while (i != y.end () && i->first < j->first) {
			res.push_back (*i);
			i++;
		}

		if (i != y.end () && i->first == j->first) {
			VectorDomainBase<Field>::_F.sub (tmp, i->second, j->second);
			if (!VectorDomainBase<Field>::_F.isZero (tmp))
				res.push_back (std::pair <size_t, Element> (j->first, tmp));
			i++;
		} else {
			res.push_back (std::pair <size_t, Element> (j->first, VectorDomainBase<Field>::_F.neg (tmp, j->second)));
		}
	}

	while (i != y.end ()) {
		res.push_back (*i);
		i++;
	}

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::subinSpecialized (Vector1 &y, const Vector2 &x,
						VectorCategories::DenseVectorTag,
						VectorCategories::DenseVectorTag) const
{
	linbox_check (y.size () == x.size ());

	typename Vector1::iterator i;
	typename Vector2::const_iterator j;

	for (i = y.begin (), j = x.begin (); i != y.end (); i++, j++)
		VectorDomainBase<Field>::_F.subin (*i, *j);

	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::subinSpecialized (Vector1 &y, const Vector2 &x,
						VectorCategories::SparseVectorTag,
						VectorCategories::SparseVectorTag) const
{
	Vector1 res;

	sub (res, y, x);
	copy (y, res);
	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::negSpecialized (Vector1 &res, const Vector2 &x,
					      VectorCategories::DenseVectorTag,
					      VectorCategories::DenseVectorTag) const
{
	linbox_check (res.size () == x.size ());

	typename Vector2::const_iterator j;
	typename Vector1::iterator k;

	for (j = x.begin (), k = res.begin (); j != x.end (); ++j, ++k)
		VectorDomainBase<Field>::_F.neg (*k, *j);

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::negSpecialized (Vector1 &res, const Vector2 &x,
					      VectorCategories::SparseVectorTag,
					      VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator j;
	Element tmp;

	res.clear ();

	for (j = x.begin (); j != x.end (); ++j)
		res.push_back (std::pair <size_t, Element> (j->first, VectorDomainBase<Field>::_F.neg (tmp, j->second)));

	return res;
}

template <class Field>
template <class Vector>
Vector &VectorDomain<Field>::neginSpecialized (Vector &y,
					       VectorCategories::DenseVectorTag) const
{
	typename Vector::iterator i;

	for (i = y.begin (); i != y.end (); ++i)
		VectorDomainBase<Field>::_F.negin (*i);

	return y;
}

template <class Field>
template <class Vector>
Vector &VectorDomain<Field>::neginSpecialized (Vector &y,
					       VectorCategories::SparseVectorTag) const
{
	typename Vector::iterator i;

	for (i = y.begin (); i != y.end (); ++i)
		VectorDomainBase<Field>::_F.negin (i->second);

	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::mulSpecialized
	(Vector1                       &res,
	 const Vector2                 &x,
	 const typename Field::Element &a,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (res.size () == x.size ());
	
	typename Vector2::const_iterator i;
	typename Vector1::iterator j;

	for (i = x.begin (), j = res.begin (); i != x.end (); ++i, ++j)
		VectorDomainBase<Field>::_F.mul (*j, *i, a);

	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::mulSpecialized
	(Vector1                       &res,
	 const Vector2                 &x,
	 const typename Field::Element &a,
	 VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator i;
	Element tmp;

	res.clear ();

	if (VectorDomainBase<Field>::_F.isZero (a))
		return res;

	for (i = x.begin (); i != x.end (); i++)
		res.push_back (std::pair <size_t, Element> (i->first, VectorDomainBase<Field>::_F.mul (tmp, i->second, a)));

	return res;
}

template <class Field>
template <class Vector>
Vector &VectorDomain<Field>::mulinSpecialized
	(Vector                        &x,
	 const typename Field::Element &a,
	 VectorCategories::DenseVectorTag) const
{
	typename Vector::iterator i;

	for (i = x.begin (); i != x.end (); i++)
		VectorDomainBase<Field>::_F.mulin (*i, a);

	return x;
}

template <class Field>
template <class Vector>
Vector &VectorDomain<Field>::mulinSpecialized
	(Vector                        &x,
	 const typename Field::Element &a,
	 VectorCategories::SparseVectorTag) const
{
	typename Vector::iterator i;

	if (VectorDomainBase<Field>::_F.isZero (a)) {
		x.clear ();
		return x;
	}

	for (i = x.begin (); i != x.end (); i++)
		VectorDomainBase<Field>::_F.mulin (i->second, a);

	return x;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::axpySpecialized
	(Vector1                       &res,
	 const Vector2                 &y,
	 const typename Field::Element &a,
	 const Vector3                 &x,
	 VectorCategories::DenseVectorTag,
	 VectorCategories::DenseVectorTag,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (y.size () == x.size ());
	linbox_check (res.size () == x.size ());

	typename Vector2::const_iterator i;
	typename Vector3::const_iterator j;
	typename Vector1::iterator k;

	for (i = y.begin (), j = x.begin (), k = res.begin (); i != y.end (); i++, j++, k++)
		VectorDomainBase<Field>::_F.axpy (*k, a, *j, *i);

	return res;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::axpySpecialized
	(Vector1                       &res,
	 const Vector2                 &y,
	 const typename Field::Element &a,
	 const Vector3                 &x,
	 VectorCategories::SparseVectorTag,
	 VectorCategories::SparseVectorTag,
	 VectorCategories::SparseVectorTag) const
{
	typename Vector2::const_iterator i;
	typename Vector3::const_iterator j;
	Element tmp;

	res.clear ();

	for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
		while (i != y.end () && i->first < j->first) {
			res.push_back (*i);
			i++;
		}

		if (i != y.end () && i->first == j->first) {
			VectorDomainBase<Field>::_F.axpy (tmp, a, j->second, i->second);
			i++;
		}
		else
			VectorDomainBase<Field>::_F.mul (tmp, a, j->second);

		if (!VectorDomainBase<Field>::_F.isZero (tmp))
			res.push_back (std::pair <size_t, Element> (j->first, tmp));
	}

	while (i != y.end ()) {
		res.push_back (*i);
		i++;
	}

	return res;
}

template <class Field>
template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<Field>::axpySpecialized
	(Vector1                       &res,
	 const Vector2                 &y,
	 const typename Field::Element &a,
	 const Vector3                 &x,
	 VectorCategories::DenseVectorTag,
	 VectorCategories::SparseVectorTag,
	 VectorCategories::DenseVectorTag) const
{
	copy (res, y);
	axpyin (res, a, x);
	return res;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::axpyinSpecialized
	(Vector1                       &y,
	 const typename Field::Element &a,
	 const Vector2                 &x,
	 VectorCategories::DenseVectorTag,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (y.size () == x.size ());

	typename Vector1::iterator i;
	typename Vector2::const_iterator j;

	for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
		VectorDomainBase<Field>::_F.axpyin (*i, a, *j);

	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::axpyinSpecialized
	(Vector1                       &y,
	 const typename Field::Element &a,
	 const Vector2                 &x,
	 VectorCategories::SparseVectorTag,
	 VectorCategories::DenseVectorTag) const
{
	Vector2 res (x.size ());

	axpy (res, a, x, y);
	copy (y, res);

	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::axpyinSpecialized
	(Vector1                       &y,
	 const typename Field::Element &a,
	 const Vector2                 &x,
	 VectorCategories::DenseVectorTag,
	 VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, y.size ()));

	typename Vector2::const_iterator j;

	for (j = x.begin (); j != x.end (); ++j)
		VectorDomainBase<Field>::_F.axpyin (y[j->first], a, j->second);

	return y;
}

template <class Field>
template <class Vector1, class Vector2>
Vector1 &VectorDomain<Field>::axpyinSpecialized
	(Vector1                       &y,
	 const typename Field::Element &a,
	 const Vector2                 &x,
	 VectorCategories::SparseVectorTag,
	 VectorCategories::SparseVectorTag) const
{
	Vector1 res;

	axpy (res, a, x, y);
	copy (y, res);
	return y;
}

template <class Field>
template <class Vector1, class Vector2>
inline typename Field::Element &DotProductDomain<Field>::dotSpecializedDD
	(Element       &res,
	 const Vector1 &v1,
	 const Vector2 &v2,
	 size_t         start_idx,
	 size_t         end_idx) const
{
	linbox_check (v1.size () == v2.size ());

	typename Vector1::const_iterator i, i_end = v1.begin () + std::min (v1.size (), end_idx);
	typename Vector2::const_iterator j;
	VectorDomainBase<Field>::accu.reset();

	for (i = v1.begin () + start_idx, j = v2.begin () + start_idx; i != i_end; ++i, ++j)
		VectorDomainBase<Field>::accu.mulacc (*i, *j);

	return VectorDomainBase<Field>::accu.get (res);
}

template <class Field>
template <class Vector1, class Vector2>
inline typename Field::Element &DotProductDomain<Field>::dotSpecializedDS
	(Element       &res,
	 const Vector1 &v1,
	 const Vector2 &v2,
	 size_t         start_idx,
	 size_t         end_idx) const
{
	linbox_check (VectorWrapper::hasDim<Field> (v1, v2.size ()));

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());
		
	VectorDomainBase<Field>::accu.reset();

	for (; i != i_end; ++i)
		VectorDomainBase<Field>::accu.mulacc (i->second, v2[i->first]);

	return VectorDomainBase<Field>::accu.get (res);
}

template <class Field>
template <class Vector1, class Vector2>
inline typename Field::Element &DotProductDomain<Field>::dotSpecializedSS
	(Element       &res,
	 const Vector1 &v1,
	 const Vector2 &v2,
	 size_t         start_idx,
	 size_t         end_idx) const
{
	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator j = (start_idx == 0) ? v2.begin () : std::lower_bound (v2.begin (), v2.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		v2.end () : std::lower_bound (v2.begin (), v2.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	VectorDomainBase<Field>::accu.reset();

	for (; i != i_end && j != j_end; ++i) {
		while (j != j_end && j->first < i->first) ++j;

		if (j != j_end && j->first == i->first)
			VectorDomainBase<Field>::accu.mulacc (i->second, j->second);
	}

	return VectorDomainBase<Field>::accu.get (res);
}

template <class Field>
template <class Vector>
inline void VectorDomain<Field>::swapSpecialized
	(Vector &v1, Vector &v2,
	 VectorCategories::DenseVectorTag) const 
{
	linbox_check (v1.size () == v2.size ());

	typename Vector::iterator j, k;

	for (j = v1.begin (), k = v2.begin (); j != v1.end (); ++j, ++k)
		std::swap (*j, *k);
}

template <class Field>
template <class Vector, class Iterator>
inline Vector &VectorDomain<Field>::permuteSpecialized
	(Vector &v, Iterator P_start, Iterator P_end,
	 VectorCategories::DenseVectorTag) const 
{
	Iterator i;

	for (i = P_start; i != P_end; ++i)
		std::swap (v[i->first], v[i->second]);

	return v;
}

template <class Field>
template <class Iterator>
typename Iterator::value_type::first_type VectorDomain<Field>::permutationImage (typename Iterator::value_type::first_type x, Iterator P_start, Iterator P_end) const
{
	for (Iterator i = P_start; i != P_end; ++i) {
		if (i->first == x)
			x = i->second;
		else if (i->second == x)
			x = i->first;
	}

	return x;
}

template <class Field>
template <class _Vector, class Iterator>
inline _Vector &VectorDomain<Field>::permuteSpecialized
	(_Vector &v, Iterator P_start, Iterator P_end,
	 VectorCategories::SparseVectorTag) const 
{
	typename _Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		j->first = permutationImage (j->first, P_start, P_end);

	std::sort (v.begin (), v.end (), VectorWrapper::CompareSparseEntries ());

	return v;
}

template <class Field>
template <class Vector>
inline int VectorDomain<Field>::firstNonzeroEntrySpecialized (typename Field::Element &a, const Vector &v,
							      VectorCategories::DenseVectorTag) const
{
	typename Vector::const_iterator i_v;

	for (i_v = v.begin (); i_v != v.end () && VectorDomainBase<Field>::_F.isZero (*i_v); ++i_v);

	if (i_v == v.end ())
		return -1;
	else {
		VectorDomainBase<Field>::_F.assign (a, *i_v);
		return i_v - v.begin ();
	}
}

template <class Field>
template <class Vector>
inline int VectorDomain<Field>::firstNonzeroEntrySpecialized (typename Field::Element &a, const Vector &v,
							      VectorCategories::SparseVectorTag) const
{
	if (v.empty ())
		return -1;
	else {
		VectorDomainBase<Field>::_F.assign (a, v.front ().second);
		return v.front ().first;
	}
}

} // namespace LinBox

#endif // __LINBOX_field_vector_domain_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
