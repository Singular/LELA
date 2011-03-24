/* linbox/field/modular.inl
 * Copyright (C) 2002 Bradford Hovinen
 * Copyright (C) 2002 Ahmet Duran
 * Copyright (C) 2002 B. David Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Ahmet Duran <duran@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_modular_INL
#define __LINBOX_field_modular_INL

//Dan Roche 7-2-04
#ifndef __LINBOX_MIN
#define __LINBOX_MIN(a,b) ( (a) < (b) ? (a) : (b) )
#endif

#include <iostream>

namespace LinBox
{

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDD
	(uint8 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();

	typename Vector1::const_iterator iterend = v1.begin () + v1.size() % _F._k;

	uint64 y = 0;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != v1.end (); j += _F._k) {
		typename Vector1::const_iterator iter_i = iterend;
		typename Vector2::const_iterator iter_j;

		iterend += _F._k;

		for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
			y += (uint64) *iter_i * (uint64) *j;

		y %= (uint64) _F._modulus;
	}

	return res = y;
}

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDSP
	(uint8 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
	typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();

	uint64 y = 0;

	if (v1.first.size () < _F._k) {
		for (; i_idx != v1.first.end (); ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		return res = y % (uint64) _F._modulus;
	} else {
		typename Vector1::first_type::const_iterator iterend = v1.first.begin () + v1.first.size() % _F._k;

		for (; i_idx != iterend; ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		y %= (uint64) _F._modulus;

		while (iterend != v1.first.end ()) {
			typename Vector1::first_type::const_iterator iter_i_idx = iterend;
			typename Vector1::second_type::const_iterator iter_i_elt = i_elt;

			iterend += _F._k;
			i_elt += _F._k;

			for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
				y += (uint64) *iter_i_elt * (uint64) v2[*iter_i_idx];

			y %= (uint64) _F._modulus;
		}

		return res = y;
	}
}

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDS
	(uint8 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();

	uint64 y = 0;

	if (v1.size () < _F._k) {
		for (; i != v1.end (); ++i)
			y += (uint64) i->second * (uint64) v2[i->first];

		return res = y % (uint64) _F._modulus;
	} else {
		typename Vector1::const_iterator iterend = v1.begin () + v1.size() % _F._k;

		for (; i != iterend; ++i)
			y += (uint64) i->second * (uint64) v2[i->first];

		y %= (uint64) _F._modulus;

		while (iterend != v1.end ()) {
			typename Vector1::const_iterator iter_i = iterend;

			iterend += _F._k;
			i += _F._k;

			for (; iter_i != iterend; ++iter_i)
				y += (uint64) iter_i->second * (uint64) v2[iter_i->first];

			y %= (uint64) _F._modulus;
		}

		return res = y;
	}
}

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedSS
	(uint8 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();

	uint64 y = 0, count;

	while (i != v1.end () && j != v2.end ()) {
		for (count = 0; count < _F._k && i != v1.end () && j != v2.end (); ++i) {
			while (j != v2.end () && j->first < i->first) ++j;

			if (j != v2.end () && i->first == j->first) {
				y += (uint64) i->second * (uint64) j->second;
				++count;
			}
		}

		y %= (uint64) _F._modulus;
	}

	return res = y;
}

template <class Vector1, class Vector2>
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedDD
	(uint16 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();

	typename Vector1::const_iterator iterend = v1.begin () + v1.size() % _F._k;

	uint64 y = 0;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != v1.end (); j += _F._k) {
		typename Vector1::const_iterator iter_i = iterend;
		typename Vector2::const_iterator iter_j;

		iterend += _F._k;

		for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
			y += (uint64) *iter_i * (uint64) *j;

		y %= (uint64) _F._modulus;
	}

	return res = y;
}

template <class Vector1, class Vector2>
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedDSP
	(uint16 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
	typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();

	uint64 y = 0;

	if (v1.first.size () < _F._k) {
		for (; i_idx != v1.first.end (); ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		return res = y % (uint64) _F._modulus;
	} else {
		typename Vector1::first_type::const_iterator iterend = v1.first.begin () + v1.first.size() % _F._k;

		for (; i_idx != iterend; ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		y %= (uint64) _F._modulus;

		while (iterend != v1.first.end ()) {
			typename Vector1::first_type::const_iterator iter_i_idx = iterend;
			typename Vector1::second_type::const_iterator iter_i_elt = i_elt;

			iterend += _F._k;
			i_elt += _F._k;

			for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
				y += (uint64) *iter_i_elt * (uint64) v2[*iter_i_idx];

			y %= (uint64) _F._modulus;
		}

		return res = y;
	}
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDD
	(uint32 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;
  
	uint64 y = 0;
	uint64 t;

	for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
		t = (uint64) *i * (uint64) *j;
		y += t;

		if (y < t)
			y += _F._two_64;
	}
  
	y %= (uint64) _F._modulus;

	return res = y;
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDSP
	(uint32 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::first_type::const_iterator i_idx;
	typename Vector1::second_type::const_iterator i_elt;
  
	uint64 y = 0;
	uint64 t;

	for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
		t = (uint64) *i_elt * (uint64) v2[*i_idx];
		y += t;

		if (y < t)
			y += _F._two_64;
	}
  
	y %= (uint64) _F._modulus;

	return res = y;
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDS
	(uint32 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();

	uint64 y = 0, t;

	for (i; i != v1.end (); ++i) {
		t = (uint64) i->second * (uint64) v2[i->first];
		y += t;

		if (y < t)
			y += _F._two_64;
	}
  
	y %= (uint64) _F._modulus;

	return res = y;
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedSS
	(uint32 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();

	uint64 y = 0, t;

	for (; i != v1.end () && j != v2.end (); ++i) {
		while (j != v2.end () && j->first < i->first) ++j;

		if (j != v2.end () && i->first == j->first) {
			t = (uint64) i->second * (uint64) j->second;
			y += t;

			if (y < t)
				y += _F._two_64;
		}
	}

	return res = y % (uint64) _F._modulus;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint8> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint8> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint8> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint8> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::first_type::const_iterator k_idx;
	typename Matrix::Column::second_type::const_iterator k_elt;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (uint64 (A.coldim ()), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
				_tmp[*k_idx] += *k_elt * *j;

		j_end += __LINBOX_MIN (uint64 (A.coldim () - (j_end - x.begin ())), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint16> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = x.begin (), j_end;
	typename Matrix::Column::const_iterator k;
	// Dan Roche, 7-1-04
	// std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint16> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint16> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint16> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::first_type::const_iterator k_idx;
	typename Matrix::Column::second_type::const_iterator k_elt;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin ();
		//Dan Roche, 7-2-04
		//j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);
		j_end = j + __LINBOX_MIN (A.coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
				_tmp[*k_idx] += *k_elt * *j;

		//j_end += __LINBOX_MIN (A->coldim () - (j_end - x.begin ()), VD.field ()._k);
		j_end += __LINBOX_MIN (A.coldim () - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint32> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	for (j = x.begin (); j != x.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) *k) * ((uint64) *j);

			*l += t;

			if (*l < t)
				*l += VD.field ()._two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint32> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	for (j = x.begin (); j != x.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) k->second) * ((uint64) *j);

			_tmp[k->first] += t;

			if (_tmp[k->first] < t)
				_tmp[k->first] += VD.field ()._two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint32> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	for (j = x.begin (); j != x.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) k->second) * ((uint64) *j);

			_tmp[k->first] += t;

			if (_tmp[k->first] < t)
				_tmp[k->first] += VD.field ()._two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint32> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.coldim () == x.size ());
	linbox_check (A.rowdim () == y.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::first_type::const_iterator k_idx;
	typename Matrix::Column::second_type::const_iterator k_elt;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	for (j = x.begin (); j != x.end (); ++j, ++i) {
		for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
		     k_idx != i->first.end ();
		     ++k_idx, ++k_elt, ++l)
		{
			t = ((uint64) *k_elt) * ((uint64) *j);

			_tmp[*k_idx] += t;

			if (_tmp[*k_idx] < t)
				_tmp[*k_idx] += VD.field ()._two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

} // namespace LinBox

#endif // __LINBOX_field_modular_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
