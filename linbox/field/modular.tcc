/* linbox/field/modular.tcc
 * Copyright 2002 Bradford Hovinen
 * Copyright 2002 Ahmet Duran
 * Copyright 2002 B. David Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Ahmet Duran <duran@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_modular_TCC
#define __LINBOX_field_modular_TCC

//Dan Roche 7-2-04
#ifndef __LINBOX_MIN
#define __LINBOX_MIN(a,b) ( (a) < (b) ? (a) : (b) )
#endif

#include <iostream>

namespace LinBox
{

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDD
	(uint8 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (v1.size () == v2.size ());
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = v1.begin () + start_idx, i_end = v1.begin () + std::min (v1.size (), end_idx);
	typename Vector2::const_iterator j = v2.begin () + start_idx;

	typename Vector1::const_iterator iterend = v1.begin () + (start_idx + std::min (v1.size () - start_idx, end_idx - start_idx) % _F._k);

	uint64 y = 0;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != i_end; j += _F._k) {
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
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDS
	(uint8 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (v1, v2.size ()));
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	uint64 y = 0;

	if (i_end - i < (long) _F._k) {
		for (; i != i_end; ++i)
			y += (uint64) i->second * (uint64) v2[i->first];

		return res = y % (uint64) _F._modulus;
	} else {
		// i still points to the beginning
		typename Vector1::const_iterator iterend = i + (i_end - i) % _F._k;

		for (; i != iterend; ++i)
			y += (uint64) i->second * (uint64) v2[i->first];

		y %= (uint64) _F._modulus;

		while (iterend != i_end) {
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
	(uint8 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator j = (start_idx == 0) ? v2.begin () : std::lower_bound (v2.begin (), v2.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		v2.end () : std::lower_bound (v2.begin (), v2.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	uint64 y = 0, count;

	while (i != i_end && j != j_end) {
		for (count = 0; count < _F._k && i != i_end && j != j_end; ++i) {
			while (j != j_end && j->first < i->first) ++j;

			if (j != j_end && i->first == j->first) {
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
	(uint16 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (v1.size () == v2.size ());
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = v1.begin () + start_idx, i_end = v1.begin () + std::min (v1.size (), end_idx);
	typename Vector2::const_iterator j = v2.begin () + start_idx;

	typename Vector1::const_iterator iterend = v1.begin () + (start_idx + std::min (v1.size () - start_idx, end_idx - start_idx) % _F._k);

	uint64 y = 0;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != i_end; j += _F._k) {
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
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedDS
	(uint16 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (v1, v2.size ()));
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	uint64 y = 0;

	if (i_end - i < (long) _F._k) {
		for (; i != i_end; ++i)
			y += (uint64) i->second * (uint64) v2[i->first];

		return res = y % (uint64) _F._modulus;
	} else {
		// i still points to the beginning
		typename Vector1::const_iterator iterend = i + (i_end - i) % _F._k;

		for (; i != iterend; ++i)
			y += (uint64) i->second * (uint64) v2[i->first];

		y %= (uint64) _F._modulus;

		while (iterend != i_end) {
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
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedSS
	(uint16 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator j = (start_idx == 0) ? v2.begin () : std::lower_bound (v2.begin (), v2.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		v2.end () : std::lower_bound (v2.begin (), v2.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	uint64 y = 0, count;

	while (i != i_end && j != j_end) {
		for (count = 0; count < _F._k && i != i_end && j != j_end; ++i) {
			while (j != j_end && j->first < i->first) ++j;

			if (j != j_end && i->first == j->first) {
				y += (uint64) i->second * (uint64) j->second;
				++count;
			}
		}

		y %= (uint64) _F._modulus;
	}

	return res = y;
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDD
	(uint32 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (v1.size () == v2.size ());
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = v1.begin () + start_idx, i_end = v1.begin () + std::min (v1.size (), end_idx);
	typename Vector2::const_iterator j = v2.begin () + start_idx;
  
	uint64 y = 0;
	uint64 t;

	for (; i != i_end; ++i, ++j) {
		t = (uint64) *i * (uint64) *j;
		y += t;

		if (y < t)
			y += _F._two_64;
	}
  
	y %= (uint64) _F._modulus;

	return res = y;
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDS
	(uint32 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (v1, v2.size ()));
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	uint64 y = 0, t;

	for (; i != i_end; ++i) {
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
	(uint32 &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const
{
	linbox_check (start_idx <= end_idx);

	typename Vector1::const_iterator i = (start_idx == 0) ? v1.begin () : std::lower_bound (v1.begin (), v1.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector2::const_iterator j = (start_idx == 0) ? v2.begin () : std::lower_bound (v2.begin (), v2.end (), start_idx, VectorWrapper::CompareSparseEntries ());

	typename Vector1::const_iterator i_end = (end_idx == static_cast<size_t> (-1)) ?
		v1.end () : std::lower_bound (v1.begin (), v1.end (), end_idx, VectorWrapper::CompareSparseEntries ());
	typename Vector1::const_iterator j_end = (end_idx == static_cast<size_t> (-1)) ?
		v2.end () : std::lower_bound (v2.begin (), v2.end (), end_idx, VectorWrapper::CompareSparseEntries ());

	uint64 y = 0, t;

	for (; i != i_end && j != j_end; ++i) {
		while (j != j_end && j->first < i->first) ++j;

		if (j != j_end && i->first == j->first) {
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
	 size_t start_idx, size_t end_idx,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + __LINBOX_MIN (end_idx - start_idx, VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += __LINBOX_MIN (end_idx - start_idx - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint8> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 size_t start_idx, size_t end_idx,
	 VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + __LINBOX_MIN (end_idx - start_idx, VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (end_idx - start_idx - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint16> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 size_t start_idx, size_t end_idx,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j = x.begin () + start_idx, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	// Dan Roche, 7-1-04
	// std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + __LINBOX_MIN (end_idx - start_idx, VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += __LINBOX_MIN (end_idx - start_idx - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint16> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 size_t start_idx, size_t end_idx,
	 VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	l_end = _tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + __LINBOX_MIN (end_idx - start_idx, VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (end_idx - start_idx - (j_end - x.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = _tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (alpha * *l + beta * *y_j) % VD.field ()._modulus;

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Modular<uint32> >::gemvColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, const Element &alpha, const Matrix &A, const Vector1 &x, const Element &beta, Vector2 &y,
	 size_t start_idx, size_t end_idx,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	for (j = x.begin () + start_idx; j != j_stop; ++j, ++i) {
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
	 size_t start_idx, size_t end_idx,
	 VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < y.size ())
		_tmp.resize (y.size ());

	std::fill (_tmp.begin (), _tmp.begin () + y.size (), 0);

	for (j = x.begin () + start_idx; j != j_stop; ++j, ++i) {
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

} // namespace LinBox

#endif // __LINBOX_field_modular_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
