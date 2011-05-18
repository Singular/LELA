/* linbox/blas/level2-modular.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for Z/p
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_MODULAR_TCC
#define __BLAS_LEVEL2_MODULAR_TCC

#include <algorithm>

#include "linbox/blas/level2-modular.h"

namespace LinBox
{

namespace BLAS2
{

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
			 uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l, l_end;

	M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + std::min (end_idx - start_idx, F._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += std::min (end_idx - start_idx - (j_end - x.begin ()), F._k);

		for (l = M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
			 uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::SparseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint8> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l, l_end;

	M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + std::min (end_idx - start_idx, F._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				M._tmp[k->first] += k->second * *j;

		j_end += std::min (end_idx - start_idx - (j_end - x.begin ()), F._k);

		for (l =M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
			 uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint16> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j = x.begin () + start_idx, j_end, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l, l_end;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + std::min (end_idx - start_idx, F._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += std::min (end_idx - start_idx - (j_end - x.begin ()), F._k);

		for (l = M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
			 uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::SparseVectorTag)
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

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin () + start_idx;
		j_end = j + std::min (end_idx - start_idx, F._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				M._tmp[k->first] += k->second * *j;

		j_end += std::min (end_idx - start_idx - (j_end - x.begin ()), F._k);

		for (l =M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != j_stop);

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
			 uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	for (j = x.begin () + start_idx; j != j_stop; ++j, ++i) {
		for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) *k) * ((uint64) *j);

			*l += t;

			if (*l < t)
				*l += F._two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
			 uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
			 size_t start_idx, size_t end_idx,
			 VectorCategories::SparseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Modular<uint32> > (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector2::const_iterator j, j_stop = x.begin () + end_idx;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	for (j = x.begin () + start_idx; j != j_stop; ++j, ++i) {
		for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) k->second) * ((uint64) *j);

			M._tmp[k->first] += t;

			if (M._tmp[k->first] < t)
				M._tmp[k->first] += F._two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_MODULAR_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
