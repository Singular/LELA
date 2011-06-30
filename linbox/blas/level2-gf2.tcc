/* linbox/blas/level2-gf2.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for GF2
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_GF2_TCC
#define __BLAS_LEVEL2_GF2_TCC

#include <algorithm>
#include <iostream>

#include "linbox/blas/level2-gf2.h"
#include "linbox/blas/level1-ll.h"
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

namespace BLAS2
{

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 size_t start_idx, size_t end_idx,
							 MatrixIteratorTypes::Row,
							 VectorRepresentationTypes::Generic,
							 VectorRepresentationTypes::Dense01)
{
	linbox_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d = false;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		BLAS1::_dot<GF2, typename Modules::Tag>::op (F, M, d, *i_A, x, start_idx, end_idx);
		F.addin (y[idx], d);
	}
		
	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 size_t start_idx, size_t end_idx,
							 MatrixIteratorTypes::Row,
							 VectorRepresentationTypes::Generic,
							 VectorRepresentationTypes::Sparse01)
{
	linbox_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	Vector2 t;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		BLAS1::_dot<GF2, typename Modules::Tag>::op (F, M, d, *i_A, x, start_idx, end_idx);

		if (d)
			t.push_back (idx);
	}

	BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, t, y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 size_t start_idx, size_t end_idx,
							 MatrixIteratorTypes::Row,
							 VectorRepresentationTypes::Generic,
							 VectorRepresentationTypes::Hybrid01)
{
	linbox_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	Vector2 t;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		BLAS1::_dot<GF2, typename Modules::Tag>::op (F, M, d, *i_A, x, start_idx, end_idx);

		if (d) {
			if (t.first.empty () || t.first.back () != idx & ~WordTraits<typename Vector2::word_type>::pos_mask)
				t.push_back (typename Vector2::value_type (idx >> WordTraits<typename Vector2::word_type>::logof_size,
									   1ULL << (idx & WordTraits<typename Vector2::word_type>::pos_mask)));
			else
				t.back ().second |= 1ULL << (idx & WordTraits<typename Vector2::word_type>::pos_mask);
		}
	}

	BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, t, y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 size_t start_idx, size_t end_idx,
							 MatrixIteratorTypes::Col,
							 VectorRepresentationTypes::Dense01,
							 VectorRepresentationTypes::Generic)
{
	linbox_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	linbox_check (end_idx <= A.coldim ());

	typename Matrix::ConstColIterator i_A, i_A_end = A.colBegin () + end_idx;
	typename Vector1::const_iterator i_x;

	for (i_x = x.begin () + start_idx, i_A = A.colBegin () + start_idx; i_A != i_A_end; ++i_x, ++i_A)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, *i_x, *i_A, y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 size_t start_idx, size_t end_idx,
							 MatrixIteratorTypes::Col,
							 VectorRepresentationTypes::Sparse01,
							 VectorRepresentationTypes::Generic)
{
	linbox_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Vector1::const_iterator i_x = std::lower_bound (x.begin (), x.end (), start_idx);

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	linbox_check (end_idx <= A.coldim ());

	for (; i_x != x.end () && *i_x < end_idx; ++i_x)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, *(A.colBegin () + *i_x), y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 size_t start_idx, size_t end_idx,
							 MatrixIteratorTypes::Col,
							 VectorRepresentationTypes::Hybrid01,
							 VectorRepresentationTypes::Generic)
{
	// linbox_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	// linbox_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstColIterator i_A;
	typename Vector1::const_iterator i_x;
	typename Vector1::word_type t;
	typename Vector1::index_type idx;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	linbox_check (end_idx <= A.coldim ());

	i_x = std::lower_bound (x.begin (), x.end (), start_idx >> WordTraits<typename Vector1::word_type>::logof_size, VectorUtils::CompareSparseEntries ());

	for (; i_x != x.end () && (i_x->first << WordTraits<typename Vector1::word_type>::logof_size) < (long) end_idx; ++i_x) {
		if (start_idx >> WordTraits<typename Vector1::word_type>::logof_size == i_x->first) {
			t = Vector1::Endianness::e_j (start_idx & WordTraits<typename Vector1::word_type>::pos_mask);
			idx = start_idx;
		} else {
			t = Vector1::Endianness::e_0;
			idx = i_x->first << WordTraits<typename Vector1::word_type>::logof_size;
		}

		i_A = A.colBegin () + idx;

		for (; t != 0 && idx < end_idx; t = Vector1::Endianness::shift_right (t, 1), ++i_A, ++idx)
			BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, ((i_x->second & t) != 0), *i_A, y);
	}

	return y;
}

#if 0 // Not implemented

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Dense01,
						      VectorRepresentationTypes::Dense01,
						      MatrixIteratorTypes::Row);

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Sparse01,
						      VectorRepresentationTypes::Sparse01,
						      MatrixIteratorTypes::Row);

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Hybrid01,
						      VectorRepresentationTypes::Hybrid01,
						      MatrixIteratorTypes::Row);

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Dense01,
						      VectorRepresentationTypes::Dense01,
						      MatrixIteratorTypes::Col);

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Sparse01,
						      VectorRepresentationTypes::Sparse01,
						      MatrixIteratorTypes::Col);

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Hybrid01,
						      VectorRepresentationTypes::Hybrid01,
						      MatrixIteratorTypes::Col);

#endif // Not implemented

} // namespace BLAS2

} // namespace LinBox

#include "linbox/blas/level2-gf2.tcc"

#endif // __BLAS_LEVEL2_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
