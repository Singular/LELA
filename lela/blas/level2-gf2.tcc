/* lela/blas/level2-gf2.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for GF2
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_GF2_TCC
#define __BLAS_LEVEL2_GF2_TCC

#include <algorithm>
#include <iostream>

#include "lela/blas/level2-gf2.h"
#include "lela/blas/level1-ll.h"
#include "lela/blas/level2-ll.h"
#include "lela/blas/level3-ll.h"

namespace LELA
{

namespace BLAS2
{

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 MatrixIteratorTypes::Row,
							 VectorRepresentationTypes::Generic,
							 VectorRepresentationTypes::Dense01)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d = false;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		BLAS1::_dot<GF2, typename Modules::Tag>::op (F, M, d, *i_A, x);
		F.addin (y[idx], d);
	}
		
	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 MatrixIteratorTypes::Row,
							 VectorRepresentationTypes::Generic,
							 VectorRepresentationTypes::Sparse01)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	Vector2 t;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		BLAS1::_dot<GF2, typename Modules::Tag>::op (F, M, d, *i_A, x);

		if (d)
			t.push_back (idx);
	}

	BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, t, y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 MatrixIteratorTypes::Row,
							 VectorRepresentationTypes::Generic,
							 VectorRepresentationTypes::Hybrid01)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	Vector2 t;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		BLAS1::_dot<GF2, typename Modules::Tag>::op (F, M, d, *i_A, x);

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
							 MatrixIteratorTypes::Col,
							 VectorRepresentationTypes::Dense01,
							 VectorRepresentationTypes::Generic)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	typename Matrix::ConstColIterator i_A;
	typename Vector1::const_iterator i_x;

	for (i_x = x.begin (), i_A = A.colBegin (); i_A != A.colEnd (); ++i_x, ++i_A)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, *i_x, *i_A, y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 MatrixIteratorTypes::Col,
							 VectorRepresentationTypes::Sparse01,
							 VectorRepresentationTypes::Generic)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	typename Vector1::const_iterator i_x;

	for (i_x = x.begin (); i_x != x.end (); ++i_x)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, *(A.colBegin () + *i_x), y);

	return y;
}

template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<GF2, GenericModule<GF2>::Tag>::gemv_impl (const GF2 &F, Modules &M,
							 bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
							 MatrixIteratorTypes::Col,
							 VectorRepresentationTypes::Hybrid01,
							 VectorRepresentationTypes::Generic)
{
	// lela_check (VectorUtils::hasDim<GF2> (x, A.coldim ()));
	// lela_check (VectorUtils::hasDim<GF2> (y, A.rowdim ()));

	typename Matrix::ConstColIterator i_A;
	typename Vector1::const_iterator i_x;
	typename Vector1::word_type t;
	typename Vector1::index_type idx;

	if (!b)
		BLAS1::_scal<GF2, typename Modules::Tag>::op (F, M, false, y);

	if (!a)
		return y;

	for (i_x = x.begin (); i_x != x.end (); ++i_x) {
		t = Vector1::Endianness::e_0;
		idx = i_x->first << WordTraits<typename Vector1::word_type>::logof_size;

		i_A = A.colBegin () + idx;

		for (; t != 0 && idx < A.coldim (); t = Vector1::Endianness::shift_right (t, 1), ++i_A, ++idx)
			BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, ((i_x->second & t) != 0), *i_A, y);
	}

	return y;
}

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Dense01,
						      VectorRepresentationTypes::Generic,
						      MatrixIteratorTypes::Row)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.coldim ()));

	if (!a)
		return BLAS3::_scal<GF2, typename Modules::Tag>::op (F, M, false, A);

	typename Matrix::RowIterator i_A;
	typename Vector1::const_iterator i_x;

	for (i_A = A.rowBegin (), i_x = x.begin (); i_A != A.rowEnd (); ++i_A, ++i_x)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, *i_x, y, *i_A);

	return A;
}

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Sparse01,
						      VectorRepresentationTypes::Generic,
						      MatrixIteratorTypes::Row)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.coldim ()));

	if (!a)
		return BLAS3::_scal<GF2, typename Modules::Tag>::op (F, M, false, A);

	typename Vector1::const_iterator i_x;

	for (i_x = x.begin (); i_x != x.end (); ++i_x)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, y, *(A.rowBegin () + *i_x));

	return A;
}

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Hybrid01,
						      VectorRepresentationTypes::Generic,
						      MatrixIteratorTypes::Row)
{
	if (!a)
		return BLAS3::_scal<GF2, typename Modules::Tag>::op (F, M, false, A);

	typename Vector1::const_iterator i_x;
	typename Vector1::word_type t;
	size_t row;

	for (i_x = x.begin (); i_x != x.end (); ++i_x) {
		row = i_x->first << WordTraits<typename Vector1::word_type>::logof_size;

		for (t = Vector1::Endianness::e_0; t != 0 && row < A.rowdim (); t = Vector1::Endianness::shift_right (t, 1), ++row)
			BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, i_x->second & t, y, *(A.rowBegin () + row));
	}

	return A;
}

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Generic,
						      VectorRepresentationTypes::Dense01,
						      MatrixIteratorTypes::Col)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.coldim ()));

	if (!a)
		return BLAS3::_scal<GF2, typename Modules::Tag>::op (F, M, false, A);

	typename Matrix::ColIterator i_A;
	typename Vector2::const_iterator i_y;

	for (i_A = A.colBegin (), i_y = y.begin (); i_A != A.colEnd (); ++i_A, ++i_y)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, *i_y, x, *i_A);

	return A;
}

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Generic,
						      VectorRepresentationTypes::Sparse01,
						      MatrixIteratorTypes::Col)
{
	lela_check (VectorUtils::hasDim<GF2> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<GF2> (y, A.coldim ()));

	if (!a)
		return BLAS3::_scal<GF2, typename Modules::Tag>::op (F, M, false, A);

	typename Vector2::const_iterator i_y;

	for (i_y = y.begin (); i_y != y.end (); ++i_y)
		BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, true, y, *(A.colBegin () + *i_y));

	return A;
}

template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<GF2, GenericModule<GF2>::Tag>::ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
						      VectorRepresentationTypes::Generic,
						      VectorRepresentationTypes::Hybrid01,
						      MatrixIteratorTypes::Col)
{
	if (!a)
		return BLAS3::_scal<GF2, typename Modules::Tag>::op (F, M, false, A);

	typename Vector2::const_iterator i_y;
	typename Vector2::word_type t;
	size_t col;

	for (i_y = y.begin (); i_y != y.end (); ++i_y) {
		col = i_y->first << WordTraits<typename Vector1::word_type>::logof_size;

		for (t = Vector1::Endianness::e_0; t != 0 && col < A.coldim (); t = Vector1::Endianness::shift_right (t, 1), ++col)
			BLAS1::_axpy<GF2, typename Modules::Tag>::op (F, M, i_y->second & t, y, *(A.colBegin () + col));
	}

	return A;
}

} // namespace BLAS2

} // namespace LELA

#include "lela/blas/level2-gf2.tcc"

#endif // __BLAS_LEVEL2_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
