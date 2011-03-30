/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/matrix-domain-gf2.inl
 * Copyright 2010 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __MATRIX_MATRIX_DOMAIN_GF2_INL
#define __MATRIX_MATRIX_DOMAIN_GF2_INL

#include "linbox/util/error.h"

namespace LinBox
{

template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomainSupportGF2::copyRow (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.copy (*i, *j);

	return A;
}

template <class Matrix1, class Matrix2>
inline bool MatrixDomainSupportGF2::areEqualRow (const Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		if (!_VD.areEqual (*i, *j))
			return false;

	return true;
}

template <class Matrix>
bool MatrixDomainSupportGF2::isZeroRow (const Matrix &A) const
{
	typename Matrix::ConstRowIterator i;

	i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Matrix>
Matrix &MatrixDomainSupportGF2::scalRow (Matrix &A, const bool &a) const
{
	if (!a) {
		typename Matrix::RowIterator i;

		for (i = A.rowBegin (); i != A.rowEnd (); ++i)
			_VD.mulin (*i, a);
	}

	return A;
}

template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGF2::axpyRow (const bool &a, const Matrix1 &A, Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::RowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		_VD.axpyin (*j, a, *i);

	return B;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGF2::gemvColDense (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y) const
{
	typename Matrix::ConstColIterator i_A;
	typename Vector1::const_iterator i_x;

	if (!b)
		_VD.subin (y, y);

	if (!a)
		return y;

	for (i_x = x.begin (), i_A = A.colBegin (); i_x != x.end (); ++i_x, ++i_A)
		if (*i_x)
			_VD.addin (y, *i_A);

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGF2::gemvColSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
						     VectorCategories::SparseZeroOneVectorTag,
						     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector1::const_iterator i_x;

	if (!b)
		_VD.subin (y, y);

	if (!a)
		return y;

	for (i_x = x.begin (); i_x != x.end (); ++i_x)
		_VD.addin (y, *(A.colBegin () + *i_x));

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGF2::gemvColSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
						     VectorCategories::HybridZeroOneVectorTag,
						     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Matrix::ConstColIterator i_A;
	typename Vector1::const_iterator i_x;
	typename Vector1::word_type t;

	if (!b)
		_VD.subin (y, y);

	if (!a)
		return y;

	for (i_x = x.begin (), i_A = A.colBegin (); i_x != x.end (); ++i_x)
		for (t = Vector1::Endianness::e_0, i_A = A.colBegin () + i_x->first; t != 0; t = Vector1::Endianness::shift_right (t, 1), ++i_A)
			if (i_x->second & t)
				_VD.addin (y, *i_A);

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGF2::gemvRowSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
						     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	if (!b)
		_VD.subin (y, y);

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		_VD.dot (d, *i_A, x);
		_F.addin (y[idx], d);
	}
		
	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGF2::gemvRowSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
						     VectorCategories::SparseZeroOneVectorTag) const
{
	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	Vector2 t;

	if (!b)
		y.clear ();

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		_VD.dot (d, *i_A, x);

		if (d)
			t.push_back (idx);
	}

	_VD.addin (y, t);

	return y;
}

template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGF2::gemvRowSpecialized (const bool &a, const Matrix &A, const Vector1 &x, const bool &b, Vector2 &y,
						     VectorCategories::HybridZeroOneVectorTag) const
{
	typename Matrix::ConstRowIterator i_A;
	size_t idx;

	bool d;

	Vector2 t;

	if (!b) {
		y.clear ();
	}

	if (!a)
		return y;

	for (i_A = A.rowBegin (), idx = 0; i_A != A.rowEnd (); ++i_A, ++idx) {
		_VD.dot (d, *i_A, x);

		if (d) {
			if (t.first.empty () || t.first.back () != idx & ~WordTraits<typename Vector2::word_type>::pos_mask)
				t.push_back (typename Vector2::value_type (idx >> WordTraits<typename Vector2::word_type>::logof_size,
									   1ULL << (idx & WordTraits<typename Vector2::word_type>::pos_mask)));
			else
				t.back ().second |= 1ULL << (idx & WordTraits<typename Vector2::word_type>::pos_mask);
		}
	}

	_VD.addin (y, t);

	return y;
}

template <class Matrix, class Vector>
Vector &MatrixDomainSupportGF2::trsvSpecialized (const Matrix &A, Vector &x,
						 MatrixCategories::RowMatrixTag,
						 VectorCategories::DenseZeroOneVectorTag) const
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (A.rowdim () == x.size ());

	bool ai;
	int i = A.rowdim () - 1;

	while (--i >= 0) {
		if (_VD.firstNonzeroEntry (ai, *(A.rowBegin () + i)) == -1)
			continue;

		_VD.dot (x[i], *(A.rowBegin () + i), x);
	}

	return x;
}

template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGF2::gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
							   VectorCategories::DenseZeroOneVectorTag) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;
	typename Matrix3::RowIterator k;
	typename Matrix1::ConstRow::const_iterator ip;

	for (i = A.rowBegin (), k = C.rowBegin (); i != A.rowEnd (); ++i, ++k) {
		if (!b)
			_VD.mulin (*k, false);

		if (a)
			for (j = B.rowBegin (), ip = i->begin (); j != B.rowEnd (); ++j, ++ip)
				if (*ip)
					_VD.addin (*k, *j);
	}

	return C;
}

template <class Matrix1, class Matrix3>
Matrix3 &MatrixDomainSupportGF2::gemmRowRowRow (const bool &a, const Matrix1 &A, const Submatrix<DenseZeroOneMatrix<> > &B, const bool &b, Matrix3 &C) const
{
	typedef typename DenseZeroOneMatrix<>::word_type word_type;

	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	if (!a) {
		if (!b)
			scal (C, false);

		return C;
	}

	size_t raw_coldim = (B.startCol () & WordTraits<word_type>::pos_mask) + B.coldim ();
	size_t coldim = (!(raw_coldim & WordTraits<word_type>::pos_mask)) ? raw_coldim : ((raw_coldim & ~WordTraits<word_type>::pos_mask) + WordTraits<word_type>::bits);

	DenseZeroOneMatrix<> &Cp = *const_cast<DenseZeroOneMatrix<> *> (&_Cp);
	Cp.resize (C.rowdim (), coldim);

	DenseZeroOneMatrix<> Bp (B.parent (), B.startRow (), B.startCol () & ~WordTraits<word_type>::pos_mask, B.rowdim (), coldim);

	gemm (a, A, Bp, b, Cp);

	Submatrix<DenseZeroOneMatrix<> > Cpp (Cp, 0, B.startCol () & WordTraits<word_type>::pos_mask, C.rowdim (), B.coldim ());
	copy (C, Cpp);

	return C;
}

template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGF2::gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
							   VectorCategories::SparseZeroOneVectorTag) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;
	typename Matrix3::RowIterator k;
	typename Matrix1::ConstRow::const_iterator ip;

	for (i = A.rowBegin (), k = C.rowBegin (); i != A.rowEnd (); ++i, ++k) {
		if (_F.isZero (b))
			_VD.mulin (*k, false);

		if (_F.isOne (a))
			for (ip = i->begin (); ip != i->end (); ++ip)
				_VD.addin (*k, *(B.rowBegin () + *ip));
	}

	return C;
}

template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGF2::gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C,
							   VectorCategories::HybridZeroOneVectorTag) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;
	typename Matrix3::RowIterator k;
	typename Matrix1::ConstRow::const_iterator ip;

	typename Matrix1::ConstRow::word_type mask;

	for (i = A.rowBegin (), k = C.rowBegin (); i != A.rowEnd (); ++i, ++k) {
		if (_F.isZero (b))
			_VD.mulin (*k, false);

		if (_F.isOne (a)) {
			for (ip = i->begin (); ip != i->end (); ++ip) {
				j = B.rowBegin () + (static_cast<size_t> (ip->first) << WordTraits<typename Matrix1::ConstRow::word_type>::logof_size);

				for (mask = Matrix1::ConstRow::Endianness::e_0; mask != 0 && j != B.rowEnd (); mask = Matrix1::ConstRow::Endianness::shift_right (mask, 1), ++j)
					if (ip->second & mask)
						_VD.addin (*k, *j);
			}
		}
	}

	return C;
}

template <class Matrix1, class Matrix3>
Matrix3 &MatrixDomainSupportGF2::gemmRowRowRowSpecialised (const bool &a, const Matrix1 &A, const Submatrix<DenseZeroOneMatrix<> > &B, const bool &b, Matrix3 &C,
							   VectorCategories::HybridZeroOneVectorTag) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Submatrix<DenseZeroOneMatrix<> >::ConstRowIterator j;
	typename Matrix3::RowIterator k;
	typename Matrix1::ConstRow::const_iterator ip;
	size_t t;

	typename Matrix1::ConstRow::word_type mask;

	for (i = A.rowBegin (), k = C.rowBegin (); i != A.rowEnd (); ++i, ++k) {
		if (_F.isZero (b))
			_VD.mulin (*k, false);

		if (_F.isOne (a)) {
			for (ip = i->begin (); ip != i->end (); ++ip) {
				t = static_cast<size_t> (ip->first) << WordTraits<typename Matrix1::ConstRow::word_type>::logof_size;
				j = B.rowBegin () + t;

				for (mask = Matrix1::ConstRow::Endianness::e_0; mask != 0 && t < B.rowdim (); mask = Matrix1::ConstRow::Endianness::shift_right (mask, 1), ++j, ++t)
					if (ip->second & mask)
						_VD.addin (*k, *j);
			}
		}
	}

	return C;
}

template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGF2::trsmSpecialized (const bool &a, const Matrix1 &A, Matrix2 &B,
						  MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag) const
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	bool ai;
	int i = A.rowdim () - 1;

	TransposeMatrix<const Matrix1> AT (A);

	while (--i >= 0) {
		if (_VD.firstNonzeroEntry (ai, *(A.rowBegin () + i)) == -1)
			continue;

		gemv (true, AT, *(A.rowBegin () + i), false, *(B.rowBegin () + i));
	}

	return B;
}

template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGF2::trsmSpecialized (const bool &a, const Matrix1 &A, Matrix2 &B,
						  MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag) const
{
	typename Matrix2::ColIterator i_B;

	if (!a) {
		subin (B, B);
		return B;
	}

	for (i_B = B.colBegin (); i_B != B.colEnd (); ++i_B)
		trsv (A, *i_B);

	return B;
}

template <class Matrix, class Iterator>
inline Matrix &MatrixDomainSupportGF2::permuteRowsByRow (Matrix   &A,
							 Iterator  P_start,
							 Iterator  P_end) const
{
	Iterator i;
	typename Matrix::RowIterator j, k;

	for (i = P_start; i != P_end; ++i) {
		j = A.rowBegin () + i->first;
		k = A.rowBegin () + i->second;

		_VD.swap (*j, *k);
	}

	return A;
}

template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGF2::permuteRowsByCol (Matrix   &A,
						  Iterator  P_start,
						  Iterator  P_end) const
{
	throw LinboxNotImplemented (__FUNCTION__);

	typename Matrix::ColIterator j;

//	for (j = A.colBegin (); j != A.colEnd (); ++j)
//		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGF2::permuteColsByRow (Matrix   &A,
						  Iterator  P_start,
						  Iterator  P_end) const
{
	typename Matrix::RowIterator j;

	for (j = A.rowBegin (); j != A.rowEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

}

#endif // __MATRIX_MATRIX_DOMAIN_GF2_INL

