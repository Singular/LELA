/* linbox/matrix/matrix-domain.inl
 * Copyright 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_domain_INL
#define __LINBOX_matrix_domain_INL

#include "linbox/matrix/transpose-matrix.h"

namespace LinBox
{

// FIXME: Currently ignores start_idx, end_idx!!!!
template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvRowSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::DenseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector2::iterator j = y.begin ();

	typename Field::Element d;

	for (; i != A.rowEnd (); ++j, ++i) {
		_VD.dot (d, x, *i);
		_F.mulin (*j, beta);
		_F.axpyin (*j, alpha, d);
	}

	return y;
}

// FIXME: Currently ignores start_idx, end_idx!!!!
template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvRowSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	typename Vector<Field>::Sparse yp;

	if (_F.isZero (beta))
		y.clear ();
	else
		_VD.mulin (y, beta);

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, x, *i);
		_F.mulin (t, alpha);

		if (!_F.isZero (t))
			yp.push_back (typename Vector<Field>::Sparse::value_type (idx, t));
	}

	return _VD.addin (y, yp);
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MVProductDomain<Field>::gemvColDense (const VectorDomain<Field>     &VD,
					       const typename Field::Element &alpha,
					       const Matrix                  &A,
					       const Vector1                 &x,
					       const typename Field::Element &beta,
					       Vector2                       &y,
					       size_t                         start_idx,
					       size_t                         end_idx) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector1::const_iterator j = x.begin () + start_idx, j_end = x.begin () + end_idx;
	typename Field::Element d;

	VD.mulin (y, beta);

	for (; j != j_end; ++j, ++i) {
		VD.field ().mul (d, alpha, *j);
		VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Vector1, class Matrix, class Vector2>
Vector2 &MatrixDomainSupportGeneric<Field>::gemvColSpecialized (const typename Field::Element &alpha,
								const Matrix                  &A,
								const Vector1                 &x,
								const typename Field::Element &beta,
								Vector2                       &y,
								size_t                         start_idx,
								size_t                         end_idx,
								VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	typename Vector1::const_iterator j = (start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Field::Element d;

	_VD.mulin (y, beta);

	for (; j != x.end () && j->first < end_idx; ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		_F.mul (d, alpha, j->second);
		_VD.axpyin (y, d, *i);
	}

	return y;
}

template <class Field>
template <class Vector1, class Vector2, class Matrix>
inline Matrix &MatrixDomainSupportGeneric<Field>::gerRowSpecialised (const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
								     VectorCategories::DenseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Matrix::RowIterator i_A;
	typename Vector1::const_iterator i_x;

	typename Field::Element axi;

	for (i_A = A.rowBegin (), i_x = x.begin (); i_A != A.rowEnd (); ++i_A, ++i_x) {
		_F.mul (axi, a, *i_x);
		_VD.axpyin (*i_A, axi, y);
	}

	return A;
}

template <class Field>
template <class Vector1, class Vector2, class Matrix>
inline Matrix &MatrixDomainSupportGeneric<Field>::gerRowSpecialised (const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
								     VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Vector1::const_iterator i_x;

	typename Field::Element axi;

	for (i_x = x.begin (); i_x != x.end (); ++i_x) {
		_F.mul (axi, a, i_x->second);
		_VD.axpyin (*(A.rowBegin () + i_x->first), axi, y);
	}

	return A;
}

template <class Field>
template <class Vector1, class Vector2, class Matrix>
inline Matrix &MatrixDomainSupportGeneric<Field>::gerColSpecialised (const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
								     VectorCategories::DenseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Matrix::ColIterator i_A;
	typename Vector2::const_iterator i_y;

	typename Field::Element ayi;

	for (i_A = A.colBegin (), i_y = y.begin (); i_A != A.colEnd (); ++i_A, ++i_y) {
		_F.mul (ayi, a, *i_y);
		_VD.axpyin (*i_A, ayi, x);
	}

	return A;
}

template <class Field>
template <class Vector1, class Vector2, class Matrix>
inline Matrix &MatrixDomainSupportGeneric<Field>::gerColSpecialised (const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
								     VectorCategories::SparseVectorTag) const
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Vector2::const_iterator i_y;

	typename Field::Element ayi;

	for (i_y = y.begin (); i_y != y.end (); ++i_y) {
		_F.mul (ayi, a, i_y->second);
		_VD.axpyin (*(A.colBegin () + i_y->first), ayi, x);
	}

	return A;
}

template <class Field>
template <class Matrix, class Vector>
Vector &MatrixDomainSupportGeneric<Field>::trsvSpecialized (const Matrix &A, Vector &x,
							    TriangularMatrixType type, bool diagIsOne,
							    MatrixCategories::RowMatrixTag,
							    VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));

	typename Field::Element ai, ai_p_1, neg_ai_inv, d;

	if (diagIsOne) {
		_F.assign (neg_ai_inv, _F.minusOne ());
		_F.init (ai_p_1, 2);
	}

	if (type == LowerTriangular) {
		typename Matrix::ConstRowIterator i_A;
		size_t idx = 0;

		for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A, ++idx) {
			if (!diagIsOne) {
				if (!VectorWrapper::getEntry (*i_A, ai, idx))
					continue; // FIXME This should throw an error

				_F.invin (ai);
				_F.add (ai_p_1, ai, _F.one ());
				_F.neg (neg_ai_inv, ai);
			}

			_VD.dot (d, *i_A, x);

			_F.mulin (x[idx], ai_p_1);
			_F.axpyin (x[idx], neg_ai_inv, d);
		}
	}
	else if (type == UpperTriangular) {
		typename Matrix::ConstRowIterator i_A = A.rowEnd ();
		size_t idx = A.rowdim ();

		do {
			--i_A; --idx;

			if (!diagIsOne) {
				if (!VectorWrapper::getEntry (*i_A, ai, idx))
					continue; // FIXME This should throw an error

				_F.invin (ai);
				_F.add (ai_p_1, ai, _F.one ());
				_F.neg (neg_ai_inv, ai);
			}

			_VD.dot (d, *i_A, x);

			_F.mulin (x[idx], ai_p_1);
			_F.axpyin (x[idx], neg_ai_inv, d);
		} while (i_A != A.rowBegin ());
	}

	return x;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomainSupportGeneric<Field>::copyRow (Matrix1 &A, const Matrix2 &B) const
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

template <class Field>
template <class Matrix1, class Matrix2>
Matrix1 &MatrixDomainSupportGeneric<Field>::copyCol (Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.copy (*i, *j);

	return A;
}

template <class Field>
template <class Matrix1, class Matrix2>
bool MatrixDomainSupportGeneric<Field>::areEqualRow (const Matrix1 &A, const Matrix2 &B) const
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

template <class Field>
template <class Matrix1, class Matrix2>
bool MatrixDomainSupportGeneric<Field>::areEqualCol (const Matrix1 &A, const Matrix2 &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		if (!_VD.areEqual (*i, *j))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
bool MatrixDomainSupportGeneric<Field>::isZeroRow (const Matrix &A) const
{
	typename Matrix::ConstRowIterator i;

	i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
bool MatrixDomainSupportGeneric<Field>::isZeroCol (const Matrix &A) const
{
	typename Matrix::ConstColIterator i;

	i = A.colBegin ();

	for (; i != A.colEnd (); ++i)
		if (!_VD.isZero (*i))
			return false;

	return true;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomainSupportGeneric<Field>::scalRow (Matrix &B, const typename Field::Element &a) const
{
	typename Matrix::RowIterator i;

	for (i = B.rowBegin (); i != B.rowEnd (); ++i)
		_VD.mulin (*i, a);

	return B;
}

template <class Field>
template <class Matrix>
Matrix &MatrixDomainSupportGeneric<Field>::scalCol (Matrix &B, const typename Field::Element &a) const
{
	typename Matrix::ColIterator i;

	for (i = B.colBegin (); i != B.colEnd (); ++i)
		_VD.mulin (*i, a);

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::axpyRow (const typename Field::Element &a,
						     const Matrix1                 &A,
						     Matrix2                       &B) const
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

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::axpyCol (const typename Field::Element &a,
						     const Matrix1                 &A,
						     Matrix2                       &B) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		_VD.axpyin (*j, a, *i);

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmRowRowCol (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::RowIterator l1;
	typename Matrix3::Row::iterator l2;
	typename Field::Element d;

	for (i = A.rowBegin (), l1 = C.rowBegin (); i != A.rowEnd (); ++i, ++l1) {
		for (j = B.colBegin (), l2 = l1->begin (); j != B.colEnd (); ++j, ++l2) {
			_VD.dot (d, *i, *j);
			_F.mulin (*l2, beta);
			_F.axpyin (*l2, alpha, d);
		}
	}

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmColRowCol (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::ColIterator l1;
	typename Matrix3::Col::iterator l2;
	typename Field::Element d;

	for (j = B.colBegin (), l1 = C.colBegin (); j != B.colEnd (); ++j, ++l1) {
		for (i = A.rowBegin (), l2 = l1->begin (); i != A.rowEnd (); ++i, ++l2) {
			_VD.dot (d, *i, *j);
			_F.mulin (*l2, beta);
			_F.axpyin (*l2, alpha, d);
		}
	}

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmRowRowRow (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i = A.rowBegin ();
	typename Matrix3::RowIterator j = C.rowBegin ();

	TransposeMatrix<const Matrix2> BT (B);

	for (; i != A.rowEnd (); ++i, ++j)
		gemv (alpha, BT, *i, beta, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix3 &MatrixDomainSupportGeneric<Field>::gemmColColCol (const typename Field::Element &alpha,
							   const Matrix1                 &A,
							   const Matrix2                 &B,
							   const typename Field::Element &beta,
							   Matrix3                       &C) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstColIterator i = B.colBegin ();
	typename Matrix3::ColIterator j = C.colBegin ();

	for (; i != B.colEnd (); ++i, ++j)
		gemv (alpha, A, *i, beta, *j);

	return C;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::trmmSpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
							     TriangularMatrixType type, bool diagIsOne,
							     MatrixCategories::RowMatrixTag,
							     MatrixCategories::RowMatrixTag) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	if (_F.isZero (a))
		return scal (B, a);

	typename Field::Element ai;

	if (diagIsOne)
		_F.assign (ai, a);

	TransposeMatrix<Matrix2> BT (B);

	if (type == LowerTriangular) {
		typename Matrix2::RowIterator i_B = B.rowEnd ();
		typename Matrix1::ConstRowIterator i_A = A.rowEnd ();

		size_t idx = B.rowdim ();

		do {
			--i_B; --i_A; --idx;

			if (!diagIsOne) {
				if (VectorWrapper::getEntry (*i_A, ai, idx))
					_F.mulin (ai, a);
				else
					_F.init (ai, 0);
			}

			gemvSpecialized (a, BT, *i_A, ai, *i_B, 0, idx, MatrixCategories::ColMatrixTag ());
		} while (i_B != B.rowBegin ());
	}
	else if (type == UpperTriangular) {
		typename Matrix2::RowIterator i_B;
		typename Matrix1::ConstRowIterator i_A;

		size_t idx = 0;

		for (i_B = B.rowBegin (), i_A = A.rowBegin (); i_B != B.rowEnd (); ++i_B, ++i_A, ++idx) {
			if (!diagIsOne) {
				if (VectorWrapper::getEntry (*i_A, ai, idx))
					_F.mulin (ai, a);
				else
					_F.init (ai, 0);
			}

			gemvSpecialized (a, BT, *i_A, ai, *i_B, idx + 1, B.rowdim (), MatrixCategories::ColMatrixTag ());
		}
	}

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::trsmSpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
							     TriangularMatrixType type, bool diagIsOne,
							     MatrixCategories::RowMatrixTag,
							     MatrixCategories::RowMatrixTag) const
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	typename Field::Element ai, ai_inv, neg_ai_inv;

	if (diagIsOne) {
		_F.assign (neg_ai_inv, _F.minusOne ());
		_F.assign (ai_inv, a);
	}

	TransposeMatrix<Matrix2> BT (B);

	if (type == LowerTriangular) {
		typename Matrix1::ConstRowIterator i_A;
		typename Matrix2::RowIterator i_B;
		size_t idx = 0;

		for (i_A = A.rowBegin (), i_B = B.rowBegin (); i_A != A.rowEnd (); ++i_A, ++i_B, ++idx) {
			if (!diagIsOne) {
				if (!VectorWrapper::getEntry (*i_A, ai, idx))
					continue; // FIXME This should throw an error

				_F.inv (ai_inv, ai);
				_F.neg (neg_ai_inv, ai_inv);
				_F.mulin (ai_inv, a);
			}

			gemvSpecialized (neg_ai_inv, BT, *i_A, ai_inv, *i_B, 0, idx, MatrixCategories::ColMatrixTag ());
		}
	}
	else if (type == UpperTriangular) {
		typename Matrix1::ConstRowIterator i_A = A.rowEnd ();
		typename Matrix2::RowIterator i_B = B.rowEnd ();
		size_t idx = A.rowdim ();

		do {
			--i_A; --i_B; --idx;

			if (!diagIsOne) {
				if (!VectorWrapper::getEntry (*i_A, ai, idx))
					continue; // FIXME This should throw an error

				_F.inv (ai_inv, ai);
				_F.neg (neg_ai_inv, ai_inv);
				_F.mulin (ai_inv, a);
			}

			gemvSpecialized (neg_ai_inv, BT, *i_A, ai_inv, *i_B, idx + 1, A.coldim (), MatrixCategories::ColMatrixTag ());
		} while (i_A != A.rowBegin ());
	}

	return B;
}

template <class Field>
template <class Matrix1, class Matrix2>
Matrix2 &MatrixDomainSupportGeneric<Field>::trsmSpecialized (const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
							     TriangularMatrixType type, bool diagIsOne,
							     MatrixCategories::RowMatrixTag,
							     MatrixCategories::ColMatrixTag) const
{
	typename Matrix2::ColIterator i_B;

	for (i_B = B.colBegin (); i_B != B.colEnd (); ++i_B) {
		trsv (A, *i_B, type, diagIsOne);
		_VD.mulin (*i_B, a);
	}

	return B;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteRowsByRow (Matrix   &A,
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

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteRowsByCol (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	typename Matrix::ColIterator j;

	for (j = A.colBegin (); j != A.colEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteColsByRow (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	typename Matrix::RowIterator j;

	for (j = A.rowBegin (); j != A.rowEnd (); ++j)
		_VD.permute (*j, P_start, P_end);

	return A;
}

template <class Field>
template <class Matrix, class Iterator>
Matrix &MatrixDomainSupportGeneric<Field>::permuteColsByCol (Matrix   &A,
							     Iterator  P_start,
							     Iterator  P_end) const
{
	Iterator i;
	typename Matrix::ColIterator j, k;

	for (i = P_start; i != P_end; ++i) {
		j = A.colBegin () + i->first;
		k = A.colBegin () + i->second;

		_VD.swap (*j, *k);
	}

	return A;
}

} // namespace LinBox

#endif // __LINBOX_matrix_domain_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
