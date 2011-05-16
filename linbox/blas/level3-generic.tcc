/* linbox/blas/level3-generic.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_GENERIC_TCC
#define __BLAS_LEVEL3_GENERIC_TCC

#include <algorithm>

#include "linbox/blas/level1.h"
#include "linbox/blas/level2.h"
#include "linbox/blas/level3-generic.h"
#include "linbox/matrix/transpose.h"

namespace LinBox
{

namespace BLAS3
{

template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::RowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		BLAS1::_copy (F, M, *i, *j);

	return B;
}

template <class Field, class Matrix1, class Matrix2>
Matrix2 &copy_impl (const Field &F, GenericModule &M, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		BLAS1::_copy (F, M, *i, *j);

	return B;
}

template <class Field, class Modules, class Matrix>
Matrix &scal_impl (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A, MatrixCategories::RowMatrixTag)
{
	typename Matrix::RowIterator i;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i)
		BLAS1::_scal (F, M, a, *i);

	return A;
}

template <class Field, class Modules, class Matrix>
Matrix &scal_impl (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A, MatrixCategories::ColMatrixTag)
{
	typename Matrix::ColIterator i;

	for (i = A.colBegin (); i != A.colEnd (); ++i)
		BLAS1::_scal (F, M, a, *i);

	return A;
}

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy_impl (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::RowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		BLAS1::_axpy (F, M, a, *i, *j);

	return B;
}

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy_impl (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		BLAS1::_axpy (F, M, a, *i, *j);

	return B;
}

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i = A.rowBegin ();
	typename Matrix3::RowIterator j = C.rowBegin ();

	TransposeMatrix<const Matrix2> BT (B);

	for (; i != A.rowEnd (); ++i, ++j)
		BLAS2::_gemv (F, M, a, BT, *i, b, *j);

	return C;
}

template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstColIterator i = B.colBegin ();
	typename Matrix3::ColIterator j = C.colBegin ();

	for (; i != B.colEnd (); ++i, ++j)
		BLAS2::_gemv (F, M, a, A, *i, b, *j);

	return C;
}

// FIXME: This assumes that the rows of C are dense!!!
template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::RowMatrixTag)
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
			BLAS1::_dot (F, M, d, *i, *j);
			F.mulin (*l2, b);
			F.axpyin (*l2, a, d);
		}
	}

	return C;
}

// FIXME: This assumes that the rows of C are dense!!!
template <class Field, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C,
		    MatrixCategories::RowMatrixTag, MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag)
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
			BLAS1::_dot (F, M, d, *i, *j);
			F.mulin (*l2, b);
			F.axpyin (*l2, a, d);
		}
	}

	return C;
}

template <class Field, class Matrix1, class Matrix2>
Matrix2 &trmm_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	if (F.isZero (a))
		return _scal (F, M, a, B);

	typename Field::Element ai;

	if (diagIsOne)
		F.assign (ai, a);

	TransposeMatrix<Matrix2> BT (B);

	if (type == LowerTriangular) {
		typename Matrix2::RowIterator i_B = B.rowEnd ();
		typename Matrix1::ConstRowIterator i_A = A.rowEnd ();

		size_t idx = B.rowdim ();

		do {
			--i_B; --i_A; --idx;

			if (!diagIsOne) {
				if (VectorWrapper::getEntry (*i_A, ai, idx))
					F.mulin (ai, a);
				else
					F.assign (ai, F.zero ());
			}

			BLAS2::_gemv (F, M, a, BT, *i_A, ai, *i_B, 0, idx);
		} while (i_B != B.rowBegin ());
	}
	else if (type == UpperTriangular) {
		typename Matrix2::RowIterator i_B;
		typename Matrix1::ConstRowIterator i_A;

		size_t idx = 0;

		for (i_B = B.rowBegin (), i_A = A.rowBegin (); i_B != B.rowEnd (); ++i_B, ++i_A, ++idx) {
			if (!diagIsOne) {
				if (VectorWrapper::getEntry (*i_A, ai, idx))
					F.mulin (ai, a);
				else
					F.assign (ai, F.zero ());
			}

			BLAS2::_gemv (F, M, a, BT, *i_A, ai, *i_B, idx + 1, B.rowdim ());
		}
	}

	return B;
}

template <class Field, class Matrix1, class Matrix2>
Matrix2 &trsm_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	typename Field::Element ai, ai_inv, neg_ai_inv;

	if (diagIsOne) {
		F.assign (neg_ai_inv, F.minusOne ());
		F.assign (ai_inv, a);
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

				F.inv (ai_inv, ai);
				F.neg (neg_ai_inv, ai_inv);
				F.mulin (ai_inv, a);
			}

			BLAS2::_gemv (F, M, neg_ai_inv, BT, *i_A, ai_inv, *i_B, 0, idx);
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

				F.inv (ai_inv, ai);
				F.neg (neg_ai_inv, ai_inv);
				F.mulin (ai_inv, a);
			}

			BLAS2::_gemv (F, M, neg_ai_inv, BT, *i_A, ai_inv, *i_B, idx + 1, A.coldim ());
		} while (i_A != A.rowBegin ());
	}

	return B;
}

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_rows_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::RowMatrixTag)
{
	typename Matrix::RowIterator j, k;

	for (; P_begin != P_end; ++P_begin) {
		j = A.rowBegin () + P_begin->first;
		k = A.rowBegin () + P_begin->second;

		BLAS1::_swap (F, M, *j, *k);
	}

	return A;
}

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_rows_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::ColMatrixTag)
{
	typename Matrix::ColIterator j;

	for (j = A.colBegin (); j != A.colEnd (); ++j)
		BLAS1::_permute (F, M, P_begin, P_end, *j);

	return A;
}

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_cols_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::RowMatrixTag)
{
	typename Matrix::RowIterator j;

	for (j = A.rowBegin (); j != A.rowEnd (); ++j)
		BLAS1::_permute (F, M, P_begin, P_end, *j);

	return A;
}

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_cols_impl (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixCategories::ColMatrixTag)
{
	typename Matrix::ColIterator j, k;

	for (; P_begin != P_end; ++P_begin) {
		j = A.colBegin () + P_begin->first;
		k = A.colBegin () + P_begin->second;

		BLAS1::_swap (F, M, *j, *k);
	}

	return A;
}

template <class Field, class Matrix1, class Matrix2>
bool equal_impl (const Field &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		if (!BLAS1::_equal (F, M, *i, *j))
			return false;

	return true;
}

template <class Field, class Matrix1, class Matrix2>
bool equal_impl (const Field &F, GenericModule &M, const Matrix1 &A, const Matrix2 &B,
		 MatrixCategories::ColMatrixTag, MatrixCategories::ColMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		if (!BLAS1::_equal (F, M, *i, *j))
			return false;

	return true;
}

template <class Field, class Matrix>
bool is_zero_impl (const Field &F, GenericModule &M, const Matrix &A, MatrixCategories::RowMatrixTag)
{
	typename Matrix::ConstRowIterator i;

	i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		if (!BLAS1::_is_zero (F, M, *i))
			return false;

	return true;
}

template <class Field, class Matrix>
bool is_zero_impl (const Field &F, GenericModule &M, const Matrix &A, MatrixCategories::ColMatrixTag)
{
	typename Matrix::ConstColIterator i;

	i = A.colBegin ();

	for (; i != A.colEnd (); ++i)
		if (!BLAS1::_is_zero (F, M, *i))
			return false;

	return true;
}

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_GENERIC_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
