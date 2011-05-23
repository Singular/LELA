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
#include "linbox/matrix/submatrix.h"

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

template <class Field, class Modules, class Matrix>
Matrix &_scal (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A);

template <class Field, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C);

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trmm (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne);

template <class Field, class Matrix1, class Matrix2>
Matrix2 &trmm_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	if (A.rowdim () == 0)
		return B;
	else if (F.isZero (a))
		return _scal (F, M, a, B);
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return _scal (F, M, a, B);
		else {
			typename Field::Element ai;

			if (!A.getEntry (ai, 0, 0))
				return _scal (F, M, F.zero (), B);
			else {
				F.mulin (ai, a);
				return _scal (F, M, ai, B);
			}
		}
	} else {
		size_t l = A.rowdim () / 2;
		Submatrix<const typename RealMatrixType<Matrix1>::Type> A11 (A, 0, 0, l, l);
		Submatrix<const typename RealMatrixType<Matrix1>::Type> A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		Submatrix<typename RealMatrixType<Matrix2>::Type> B1 (B, 0, 0, l, B.coldim ());
		Submatrix<typename RealMatrixType<Matrix2>::Type> B2 (B, l, 0, B.rowdim () - l, B.coldim ());

		if (type == LowerTriangular) {
			Submatrix<const typename RealMatrixType<Matrix1>::Type> A21 (A, l, 0, A.rowdim () - l, l);
			_trmm (F, M, a, A22, B2, type, diagIsOne);
			_gemm (F, M, a, A21, B1, F.one (), B2);
			_trmm (F, M, a, A11, B1, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			Submatrix<const typename RealMatrixType<Matrix1>::Type> A12 (A, 0, l, l, A.coldim () - l);
			_trmm (F, M, a, A11, B1, type, diagIsOne);
			_gemm (F, M, a, A12, B2, F.one (), B1);
			_trmm (F, M, a, A22, B2, type, diagIsOne);
		}

		return B;
	}
}

template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trsm (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne);

template <class Field, class Matrix1, class Matrix2>
Matrix2 &trsm_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
		    MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	if (A.rowdim () == 0)
		return B;
	else if (F.isZero (a))
		return _scal (F, M, a, B);
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return _scal (F, M, a, B);
		else {
			typename Field::Element ai;

			if (!A.getEntry (ai, 0, 0))
				// FIXME: Should return an error
				return _scal (F, M, F.zero (), B);
			else {
				F.invin (ai);
				F.mulin (ai, a);
				return _scal (F, M, ai, B);
			}
		}
	} else {
		size_t l = A.rowdim () / 2;
		Submatrix<const typename RealMatrixType<Matrix1>::Type> A11 (A, 0, 0, l, l);
		Submatrix<const typename RealMatrixType<Matrix1>::Type> A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		Submatrix<typename RealMatrixType<Matrix2>::Type> B1 (B, 0, 0, l, B.coldim ());
		Submatrix<typename RealMatrixType<Matrix2>::Type> B2 (B, l, 0, B.rowdim () - l, B.coldim ());

		if (type == LowerTriangular) {
			Submatrix<const typename RealMatrixType<Matrix1>::Type> A21 (A, l, 0, A.rowdim () - l, l);
			_trsm (F, M, a, A11, B1, type, diagIsOne);
			_gemm (F, M, F.minusOne (), A21, B1, a, B2);
			_trsm (F, M, F.one (), A22, B2, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			Submatrix<const typename RealMatrixType<Matrix1>::Type> A12 (A, 0, l, l, A.coldim () - l);
			_trsm (F, M, a, A22, B2, type, diagIsOne);
			_gemm (F, M, F.minusOne (), A12, B2, a, B1);
			_trsm (F, M, F.one (), A11, B1, type, diagIsOne);
		}

		return B;
	}
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
