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

#include "linbox/blas/level3-generic.h"
#include "linbox/blas/level1-ll.h"
#include "linbox/blas/level2-ll.h"
#include "linbox/blas/level3-ll.h"
#include "linbox/matrix/transpose.h"
#include "linbox/matrix/submatrix.h"

namespace LinBox
{

namespace BLAS3
{

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
Matrix2 &_copy<Ring, GenericModule::Tag>::copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::RowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		BLAS1::_copy<Ring, typename Modules::Tag>::op (F, M, *i, *j);

	return B;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
Matrix2 &_copy<Ring, GenericModule::Tag>::copy_impl (const Ring &F, Modules &M, const Matrix1 &A, Matrix2 &B,
						     MatrixIteratorTypes::Col, MatrixIteratorTypes::Col)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		BLAS1::_copy<Ring, typename Modules::Tag>::op (F, M, *i, *j);

	return B;
}

template <class Ring>
template <class Modules, class Matrix>
Matrix &_scal<Ring, GenericModule::Tag>::scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixIteratorTypes::Row)
{
	typename Matrix::RowIterator i;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i)
		BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, a, *i);

	return A;
}

template <class Ring>
template <class Modules, class Matrix>
Matrix &_scal<Ring, GenericModule::Tag>::scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Matrix &A, MatrixIteratorTypes::Col)
{
	typename Matrix::ColIterator i;

	for (i = A.colBegin (); i != A.colEnd (); ++i)
		BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, a, *i);

	return A;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
Matrix2 &_axpy<Ring, GenericModule::Tag>::axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::RowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, a, *i, *j);

	return B;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
Matrix2 &_axpy<Ring, GenericModule::Tag>::axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B,
						     MatrixIteratorTypes::Col, MatrixIteratorTypes::Col)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, a, *i, *j);

	return B;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm<Ring, GenericModule::Tag>::gemm_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i = A.rowBegin ();
	typename Matrix3::RowIterator j = C.rowBegin ();

	TransposeMatrix<const Matrix2> BT (B);

	for (; i != A.rowEnd (); ++i, ++j)
		BLAS2::_gemv<Ring, typename Modules::Tag>::op (F, M, a, BT, *i, b, *j);

	return C;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm<Ring, GenericModule::Tag>::gemm_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
						     MatrixIteratorTypes::Col, MatrixIteratorTypes::Col, MatrixIteratorTypes::Col)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstColIterator i = B.colBegin ();
	typename Matrix3::ColIterator j = C.colBegin ();

	for (; i != B.colEnd (); ++i, ++j)
		BLAS2::_gemv<Ring, typename Modules::Tag>::op (F, M, a, A, *i, b, *j);

	return C;
}

// FIXME: This assumes that the rows of C are dense!!!
template <class Ring>
template <class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm<Ring, GenericModule::Tag>::gemm_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Col, MatrixIteratorTypes::Row)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::RowIterator l1;
	typename Matrix3::Row::iterator l2;
	typename Ring::Element d;

	for (i = A.rowBegin (), l1 = C.rowBegin (); i != A.rowEnd (); ++i, ++l1) {
		for (j = B.colBegin (), l2 = l1->begin (); j != B.colEnd (); ++j, ++l2) {
			BLAS1::_dot<Ring, typename Modules::Tag>::op (F, M, d, *i, *j);
			F.mulin (*l2, b);
			F.axpyin (*l2, a, d);
		}
	}

	return C;
}

// FIXME: This assumes that the rows of C are dense!!!
template <class Ring>
template <class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm<Ring, GenericModule::Tag>::gemm_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Col, MatrixIteratorTypes::Col)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Matrix3::ColIterator l1;
	typename Matrix3::Col::iterator l2;
	typename Ring::Element d;

	for (j = B.colBegin (), l1 = C.colBegin (); j != B.colEnd (); ++j, ++l1) {
		for (i = A.rowBegin (), l2 = l1->begin (); i != A.rowEnd (); ++i, ++l2) {
			BLAS1::_dot<Ring, typename Modules::Tag>::op (F, M, d, *i, *j);
			F.mulin (*l2, b);
			F.axpyin (*l2, a, d);
		}
	}

	return C;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
Matrix2 &_trmm<Ring, GenericModule::Tag>::trmm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	if (A.rowdim () == 0)
		return B;
	else if (F.isZero (a))
		return _scal<Ring, typename Modules::Tag>::op (F, M, a, B);
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return _scal<Ring, typename Modules::Tag>::op (F, M, a, B);
		else {
			typename Ring::Element ai;

			if (!A.getEntry (ai, 0, 0))
				return _scal<Ring, typename Modules::Tag>::op (F, M, F.zero (), B);
			else {
				F.mulin (ai, a);
				return _scal<Ring, typename Modules::Tag>::op (F, M, ai, B);
			}
		}
	}
	else if ((((A.coldim () / 2) / Matrix1::colAlign) * Matrix1::colAlign) / Matrix1::rowAlign == 0) {
		size_t l = A.rowdim () / 2;
		typename Matrix1::ConstSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix1::ConstSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename Matrix2::AlignedSubmatrixType B1 (B, 0, 0, l, B.coldim ());
		typename Matrix2::AlignedSubmatrixType B2 (B, l, 0, B.rowdim () - l, B.coldim ());

		if (type == LowerTriangular) {
			typename Matrix1::ConstSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A22, B2, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, a, A21, B1, F.one (), B2);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A11, B1, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix1::ConstSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A11, B1, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, a, A12, B2, F.one (), B1);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A22, B2, type, diagIsOne);
		}

		return B;
	} else {
		size_t l = (((((A.rowdim () / 2) / Matrix1::colAlign) * Matrix1::colAlign) / Matrix1::rowAlign) * Matrix1::rowAlign);
		typename Matrix1::ConstAlignedSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix1::ConstAlignedSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename Matrix2::AlignedSubmatrixType B1 (B, 0, 0, l, B.coldim ());
		typename Matrix2::AlignedSubmatrixType B2 (B, l, 0, B.rowdim () - l, B.coldim ());

		if (type == LowerTriangular) {
			typename Matrix1::ConstAlignedSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A22, B2, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, a, A21, B1, F.one (), B2);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A11, B1, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix1::ConstAlignedSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A11, B1, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, a, A12, B2, F.one (), B1);
			_trmm<Ring, typename Modules::Tag>::op (F, M, a, A22, B2, type, diagIsOne);
		}

		return B;
	}
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
Matrix2 &_trsm<Ring, GenericModule::Tag>::trsm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne,
						     MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == B.rowdim ());

	if (A.rowdim () == 0)
		return B;
	else if (F.isZero (a))
		return _scal<Ring, typename Modules::Tag>::op (F, M, a, B);
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return _scal<Ring, typename Modules::Tag>::op (F, M, a, B);
		else {
			typename Ring::Element ai;

			if (!A.getEntry (ai, 0, 0))
				// FIXME: Should return an error
				return _scal<Ring, typename Modules::Tag>::op (F, M, F.zero (), B);
			else {
				F.invin (ai);
				F.mulin (ai, a);
				return _scal<Ring, typename Modules::Tag>::op (F, M, ai, B);
			}
		}
	}
	else if ((((A.coldim () / 2) / Matrix1::colAlign) * Matrix1::colAlign) / Matrix1::rowAlign == 0) {
		size_t l = A.rowdim () / 2;
		typename Matrix1::ConstSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix1::ConstSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename Matrix2::AlignedSubmatrixType B1 (B, 0, 0, l, B.coldim ());
		typename Matrix2::AlignedSubmatrixType B2 (B, l, 0, B.rowdim () - l, B.coldim ());

		if (type == LowerTriangular) {
			typename Matrix1::ConstSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trsm<Ring, typename Modules::Tag>::op (F, M, a, A11, B1, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A21, B1, a, B2);
			_trsm<Ring, typename Modules::Tag>::op (F, M, F.one (), A22, B2, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix1::ConstSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trsm<Ring, typename Modules::Tag>::op (F, M, a, A22, B2, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A12, B2, a, B1);
			_trsm<Ring, typename Modules::Tag>::op (F, M, F.one (), A11, B1, type, diagIsOne);
		}

		return B;
	} else {
		size_t l = (((((A.rowdim () / 2) / Matrix1::colAlign) * Matrix1::colAlign) / Matrix1::rowAlign) * Matrix1::rowAlign);
		typename Matrix1::ConstAlignedSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix1::ConstAlignedSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename Matrix2::AlignedSubmatrixType B1 (B, 0, 0, l, B.coldim ());
		typename Matrix2::AlignedSubmatrixType B2 (B, l, 0, B.rowdim () - l, B.coldim ());

		if (type == LowerTriangular) {
			typename Matrix1::ConstAlignedSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trsm<Ring, typename Modules::Tag>::op (F, M, a, A11, B1, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A21, B1, a, B2);
			_trsm<Ring, typename Modules::Tag>::op (F, M, F.one (), A22, B2, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix1::ConstAlignedSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trsm<Ring, typename Modules::Tag>::op (F, M, a, A22, B2, type, diagIsOne);
			_gemm<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A12, B2, a, B1);
			_trsm<Ring, typename Modules::Tag>::op (F, M, F.one (), A11, B1, type, diagIsOne);
		}

		return B;
	}
}

template <class Ring>
template <class Modules, class Iterator, class Matrix>
Matrix &_permute_rows<Ring, GenericModule::Tag>::permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Row)
{
	typename Matrix::RowIterator j, k;

	for (; P_begin != P_end; ++P_begin) {
		j = A.rowBegin () + P_begin->first;
		k = A.rowBegin () + P_begin->second;

		BLAS1::_swap<Ring, typename Modules::Tag>::op (F, M, *j, *k);
	}

	return A;
}

template <class Ring>
template <class Modules, class Iterator, class Matrix>
Matrix &_permute_rows<Ring, GenericModule::Tag>::permute_rows_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Col)
{
	typename Matrix::ColIterator j;

	for (j = A.colBegin (); j != A.colEnd (); ++j)
		BLAS1::_permute<Ring, typename Modules::Tag>::op (F, M, P_begin, P_end, *j);

	return A;
}

template <class Ring>
template <class Modules, class Iterator, class Matrix>
Matrix &_permute_cols<Ring, GenericModule::Tag>::permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Row)
{
	typename Matrix::RowIterator j;

	for (j = A.rowBegin (); j != A.rowEnd (); ++j)
		BLAS1::_permute<Ring, typename Modules::Tag>::op (F, M, P_begin, P_end, *j);

	return A;
}

template <class Ring>
template <class Modules, class Iterator, class Matrix>
Matrix &_permute_cols<Ring, GenericModule::Tag>::permute_cols_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A, MatrixIteratorTypes::Col)
{
	typename Matrix::ColIterator j, k;

	for (; P_begin != P_end; ++P_begin) {
		j = A.colBegin () + P_begin->first;
		k = A.colBegin () + P_begin->second;

		BLAS1::_swap<Ring, typename Modules::Tag>::op (F, M, *j, *k);
	}

	return A;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
bool _equal<Ring, GenericModule::Tag>::equal_impl (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
						   MatrixIteratorTypes::Row, MatrixIteratorTypes::Row)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstRowIterator i;
	typename Matrix2::ConstRowIterator j;

	i = A.rowBegin ();
	j = B.rowBegin ();

	for (; i != A.rowEnd (); ++i, ++j)
		if (!BLAS1::_equal<Ring, typename Modules::Tag>::op (F, M, *i, *j))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Matrix1, class Matrix2>
bool _equal<Ring, GenericModule::Tag>::equal_impl (const Ring &F, Modules &M, const Matrix1 &A, const Matrix2 &B,
						   MatrixIteratorTypes::Col, MatrixIteratorTypes::Col)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ConstColIterator i;
	typename Matrix2::ConstColIterator j;

	i = A.colBegin ();
	j = B.colBegin ();

	for (; i != A.colEnd (); ++i, ++j)
		if (!BLAS1::_equal<Ring, typename Modules::Tag>::op (F, M, *i, *j))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Matrix>
bool _is_zero<Ring, GenericModule::Tag>::is_zero_impl (const Ring &F, Modules &M, const Matrix &A, MatrixIteratorTypes::Row)
{
	typename Matrix::ConstRowIterator i;

	i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		if (!BLAS1::_is_zero<Ring, typename Modules::Tag>::op (F, M, *i))
			return false;

	return true;
}

template <class Ring>
template <class Modules, class Matrix>
bool _is_zero<Ring, GenericModule::Tag>::is_zero_impl (const Ring &F, Modules &M, const Matrix &A, MatrixIteratorTypes::Col)
{
	typename Matrix::ConstColIterator i;

	i = A.colBegin ();

	for (; i != A.colEnd (); ++i)
		if (!BLAS1::_is_zero<Ring, typename Modules::Tag>::op (F, M, *i))
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
