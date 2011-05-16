/* linbox/blas/level2-generic.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 2 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_GENERIC_TCC
#define __BLAS_LEVEL2_GENERIC_TCC

#include <algorithm>

#include "linbox/blas/level1.h"
#include "linbox/blas/level2-generic.h"

namespace LinBox
{

namespace BLAS2
{

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::DenseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector2::iterator j = y.begin ();

	typename Field::Element d;

	for (; i != A.rowEnd (); ++j, ++i) {
		BLAS1::_dot (F, M, d, x, *i, start_idx, end_idx);
		F.mulin (*j, b);
		F.axpyin (*j, a, d);
	}

	return y;
}

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::RowMatrixTag,
		    VectorCategories::GenericVectorTag,
		    VectorCategories::SparseVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Field::Element t;
	unsigned int idx = 0;

	typename Vector<Field>::Sparse yp;

	if (F.isZero (b))
		y.clear ();
	else
		BLAS1::_scal (F, M, y, b);

	for (; i != A.rowEnd (); ++i, ++idx) {
		BLAS1::_dot (F, M, t, x, *i, start_idx, end_idx);
		F.mulin (t, a);

		if (!F.isZero (t))
			yp.push_back (typename Vector<Field>::Sparse::value_type (idx, t));
	}

	return BLAS1::_axpy (F, M, F.one (), yp, y);
}

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::DenseVectorTag,
		    VectorCategories::GenericVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector1::const_iterator j = x.begin () + start_idx, j_end = x.begin () + end_idx;
	typename Field::Element d;

	BLAS1::_scal (F, M, b, y);

	for (; j != j_end; ++j, ++i) {
		F.mul (d, a, *j);
		BLAS1::_axpy (F, M, d, *i, y);
	}

	return y;
}

template <class Field, class Matrix, class Vector1, class Vector2>
Vector2 &gemv_impl (const Field &F, GenericModule &M,
		    const typename Field::Element &a, const Matrix &A, const Vector1 &x, const typename Field::Element &b, Vector2 &y,
		    size_t start_idx, size_t end_idx,
		    MatrixCategories::ColMatrixTag,
		    VectorCategories::SparseVectorTag,
		    VectorCategories::GenericVectorTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.coldim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.rowdim ()));

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	typename Vector1::const_iterator j =
		(start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorWrapper::CompareSparseEntries ());
	typename Field::Element d;

	BLAS1::_scal (F, M, b, y);

	for (; j != x.end () && j->first < end_idx; ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		F.mul (d, a, j->second);
		BLAS1::_axpy (F, M, d, *i, y);
	}

	return y;
}

// FIXME: Not yet implemented
template <class Field, class Matrix, class Vector>
Vector &trmv_impl (const Field &F, GenericModule &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowMatrixTag,
		   VectorCategories::DenseVectorTag);

template <class Field, class Matrix, class Vector>
Vector &trsv_impl (const Field &F, GenericModule &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
		   MatrixCategories::RowMatrixTag,
		   VectorCategories::DenseVectorTag)
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));

	typename Field::Element ai, ai_inv, neg_ai_inv, d;

	if (diagIsOne) {
		F.assign (neg_ai_inv, F.minusOne ());
	}

	if (type == LowerTriangular) {
		typename Matrix::ConstRowIterator i_A;
		size_t idx = 0;

		for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A, ++idx) {
			if (!diagIsOne) {
				if (!VectorWrapper::getEntry (*i_A, ai, idx))
					continue; // FIXME This should throw an error

				F.inv (ai_inv, ai);
				F.neg (neg_ai_inv, ai_inv);
			}

			BLAS1::_dot (F, M, d, *i_A, x, 0, idx);

			F.mulin (x[idx], ai_inv);
			F.axpyin (x[idx], d, neg_ai_inv);
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

				F.inv (ai_inv, ai);
				F.neg (neg_ai_inv, ai_inv);
			}

			BLAS1::_dot (F, M, d, *i_A, x, idx + 1, (size_t) -1);

			F.mulin (x[idx], ai_inv);
			F.axpyin (x[idx], d, neg_ai_inv);
		} while (i_A != A.rowBegin ());
	}

	return x;
}

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseVectorTag,
		  VectorCategories::DenseVectorTag,
		  MatrixCategories::RowMatrixTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Matrix::RowIterator i_A;
	typename Vector1::const_iterator i_x;

	typename Field::Element axi;

	for (i_A = A.rowBegin (), i_x = x.begin (); i_A != A.rowEnd (); ++i_A, ++i_x) {
		F.mul (axi, a, *i_x);
		BLAS1::_axpy (F, M, axi, y, *i_A);
	}

	return A;
}

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseVectorTag,
		  VectorCategories::SparseVectorTag,
		  MatrixCategories::RowMatrixTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Vector1::const_iterator i_x;

	typename Field::Element axi;

	for (i_x = x.begin (); i_x != x.end (); ++i_x) {
		F.mul (axi, a, i_x->second);
		BLAS1::_axpy (F, M, axi, y, *(A.rowBegin () + i_x->first));
	}

	return A;
}

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::DenseVectorTag,
		  VectorCategories::DenseVectorTag,
		  MatrixCategories::ColMatrixTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Matrix::ColIterator i_A;
	typename Vector2::const_iterator i_y;

	typename Field::Element ayi;

	for (i_A = A.colBegin (), i_y = y.begin (); i_A != A.colEnd (); ++i_A, ++i_y) {
		F.mul (ayi, a, *i_y);
		BLAS1::_axpy (F, M, ayi, x, *i_A);
	}

	return A;
}

template <class Field, class Vector1, class Vector2, class Matrix>
Matrix &ger_impl (const Field &F, GenericModule &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
		  VectorCategories::SparseVectorTag,
		  VectorCategories::SparseVectorTag,
		  MatrixCategories::ColMatrixTag)
{
	linbox_check (VectorWrapper::hasDim<Field> (x, A.rowdim ()));
	linbox_check (VectorWrapper::hasDim<Field> (y, A.coldim ()));

	typename Vector2::const_iterator i_y;

	typename Field::Element ayi;

	for (i_y = y.begin (); i_y != y.end (); ++i_y) {
		F.mul (ayi, a, i_y->second);
		BLAS1::_axpy (F, M, ayi, x, *(A.colBegin () + i_y->first));
	}

	return A;
}

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_GENERIC_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
