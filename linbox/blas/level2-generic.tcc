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

#include "linbox/blas/level2-generic.h"
#include "linbox/blas/level1-ll.h"
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

namespace BLAS2
{

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
						     size_t start_idx, size_t end_idx,
						     MatrixIteratorTypes::Row,
						     VectorRepresentationTypes::Generic,
						     VectorRepresentationTypes::Dense)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	linbox_check (end_idx <= A.coldim ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector2::iterator j = y.begin ();

	typename Ring::Element d;

	for (; i != A.rowEnd (); ++j, ++i) {
		BLAS1::_dot<Ring, typename Modules::Tag>::op (F, M, d, x, *i, start_idx, end_idx);
		F.mulin (*j, b);
		F.axpyin (*j, a, d);
	}

	return y;
}

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
						     size_t start_idx, size_t end_idx,
						     MatrixIteratorTypes::Row,
						     VectorRepresentationTypes::Generic,
						     VectorRepresentationTypes::Sparse)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	linbox_check (end_idx <= A.coldim ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Ring::Element t;
	unsigned int idx = 0;

	typename Vector<Ring>::Sparse yp;

	if (F.isZero (b))
		y.clear ();
	else
		BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, y, b);

	for (; i != A.rowEnd (); ++i, ++idx) {
		BLAS1::_dot<Ring, typename Modules::Tag>::op (F, M, t, x, *i, start_idx, end_idx);
		F.mulin (t, a);

		if (!F.isZero (t))
			yp.push_back (typename Vector<Ring>::Sparse::value_type (idx, t));
	}

	return BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, F.one (), yp, y);
}

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
						     size_t start_idx, size_t end_idx,
						     MatrixIteratorTypes::Col,
						     VectorRepresentationTypes::Dense,
						     VectorRepresentationTypes::Generic)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));
	linbox_check (start_idx <= end_idx);

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	// linbox_check (end_idx <= A.coldim ());

	typename Matrix::ConstColIterator i = A.colBegin () + start_idx;
	typename Vector1::const_iterator j = x.begin () + start_idx, j_end = x.begin () + end_idx;
	typename Ring::Element d;

	BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, b, y);

	for (; j != j_end; ++j, ++i) {
		F.mul (d, a, *j);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, d, *i, y);
	}

	return y;
}

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl (const Ring &F, Modules &M,
						     const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
						     size_t start_idx, size_t end_idx,
						     MatrixIteratorTypes::Col,
						     VectorRepresentationTypes::Sparse,
						     VectorRepresentationTypes::Generic)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	if (end_idx == (size_t) -1)
		end_idx = A.coldim ();

	linbox_check (end_idx <= A.coldim ());

	typename Vector1::const_iterator j =
		(start_idx == 0) ? x.begin () : std::lower_bound (x.begin (), x.end (), start_idx, VectorUtils::CompareSparseEntries ());
	typename Ring::Element d;

	BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, b, y);

	for (; j != x.end () && j->first < end_idx; ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		F.mul (d, a, j->second);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, d, *i, y);
	}

	return y;
}

template <class Ring>
template <class Modules, class Matrix, class Vector>
Vector &_trmv<Ring, typename GenericModule<Ring>::Tag>::op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));

	if (A.rowdim () == 0)
		return x;
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return x;
		else {
			typename Ring::Element ai;

			if (!A.getEntry (ai, 0, 0))
				// FIXME: Should return an error
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, F.zero (), x);
			else
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, ai, x);
		}
	}
	else if ((((A.coldim () / 2) / Matrix::colAlign) * Matrix::colAlign) / Matrix::rowAlign == 0) {
		size_t l = A.rowdim () / 2;
		typename Matrix::ConstSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix::ConstSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x1 (x, 0, l);
		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x2 (x, l, A.rowdim ());

		if (type == LowerTriangular) {
			typename Matrix::ConstSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trmv<Ring, typename Modules::Tag>::op (F, M, F.one (), A22, x2, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.one (), A21, x1, F.one (), x2);
			_trmv<Ring, typename Modules::Tag>::op (F, M, F.one (), A11, x1, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix::ConstSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trmv<Ring, typename Modules::Tag>::op (F, M, F.one (), A11, x1, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.one (), A12, x2, F.one (), x1);
			_trmv<Ring, typename Modules::Tag>::op (F, M, F.one (), A22, x2, type, diagIsOne);
		}

		return x;
	} else {
		size_t l = (((((A.rowdim () / 2) / Matrix::colAlign) * Matrix::colAlign) / Matrix::rowAlign) * Matrix::rowAlign);
		typename Matrix::ConstAlignedSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix::ConstAlignedSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x1 (x, 0, l);
		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x2 (x, l, A.rowdim ());

		if (type == LowerTriangular) {
			typename Matrix::ConstAlignedSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.one (), A21, x1, F.one (), x2);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix::ConstAlignedSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.one (), A12, x2, F.one (), x1);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
		}

		return x;
	}
}

template <class Ring>
template <class Modules, class Matrix, class Vector>
Vector &_trsv<Ring, typename GenericModule<Ring>::Tag>::op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
{
	linbox_check (A.coldim () == A.rowdim ());
	linbox_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));

	if (A.rowdim () == 0)
		return x;
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return x;
		else {
			typename Ring::Element ai;

			if (!A.getEntry (ai, 0, 0))
				// FIXME: Should return an error
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, F.zero (), x);
			else {
				F.invin (ai);
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, ai, x);
			}
		}
	}
	else if ((((A.coldim () / 2) / Matrix::colAlign) * Matrix::colAlign) / Matrix::rowAlign == 0) {
		size_t l = A.rowdim () / 2;
		typename Matrix::ConstSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix::ConstSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x1 (x, 0, l);
		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x2 (x, l, A.rowdim ());

		if (type == LowerTriangular) {
			typename Matrix::ConstSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A21, x1, F.one (), x2);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix::ConstSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A12, x2, F.one (), x1);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
		}

		return x;
	} else {
		size_t l = (((((A.rowdim () / 2) / Matrix::colAlign) * Matrix::colAlign) / Matrix::rowAlign) * Matrix::rowAlign);
		typename Matrix::ConstAlignedSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix::ConstAlignedSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x1 (x, 0, l);
		typename VectorTraits<Ring, Vector>::AlignedSubvectorType x2 (x, l, A.rowdim ());

		if (type == LowerTriangular) {
			typename Matrix::ConstAlignedSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A21, x1, F.one (), x2);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix::ConstAlignedSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.minusOne (), A12, x2, F.one (), x1);
			_trsv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
		}

		return x;
	}
}

template <class Ring>
template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
						  VectorRepresentationTypes::Dense,
						  VectorRepresentationTypes::Dense,
						  MatrixIteratorTypes::Row)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

	typename Matrix::RowIterator i_A;
	typename Vector1::const_iterator i_x;

	typename Ring::Element axi;

	for (i_A = A.rowBegin (), i_x = x.begin (); i_A != A.rowEnd (); ++i_A, ++i_x) {
		F.mul (axi, a, *i_x);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, axi, y, *i_A);
	}

	return A;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
						  VectorRepresentationTypes::Sparse,
						  VectorRepresentationTypes::Sparse,
						  MatrixIteratorTypes::Row)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

	typename Vector1::const_iterator i_x;

	typename Ring::Element axi;

	for (i_x = x.begin (); i_x != x.end (); ++i_x) {
		F.mul (axi, a, i_x->second);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, axi, y, *(A.rowBegin () + i_x->first));
	}

	return A;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
						  VectorRepresentationTypes::Dense,
						  VectorRepresentationTypes::Dense,
						  MatrixIteratorTypes::Col)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

	typename Matrix::ColIterator i_A;
	typename Vector2::const_iterator i_y;

	typename Ring::Element ayi;

	for (i_A = A.colBegin (), i_y = y.begin (); i_A != A.colEnd (); ++i_A, ++i_y) {
		F.mul (ayi, a, *i_y);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, ayi, x, *i_A);
	}

	return A;
}

template <class Ring>
template <class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
						  VectorRepresentationTypes::Sparse,
						  VectorRepresentationTypes::Sparse,
						  MatrixIteratorTypes::Col)
{
	linbox_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	linbox_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

	typename Vector2::const_iterator i_y;

	typename Ring::Element ayi;

	for (i_y = y.begin (); i_y != y.end (); ++i_y) {
		F.mul (ayi, a, i_y->second);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, ayi, x, *(A.colBegin () + i_y->first));
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
