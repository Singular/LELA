/* lela/blas/level2-generic.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 2 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_GENERIC_TCC
#define __BLAS_LEVEL2_GENERIC_TCC

#include <algorithm>

#include "lela/blas/level2-generic.h"
#include "lela/blas/level1-ll.h"
#include "lela/blas/level2-ll.h"
#include "lela/util/error.h"
#include "lela/integer.h"

namespace LELA
{

namespace BLAS2
{

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl
	(const Ring &F, Modules &M,
	 const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
	 MatrixIteratorTypes::Row,
	 VectorRepresentationTypes::Generic,
	 VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector2::iterator j = y.begin ();

	typename Ring::Element d;

	for (; i != A.rowEnd (); ++j, ++i) {
		BLAS1::_dot<Ring, typename Modules::Tag>::op (F, M, d, x, *i);
		F.mulin (*j, b);
		F.axpyin (*j, a, d);
	}

	return y;
}

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl
	(const Ring &F, Modules &M,
	 const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
	 MatrixIteratorTypes::Row,
	 VectorRepresentationTypes::Generic,
	 VectorRepresentationTypes::Sparse)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Ring::Element t;
	unsigned int idx = 0;

	typename Vector<Ring>::Sparse yp;

	BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, b, y);

	for (; i != A.rowEnd (); ++i, ++idx) {
		BLAS1::_dot<Ring, typename Modules::Tag>::op (F, M, t, x, *i);
		F.mulin (t, a);

		if (!F.isZero (t)) {
			yp.push_back (typename Vector<Ring>::Sparse::value_type (idx, typename Ring::Element ()));
			F.copy (yp.back ().second, t);
		}
	}

	return BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, F.one (), yp, y);
}

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl
	(const Ring &F, Modules &M,
	 const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
	 MatrixIteratorTypes::Col,
	 VectorRepresentationTypes::Dense,
	 VectorRepresentationTypes::Generic)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector1::const_iterator j = x.begin ();
	typename Ring::Element d;

	BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, b, y);

	for (; j != x.end (); ++j, ++i) {
		F.mul (d, a, *j);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, d, *i, y);
	}

	return y;
}

template <class Ring>
template <class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Ring, typename GenericModule<Ring>::Tag>::gemv_impl
	(const Ring &F, Modules &M,
	 const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
	 MatrixIteratorTypes::Col,
	 VectorRepresentationTypes::Sparse,
	 VectorRepresentationTypes::Generic)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.rowdim ()));

	typename Vector1::const_iterator j;
	typename Ring::Element d;

	BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, b, y);

	for (j = x.begin (); j != x.end (); ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + j->first;
		F.mul (d, a, j->second);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, d, *i, y);
	}

	return y;
}

template <class Ring>
template <class Modules, class Matrix, class Vector>
Vector &_trmv<Ring, typename GenericModule<Ring>::Tag>::trmv_impl
	(const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
	 VectorRepresentationTypes::Dense)
{
	lela_check (A.coldim () == A.rowdim ());
	lela_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));

	static const int align = const_lcm<const_lcm<Matrix::colAlign, Matrix::rowAlign>::val, VectorTraits<Ring, Vector>::align>::val;

	if (A.rowdim () == 0)
		return x;
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return x;
		else {
			typename Ring::Element ai;

			if (!A.getEntry (ai, 0, 0))
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, F.zero (), x);
			else
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, ai, x);
		}
	}
	else if (A.coldim () / (2 * align) == 0) {
		size_t l = A.rowdim () / 2;
		typename Matrix::ConstSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix::ConstSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename VectorTraits<Ring, Vector>::SubvectorType x1 (x, 0, l);
		typename VectorTraits<Ring, Vector>::SubvectorType x2 (x, l, A.rowdim ());

		if (type == LowerTriangular) {
			typename Matrix::ConstSubmatrixType A21 (A, l, 0, A.rowdim () - l, l);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.one (), A21, x1, F.one (), x2);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
		}
		else if (type == UpperTriangular) {
			typename Matrix::ConstSubmatrixType A12 (A, 0, l, l, A.coldim () - l);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A11, x1, type, diagIsOne);
			_gemv<Ring, typename Modules::Tag>::op (F, M, F.one (), A12, x2, F.one (), x1);
			_trmv<Ring, typename Modules::Tag>::op (F, M, A22, x2, type, diagIsOne);
		}

		return x;
	} else {
		size_t l = (A.rowdim () / (2 * align)) * align;
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
Vector &_trsv<Ring, typename GenericModule<Ring>::Tag>::trsv_impl
	(const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
	 VectorRepresentationTypes::Dense)
{
	lela_check (A.coldim () == A.rowdim ());
	lela_check (VectorUtils::hasDim<Ring> (x, A.coldim ()));

	static const int align = const_lcm<const_lcm<Matrix::colAlign, Matrix::rowAlign>::val, VectorTraits<Ring, Vector>::align>::val;

	if (A.rowdim () == 0)
		return x;
	else if (A.rowdim () == 1) {
		if (diagIsOne)
			return x;
		else {
			typename Ring::Element ai;

			if (!A.getEntry (ai, 0, 0) || !F.invin (ai))
				throw DiagonalEntryNotInvertible ();
			else
				return BLAS1::_scal<Ring, typename Modules::Tag>::op (F, M, ai, x);
		}
	}
	else if (A.coldim () / (2 * align) == 0) {
		size_t l = A.rowdim () / 2;
		typename Matrix::ConstSubmatrixType A11 (A, 0, 0, l, l);
		typename Matrix::ConstSubmatrixType A22 (A, l, l, A.rowdim () - l, A.coldim () - l);

		typename VectorTraits<Ring, Vector>::SubvectorType x1 (x, 0, l);
		typename VectorTraits<Ring, Vector>::SubvectorType x2 (x, l, A.rowdim ());

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
		size_t l = (A.rowdim () / (2 * align)) * align;
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
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
	 VectorRepresentationTypes::Dense,
	 VectorRepresentationTypes::Generic,
	 MatrixIteratorTypes::Row)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

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
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
	 VectorRepresentationTypes::Sparse,
	 VectorRepresentationTypes::Generic,
	 MatrixIteratorTypes::Row)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

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
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
	 VectorRepresentationTypes::Generic,
	 VectorRepresentationTypes::Dense,
	 MatrixIteratorTypes::Col)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

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
Matrix &_ger<Ring, typename GenericModule<Ring>::Tag>::ger_impl
	(const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
	 VectorRepresentationTypes::Generic,
	 VectorRepresentationTypes::Sparse,
	 MatrixIteratorTypes::Col)
{
	lela_check (VectorUtils::hasDim<Ring> (x, A.rowdim ()));
	lela_check (VectorUtils::hasDim<Ring> (y, A.coldim ()));

	typename Vector2::const_iterator i_y;

	typename Ring::Element ayi;

	for (i_y = y.begin (); i_y != y.end (); ++i_y) {
		F.mul (ayi, a, i_y->second);
		BLAS1::_axpy<Ring, typename Modules::Tag>::op (F, M, ayi, x, *(A.colBegin () + i_y->first));
	}

	return A;
}

} // namespace BLAS2

} // namespace LELA

#endif // __BLAS_LEVEL2_GENERIC_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
