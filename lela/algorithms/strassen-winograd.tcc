/* lela/algorithms/strassen-winograd.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementation of Straßen-Winograd fast matrix multiplication based on
 *
 * Boyer, B., Dumas, J.-G., Pernet, C., & Zhou, W. (2007). Memory efficient
 * scheduling of Strassen-Winogradʼs matrix multiplication algorithm.
 * Proceedings of the 2009 international symposium on Symbolic and algebraic
 * computation ISSAC 09, 55. ACM Press. Retrieved from
 * http://arxiv.org/abs/0707.2347
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ALGORITHMS_STRASSEN_WINOGRAD_TCC
#define __LELA_ALGORITHMS_STRASSEN_WINOGRAD_TCC

#include "lela/blas/level3-ll.h"
#include "lela/algorithms/strassen-winograd.h"
#include "lela/integer.h"
#include "lela/util/timer.h"

// Define to get a detailed profile of what the algorithm is doing at a cost of performance
#define __LELA_SW_DETAILED_PROFILE

#ifdef __LELA_SW_DETAILED_PROFILE
#  define SW_TIMER_DECLARE(part) LELA::UserTimer part##_timer; double part##_time = 0.0;
#  define SW_TIMER_START(part) part##_timer.start ()
#  define SW_TIMER_STOP(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define SW_TIMER_REPORT(part) \
	commentator.report (Commentator::LEVEL_NORMAL, TIMING_MEASURE) \
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
#else
#  define SW_TIMER_DECLARE(part)
#  define SW_TIMER_START(part)
#  define SW_TIMER_STOP(part)
#  define SW_TIMER_REPORT(part)
#endif

namespace LELA
{

template <class T1, class T2>
T1 round_down (T1 n, const T2 m)
	{ return m * (n / m); }

template <class Matrix1, class Matrix2>
inline size_t align_row (size_t i)
	{ return round_down (i, const_lcm <Matrix1::rowAlign, Matrix2::rowAlign>::val); }

template <class Matrix1, class Matrix2>
inline size_t align_col (size_t i)
	{ return round_down (i, const_lcm <Matrix1::colAlign, Matrix2::colAlign>::val); }

template <class RowMatrix, class ColMatrix>
inline size_t align_rowcol (size_t i)
	{ return round_down (i, const_lcm <RowMatrix::rowAlign, ColMatrix::colAlign>::val); }

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::gemm_res (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
						size_t m, size_t k, size_t n)
{
#ifdef __LELA_SW_DETAILED_PROFILE
	commentator.start ("Residual gemm", __FUNCTION__);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "Size-parameters: " << m << ", " << k << ", " << n << std::endl;
#endif // __LELA_SW_DETAILED_PROFILE

	if (2 * k < A.coldim ()) {
		typename Matrix3::AlignedSubmatrixType      C11 (C, 0,     0,     2 * m,               2 * n);
		typename Matrix1::ConstAlignedSubmatrixType A12 (A, 0,     2 * k, 2 * m,               A.coldim () - 2 * k);
		typename Matrix2::ConstAlignedSubmatrixType B21 (B, 2 * k, 0,     B.rowdim () - 2 * k, 2 * n);

		BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A12, B21, R.one (), C11);
	}

	if (2 * n < C.coldim ()) {
		typename Matrix3::AlignedSubmatrixType      C12 (C, 0,     2 * n, 2 * m, C.coldim () - 2 * n);
		typename Matrix1::ConstAlignedSubmatrixType A11 (A, 0,     0,     2 * m, 2 * k);
		typename Matrix2::ConstAlignedSubmatrixType B12 (B, 0,     2 * n, 2 * k, B.coldim () - 2 * n);

		BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A11, B12, b, C12);

		if (2 * k < A.coldim ()) {
			typename Matrix1::ConstAlignedSubmatrixType A12 (A, 0,     2 * k, 2 * m,               A.coldim () - 2 * k);
			typename Matrix2::ConstAlignedSubmatrixType B22 (B, 2 * k, 2 * n, B.rowdim () - 2 * k, B.coldim () - 2 * n);

			BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A12, B22, R.one (), C12);
		}
	}

	if (2 * m < C.rowdim ()) {
		typename Matrix3::AlignedSubmatrixType      C21 (C, 2 * m, 0,     C.rowdim () - 2 * m, 2 * n);
		typename Matrix1::ConstAlignedSubmatrixType A21 (A, 2 * m, 0,     A.rowdim () - 2 * m, 2 * k);
		typename Matrix2::ConstAlignedSubmatrixType B11 (B, 0,     0,     2 * k,               2 * n);

		BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A21, B11, b, C21);

		if (2 * k < A.coldim ()) {
			typename Matrix1::ConstAlignedSubmatrixType A22 (A, 2 * m, 2 * k, A.rowdim () - 2 * m, A.coldim () - 2 * k);
			typename Matrix2::ConstAlignedSubmatrixType B21 (B, 2 * k, 0,     B.rowdim () - 2 * k, 2 * n);

			BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A22, B21, R.one (), C21);
		}
	}

	if (2 * m < A.rowdim () && 2 * n < B.coldim ()) {
		typename Matrix3::AlignedSubmatrixType      C22 (C, 2 * m, 2 * n, C.rowdim () - 2 * m, C.coldim () - 2 * n);
		typename Matrix1::ConstAlignedSubmatrixType A21 (A, 2 * m, 0,     A.rowdim () - 2 * m, 2 * k);
		typename Matrix2::ConstAlignedSubmatrixType B12 (B, 0,     2 * n, 2 * k,               B.coldim () - 2 * n);

		BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A21, B12, b, C22);

		if (2 * k < A.coldim ()) {
			typename Matrix1::ConstAlignedSubmatrixType A22 (A, 2 * m, 2 * k, A.rowdim () - 2 * m, A.coldim () - 2 * k);
			typename Matrix2::ConstAlignedSubmatrixType B22 (B, 2 * k, 2 * n, B.rowdim () - 2 * k, B.coldim () - 2 * n);

			BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A22, B22, R.one (), C22);
		}
	}

#ifdef __LELA_SW_DETAILED_PROFILE
	commentator.stop (MSG_DONE);
#endif // __LELA_SW_DETAILED_PROFILE

	return C;	
}

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::mul (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, Matrix3 &C)
{
	lela_check (A.coldim () == B.rowdim ());
	lela_check (A.rowdim () == C.rowdim ());
	lela_check (B.coldim () == C.coldim ());

	size_t m = align_row<Matrix1, Matrix3> (C.rowdim () / 2), k = align_rowcol<Matrix1, Matrix2> (A.coldim () / 2), n = align_col<Matrix2, Matrix3> (C.coldim () / 2);

	if (C.rowdim () < _cutoff || C.coldim () < _cutoff || A.coldim () < _cutoff || m == 0 || n == 0 || k == 0)
		return BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A, B, R.zero (), C);
	else {
#ifdef __LELA_SW_DETAILED_PROFILE
		commentator.start ("StrassenWinograd::mul", __FUNCTION__);

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Sizes: " << C.rowdim () << ", " << A.coldim () << ", " << C.coldim () << std::endl;

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Size-parameters: " << m << ", " << k << ", " << n << std::endl;

#endif // __LELA_SW_DETAILED_PROFILE

		SW_TIMER_DECLARE(other);
		SW_TIMER_DECLARE(mul1);
		SW_TIMER_DECLARE(mul2);
		SW_TIMER_DECLARE(mul3);
		SW_TIMER_DECLARE(mul4);
		SW_TIMER_DECLARE(mul5);
		SW_TIMER_DECLARE(mul6);
		SW_TIMER_DECLARE(mul7);

		SW_TIMER_START(other);

		typename Matrix3::ContainerType X1 (m, std::max (k, n)), X2 (k, n);

		typename Matrix1::ConstAlignedSubmatrixType A11 (A, 0, 0, m, k);
		typename Matrix1::ConstAlignedSubmatrixType A12 (A, 0, k, m, k);
		typename Matrix1::ConstAlignedSubmatrixType A21 (A, m, 0, m, k);
		typename Matrix1::ConstAlignedSubmatrixType A22 (A, m, k, m, k);

		typename Matrix2::ConstAlignedSubmatrixType B11 (B, 0, 0, k, n);
		typename Matrix2::ConstAlignedSubmatrixType B12 (B, 0, n, k, n);
		typename Matrix2::ConstAlignedSubmatrixType B21 (B, k, 0, k, n);
		typename Matrix2::ConstAlignedSubmatrixType B22 (B, k, n, k, n);

		typename Matrix3::AlignedSubmatrixType C11 (C, 0, 0, m, n);
		typename Matrix3::AlignedSubmatrixType C12 (C, 0, n, m, n);
		typename Matrix3::AlignedSubmatrixType C21 (C, m, 0, m, n);
		typename Matrix3::AlignedSubmatrixType C22 (C, m, n, m, n);

		typename Matrix3::ContainerType::AlignedSubmatrixType X11 (X1, 0, 0, m, k);
		typename Matrix3::ContainerType::AlignedSubmatrixType X12 (X1, 0, 0, m, n);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, A11, X11);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), A21, X11);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, B22, X2);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), B12, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul1);
		mul (R, M, R.one (), X11, X2, C21);
		SW_TIMER_STOP(mul1);

		SW_TIMER_START(other);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, A21, X11);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), A22, X11);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, B12, X2);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), B11, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul2);
		mul (R, M, R.one (), X11, X2, C22);
		SW_TIMER_STOP(mul2);

		SW_TIMER_START(other);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), A11, X11);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, R.minusOne (), X2);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), B22, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul3);
		mul (R, M, R.one (), X11, X2, C12);
		SW_TIMER_STOP(mul3);

		SW_TIMER_START(other);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, R.minusOne (), X11);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), A12, X11);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul4);
		mul (R, M, R.one (), X11, B22, C11);
		SW_TIMER_STOP(mul4);
		SW_TIMER_START(mul5);
		mul (R, M, R.one (), A11, B11, X12);
		SW_TIMER_STOP(mul5);

		SW_TIMER_START(other);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X12, C12);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), C12, C21);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), C22, C12);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), C21, C22);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), C11, C12);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), B21, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul6);
		mul (R, M, R.one (), A22, X2, C11);
		SW_TIMER_STOP(mul6);

		SW_TIMER_START(other);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), C11, C21);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul7);
		mul (R, M, R.one (), A12, B21, C11);
		SW_TIMER_STOP(mul7);

		SW_TIMER_START(other);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X12, C11);

		typename Matrix3::AlignedSubmatrixType C_part (C, 0, 0, 2 * m, 2 * n);
		BLAS3::_scal<Ring, ParentTag>::op (R, M, a, C_part);

		SW_TIMER_STOP(other);

		SW_TIMER_REPORT(mul1);
		SW_TIMER_REPORT(mul2);
		SW_TIMER_REPORT(mul3);
		SW_TIMER_REPORT(mul4);
		SW_TIMER_REPORT(mul5);
		SW_TIMER_REPORT(mul6);
		SW_TIMER_REPORT(mul7);
		SW_TIMER_REPORT(other);

#ifdef __LELA_SW_DETAILED_PROFILE
		commentator.stop (MSG_DONE);
#endif // __LELA_SW_DETAILED_PROFILE

		return gemm_res (R, M, a, A, B, R.zero (), C, m, k, n);
	}
}

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::addmul (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
{
	lela_check (A.coldim () == B.rowdim ());
	lela_check (A.rowdim () == C.rowdim ());
	lela_check (B.coldim () == C.coldim ());

	size_t m = align_row<Matrix1, Matrix3> (C.rowdim () / 2), k = align_rowcol<Matrix1, Matrix2> (A.coldim () / 2), n = align_col<Matrix2, Matrix3> (C.coldim () / 2);

	if (C.rowdim () < _cutoff || C.coldim () < _cutoff || A.coldim () < _cutoff || m == 0 || n == 0 || k == 0)
		return BLAS3::_gemm<Ring, ParentTag>::op (R, M, a, A, B, b, C);
	else {
#ifdef __LELA_SW_DETAILED_PROFILE
		commentator.start ("StrassenWinograd::mul", __FUNCTION__);

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Sizes: " << C.rowdim () << ", " << A.coldim () << ", " << C.coldim () << std::endl;

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Size-parameters: " << m << ", " << k << ", " << n << std::endl;
#endif // __LELA_SW_DETAILED_PROFILE

		SW_TIMER_DECLARE(other);
		SW_TIMER_DECLARE(mul1);
		SW_TIMER_DECLARE(mul2);
		SW_TIMER_DECLARE(mul3);
		SW_TIMER_DECLARE(mul4);
		SW_TIMER_DECLARE(mul5);
		SW_TIMER_DECLARE(mul6);
		SW_TIMER_DECLARE(mul7);

		SW_TIMER_START(other);

		typename Matrix3::ContainerType X1 (m, k), X2 (k, n), X3 (m, n);

		typename Matrix1::ConstAlignedSubmatrixType A11 (A, 0, 0, m, k);
		typename Matrix1::ConstAlignedSubmatrixType A12 (A, 0, k, m, k);
		typename Matrix1::ConstAlignedSubmatrixType A21 (A, m, 0, m, k);
		typename Matrix1::ConstAlignedSubmatrixType A22 (A, m, k, m, k);

		typename Matrix2::ConstAlignedSubmatrixType B11 (B, 0, 0, k, n);
		typename Matrix2::ConstAlignedSubmatrixType B12 (B, 0, n, k, n);
		typename Matrix2::ConstAlignedSubmatrixType B21 (B, k, 0, k, n);
		typename Matrix2::ConstAlignedSubmatrixType B22 (B, k, n, k, n);

		typename Matrix3::AlignedSubmatrixType C11 (C, 0, 0, m, n);
		typename Matrix3::AlignedSubmatrixType C12 (C, 0, n, m, n);
		typename Matrix3::AlignedSubmatrixType C21 (C, m, 0, m, n);
		typename Matrix3::AlignedSubmatrixType C22 (C, m, n, m, n);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, A21, X1);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), A22, X1);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, B12, X2);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), B11, X2);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, R.zero (), X3);
		SW_TIMER_STOP(other);

		SW_TIMER_START(mul1);
		mul (R, M, a, X1, X2, X3);
		SW_TIMER_STOP(mul1);

		SW_TIMER_START(other);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, b, C22);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X3, C22);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, b, C12);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X3, C12);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), A11, X1);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, R.minusOne (), X2);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), B22, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul2);
		mul (R, M, a, A11, B11, X3);
		SW_TIMER_STOP(mul2);

		SW_TIMER_START(other);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, b, C11);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X3, C11);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul3);
		addmul (R, M, a, X1, X2, R.one (), X3);
		SW_TIMER_STOP(mul3);
		SW_TIMER_START(mul4);
		addmul (R, M, a, A12, B21, R.one (), C11);
		SW_TIMER_STOP(mul4);

		SW_TIMER_START(other);

		BLAS3::_scal<Ring, ParentTag>::op (R, M, R.minusOne (), X1);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), A12, X1);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), B21, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul5);
		addmul (R, M, a, X1, B22, R.one (), C12);
		SW_TIMER_STOP(mul5);

		SW_TIMER_START(other);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X3, C12);

		typename Ring::Element neg_b;
		R.neg (neg_b, b);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul6);
		addmul (R, M, a, A22, X2, neg_b, C21);
		SW_TIMER_STOP(mul6);

		SW_TIMER_START(other);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, A11, X1);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), A21, X1);

		BLAS3::_copy<Ring, ParentTag>::op (R, M, B22, X2);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.minusOne (), B12, X2);

		SW_TIMER_STOP(other);

		SW_TIMER_START(mul7);
		addmul (R, M, a, X1, X2, R.one (), X3);
		SW_TIMER_STOP(mul7);

		SW_TIMER_START(other);

		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X3, C22);
		BLAS3::_scal<Ring, ParentTag>::op (R, M, R.minusOne (), C21);
		BLAS3::_axpy<Ring, ParentTag>::op (R, M, R.one (), X3, C21);

		SW_TIMER_STOP(other);

		SW_TIMER_REPORT(mul1);
		SW_TIMER_REPORT(mul2);
		SW_TIMER_REPORT(mul3);
		SW_TIMER_REPORT(mul4);
		SW_TIMER_REPORT(mul5);
		SW_TIMER_REPORT(mul6);
		SW_TIMER_REPORT(mul7);
		SW_TIMER_REPORT(other);

#ifdef __LELA_SW_DETAILED_PROFILE
		commentator.stop (MSG_DONE);
#endif // __LELA_SW_DETAILED_PROFILE

		return gemm_res (R, M, a, A, B, b, C, m, k, n);
	}
}

// FIXME

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::mul_ip (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, Matrix3 &C)
{
	return C;
}

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::addmul_ip (const Ring &R, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
{
	return C;
}

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::mul_ow (const Ring &R, Modules &M, const typename Ring::Element &a, Matrix1 &A, Matrix2 &B, Matrix3 &C)
{
	return C;
}

template <class ParentTag>
template <class Ring, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &StrassenWinograd<ParentTag>::addmul_ow (const Ring &R, Modules &M, const typename Ring::Element &a, Matrix1 &A, Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
{
	return C;
}

} // namespace LELA

#endif // __LELA_ALGORITHMS_STRASSEN_WINOGRAD_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
