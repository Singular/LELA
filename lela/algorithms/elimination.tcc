/* lela/algorithms/elimination.tcc
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Naive Gaussian elimination
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ALGORITHMS_ELIMINATION_TCC
#define __LELA_ALGORITHMS_ELIMINATION_TCC

#include <iostream>
#include <iomanip>
#include <cassert>

#include "lela/algorithms/elimination.h"
#include "lela/blas/level1.h"
#include "lela/blas/level3.h"
#include "lela/vector/stream.h"

#ifdef DETAILED_PROFILE
#  define TIMER_DECLARE(part) LELA::UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_START(part) part##_timer.start ()
#  define TIMER_STOP(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define TIMER_REPORT(part) \
	commentator.report (Commentator::LEVEL_NORMAL, TIMING_MEASURE) \
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
#else
#  define TIMER_DECLARE(part)
#  define TIMER_START(part)
#  define TIMER_STOP(part)
#  define TIMER_REPORT(part)
#endif

#ifndef PROGRESS_STEP
#  define PROGRESS_STEP 1024
#endif // PROGRESS_STEP

namespace LELA
{

template <class Ring, class Modules>
template <class Matrix, class PivotStrategy>
Matrix &Elimination<Ring, Modules>::echelonize (Matrix        &A,
						Permutation   &P,
						size_t        &rank,
						Element       &det,
						PivotStrategy  PS,
						bool           compute_L) const
{
	commentator.start ("Echelonize (elimination)", __FUNCTION__, A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	TIMER_DECLARE(GetPivot);
	TIMER_DECLARE(Permute);
	TIMER_DECLARE(ElimBelow);

	typename Matrix::RowIterator i_A, j_A = A.rowBegin ();

	size_t pivot_row, pivot_col;
	size_t i, j;
	Element a, x, negxinv, negaxinv;

	ctx.F.copy (a, ctx.F.zero ());
	ctx.F.copy (x, ctx.F.zero ());

	P.clear ();
	rank = 0;
	ctx.F.init (det, 1);

	for (i_A = A.rowBegin (), i = 0, pivot_col = 0; i_A != A.rowEnd () && pivot_col < A.coldim (); ++i, ++i_A, ++pivot_col) {
		TIMER_START(GetPivot);
		pivot_row = i;
		if (!PS.getPivot (A, x, pivot_row, pivot_col))
			break;
		TIMER_STOP(GetPivot);

		lela_check (pivot_row < A.rowdim ());
		lela_check (pivot_col < A.coldim ());

		// DEBUG
		// report << "Row " << i << ", pivot is (" << pivot_row << "," << pivot_col << ")" << std::endl;
		// report << "Current A:" << std::endl;
		// BLAS3::write (ctx, report, A);

		TIMER_START(Permute);
		if (i != pivot_row) {
			Transposition t (i, pivot_row);
			P.push_back (t);
			BLAS3::permute_rows (ctx, &t, &t + 1, A);
		}
		TIMER_STOP(Permute);

		// report << "A after permutation:" << std::endl;
		// BLAS3::write (ctx, report, A);

		ctx.F.mulin (det, x);

		if (!ctx.F.inv (negxinv, x))
			throw LELAError ("Could not invert pivot-element in the ring");

		ctx.F.negin (negxinv);

		TIMER_START(ElimBelow);

		for (j_A = i_A, j = i + 1; ++j_A != A.rowEnd (); ++j) {
			if (A.getEntry (a, j, pivot_col) && !ctx.F.isZero (a)) {
				// DEBUG
				// report << "Eliminating row " << j << " from row " << i << std::endl;

				ctx.F.mul (negaxinv, a, negxinv);
				BLAS1::axpy (ctx, negaxinv, *i_A, *j_A);

				if (compute_L)
					A.setEntry (j, i, negaxinv);
			}
		}
		TIMER_STOP(ElimBelow);

		++rank;

		if (i % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	}

	TIMER_REPORT(GetPivot);
	TIMER_REPORT(Permute);
	TIMER_REPORT(ElimBelow);

	commentator.stop (MSG_DONE);

	return A;
}

template <class Ring, class Modules>
template <class Matrix1, class Matrix2, class PivotStrategy>
Matrix1 &Elimination<Ring, Modules>::echelonize_reduced (Matrix1       &A,
							 Matrix2       &L,
							 Permutation   &P,
							 size_t        &rank,
							 Element       &det,
							 PivotStrategy  PS,
							 bool           compute_L) const
{
	lela_check (!compute_L || L.rowdim () == A.rowdim ());
	lela_check (!compute_L || L.coldim () == A.rowdim ());

	commentator.start ("Echelonize (elimination)", __FUNCTION__, A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	TIMER_DECLARE(GetPivot);
	TIMER_DECLARE(Permute);
	TIMER_DECLARE(Elim);

	typename Matrix1::RowIterator i_A, j_A;
	typename Matrix2::RowIterator i_L = L.rowBegin (), j_L = L.rowBegin ();

	size_t pivot_row, pivot_col;
	size_t i, j;
	Element a, x, xinv, nega;

	ctx.F.copy (a, ctx.F.zero ());
	ctx.F.copy (x, ctx.F.zero ());

	P.clear ();
	rank = 0;
	ctx.F.init (det, 1);

	if (compute_L) {
		// Set L initially to the identity-matrix
		StandardBasisStream<Ring, typename Matrix2::Row> s (ctx.F, L.rowdim ());

		for (i_L = L.rowBegin (); i_L != L.rowEnd (); ++i_L)
			s >> *i_L;

		i_L = L.rowBegin ();
	}

	for (i_A = A.rowBegin (), i = 0, pivot_col = 0; i_A != A.rowEnd () && pivot_col < A.coldim (); ++i, ++i_A, ++pivot_col) {
		TIMER_START(GetPivot);
		pivot_row = i;
		if (!PS.getPivot (A, x, pivot_row, pivot_col))
			break;
		TIMER_STOP(GetPivot);

		lela_check (pivot_row < A.rowdim ());
		lela_check (pivot_col < A.coldim ());

		// DEBUG
		// report << "Row " << i << ", pivot is (" << pivot_row << "," << pivot_col << ")" << std::endl;
		// report << "Current A:" << std::endl;
		// BLAS3::write (ctx, report, A);

		TIMER_START(Permute);
		if (i != pivot_row) {
			Transposition t (i, pivot_row);
			P.push_back (t);
			BLAS3::permute_rows (ctx, &t, &t + 1, A);

			if (compute_L) {
				BLAS3::permute_rows (ctx, &t, &t + 1, L);
				BLAS3::permute_cols (ctx, &t, &t + 1, L);
			}
		}
		TIMER_STOP(Permute);

		// report << "A after permutation:" << std::endl;
		// BLAS3::write (ctx, report, A);

		ctx.F.mulin (det, x);

		if (!ctx.F.inv (xinv, x))
			throw LELAError ("Could not invert pivot-element in the ring");

		TIMER_START(Elim);
		BLAS1::scal (ctx, xinv, *i_A);

		if (compute_L) {
			BLAS1::scal (ctx, xinv, *i_L);
			j_L = L.rowBegin ();
		}

		for (j_A = A.rowBegin (), j = 0; j_A != A.rowEnd (); ++j, ++j_A) {
			if (j_A != i_A && A.getEntry (a, j, pivot_col) && !ctx.F.isZero (a)) {
				// DEBUG
				// report << "Eliminating row " << j << " from row " << i << std::endl;

				ctx.F.neg (nega, a);
				BLAS1::axpy (ctx, nega, *i_A, *j_A);
				
				if (compute_L)
					BLAS1::axpy (ctx, nega, *i_L, *j_L);
			}

			if (compute_L)
				++j_L;
		}

		TIMER_STOP(Elim);

		if (compute_L)
			++i_L;

		++rank;

		if (i % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	}

	TIMER_REPORT(GetPivot);
	TIMER_REPORT(Permute);
	TIMER_REPORT(Elim);

	commentator.stop (MSG_DONE);

	return A;
}

template <class Ring, class Modules>
template <class Matrix, class PivotStrategy>
Matrix &Elimination<Ring, Modules>::pluq (Matrix        &A,
					  Permutation   &P,
					  Permutation   &Q,
					  size_t        &rank,
					  Element       &det,
					  PivotStrategy  PS) const
{
	commentator.start ("PLUQ-decomposition (elimination)", __FUNCTION__, A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	TIMER_DECLARE(GetPivot);
	TIMER_DECLARE(Permute);
	TIMER_DECLARE(ElimBelow);

	typename Matrix::RowIterator i_A, j_A;

	size_t pivot_row, pivot_col;
	size_t i, j;
	Element a, x, negxinv, negaxinv;

	ctx.F.copy (a, ctx.F.zero ());
	ctx.F.copy (x, ctx.F.zero ());

	P.clear ();
	Q.clear ();
	rank = 0;
	ctx.F.init (det, 1);

	for (i_A = A.rowBegin (), i = 0; i_A != A.rowEnd () && i < A.coldim (); ++i, ++i_A) {
		TIMER_START(GetPivot);
		pivot_row = pivot_col = i;
		if (!PS.getPivot (A, x, pivot_row, pivot_col))
			break;
		TIMER_STOP(GetPivot);

		lela_check (pivot_row < A.rowdim ());
		lela_check (pivot_col < A.coldim ());

		TIMER_START(Permute);
		if (i != pivot_row) {
			Transposition t (i, pivot_row);
			P.push_back (t);
			BLAS3::permute_rows (ctx, &t, &t + 1, A);
		}

		if (i != pivot_col) {
			Transposition t (i, pivot_col);
			Q.push_back (t);
			BLAS3::permute_cols (ctx, &t, &t + 1, A);
		}
		TIMER_STOP(Permute);

		// DEBUG
		// report << "Row " << i << ", pivot-row is " << pivot << ", pivot-col is " << col << std::endl;
		// report << "Current A:" << std::endl;
		// BLAS3::write (ctx, report, A);

		ctx.F.mulin (det, x);

		if (!ctx.F.inv (negxinv, x))
			throw LELAError ("Could not invert pivot-element in the ring");

		ctx.F.negin (negxinv);

		TIMER_START(ElimBelow);

		for (j_A = i_A, j = i + 1; ++j_A != A.rowEnd (); ++j) {
			if (A.getEntry (a, j, i) && !ctx.F.isZero (a)) {
				// DEBUG
				// report << "Eliminating row " << j << " from row " << i << std::endl;

				// Use subvectors here to avoid trampling on content of L
				typename VectorTraits<Ring, typename Matrix::Row>::SubvectorType i_sub (*i_A, i + 1, A.coldim ());
				typename VectorTraits<Ring, typename Matrix::Row>::SubvectorType j_sub (*j_A, i + 1, A.coldim ());

				ctx.F.mul (negaxinv, a, negxinv);
				BLAS1::axpy (ctx, negaxinv, i_sub, j_sub);

				// Set entry in L
				ctx.F.negin (negaxinv);
				A.setEntry (j, i, negaxinv);
			}
		}
		TIMER_STOP(ElimBelow);

		++rank;

		if (i % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();
	}

	std::reverse (P.begin (), P.end ());
	std::reverse (Q.begin (), Q.end ());

	TIMER_REPORT(GetPivot);
	TIMER_REPORT(Permute);
	TIMER_REPORT(ElimBelow);

	commentator.stop (MSG_DONE);

	return A;
}

template <class Ring, class Modules>
template <class Matrix1, class Matrix2>
void Elimination<Ring, Modules>::move_L (Matrix1 &L, Matrix2 &A) const
{
	lela_check (L.rowdim () == A.rowdim ());
	lela_check (std::min (A.rowdim (), A.coldim ()) <= L.coldim ());

	size_t i, j;
	typename Ring::Element aij;

	for (i = 0; i < L.rowdim (); ++i) {
		for (j = 0; j < std::min (i, A.coldim ()); ++j) {
			if (A.getEntry (aij, i, j)) {
				L.setEntry (i, j, aij);
				A.setEntry (i, j, ctx.F.zero ());
				A.eraseEntry (i, j);
			}
		}
	}
}

} // namespace LELA

#endif // __LELA_ALGORITHMS_ELIMINATION_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
