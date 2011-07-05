/* linbox/algorithms/elimination.tcc
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Naive Gaussian elimination
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ALGORITHMS_ELIMINATION_TCC
#define __LINBOX_ALGORITHMS_ELIMINATION_TCC

#include <iostream>
#include <iomanip>
#include <cassert>

#include "linbox/algorithms/elimination.h"
#include "linbox/blas/level1.h"
#include "linbox/blas/level3.h"

#ifdef DETAILED_PROFILE
#  define TIMER_DECLARE(part) LinBox::UserTimer part##_timer; double part##_time = 0.0;
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

namespace LinBox
{

template <class Ring, class Modules>
template <class Matrix>
void Elimination<Ring, Modules>::SetIdentity (Matrix &U, size_t start_row) const
{
	StandardBasisStream<Ring, typename Matrix::Row> stream (ctx.F, U.coldim ());
	typename Matrix::RowIterator i = U.rowBegin () + start_row;

	for (; i != U.rowEnd (); ++i)
		stream >> *i;
}

template <class Ring, class Modules>
template <class Matrix>
int Elimination<Ring, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						     VectorRepresentationTypes::Dense) const
{
	typename Matrix::ConstRowIterator i;

	for (; col < A.coldim (); ++col) {
		int k = start_row;

		for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k)
			if (!ctx.F.isZero ((*i)[col]))
				return k;
	}

	return -1;
}

template <class Ring, class Modules>
template <class Matrix>
int Elimination<Ring, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						     VectorRepresentationTypes::Dense01) const
{
	typename Matrix::ConstRowIterator i;

	for (; col < A.coldim (); ++col) {
		int k = start_row;

		for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k)
			if ((*i)[col])
				return k;
	}

	return -1;
}

template <class Ring, class Modules>
template <class Matrix>
int Elimination<Ring, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						     VectorRepresentationTypes::Sparse) const
{
	typename Matrix::ConstRowIterator i;

	size_t min_nonzero = (size_t) -1, pivot = -1, k = start_row;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty ()) {
			if (i->front ().first < col) {
				col = i->front ().first;
				min_nonzero = i->size ();
				pivot = k;
			}
			else if (i->front ().first == col && i->size () < min_nonzero) {
				min_nonzero = i->size ();
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Ring, class Modules>
template <class Matrix>
int Elimination<Ring, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						     VectorRepresentationTypes::Sparse01) const
{
	typename Matrix::ConstRowIterator i;

	size_t min_nonzero = (size_t) -1, pivot = -1, k = start_row, s;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty ()) {
			if (i->front () < col) {
				col = i->front ();
				min_nonzero = i->size ();
				pivot = k;
			}
			else if (i->front () == col && (s = i->size ()) < min_nonzero) {
				min_nonzero = s;
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Ring, class Modules>
template <class Matrix>
int Elimination<Ring, Modules>::GetPivotSpecialised (const Matrix &A, int start_row, size_t &col,
						     VectorRepresentationTypes::Hybrid01) const
{
	typename Matrix::ConstRowIterator i;
	typename Ring::Element a;

	size_t min_blocks = (size_t) -1, pivot = -1, k = start_row, s;
	col = A.coldim ();

	for (i = A.rowBegin () + start_row; i != A.rowEnd (); ++i, ++k) {
		if (!i->empty () && i->front ().first << WordTraits<typename Matrix::Row::word_type>::logof_size <= (int) col) {
			int idx = BLAS1::head (ctx, a, *i);
			if ((size_t) idx < col) {
				col = idx;
				min_blocks = i->size ();
				pivot = k;
			}
			else if ((size_t) idx == col && (s = i->size ()) < min_blocks) {
				min_blocks = s;
				pivot = k;
			}
		}
	}

	return pivot;
}

template <class Ring, class Modules>
template <class Matrix1, class Matrix2>
void Elimination<Ring, Modules>::RowEchelonForm (Matrix1       &A,
						 Matrix2       &L,
						 Permutation   &P,
						 size_t        &rank,
						 Element       &det,
						 bool           reduced,
						 bool           compute_L,
						 size_t         start_row) const
{
	commentator.start ("Standard row-echelon form", __FUNCTION__, A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	TIMER_DECLARE(GetPivot);
	TIMER_DECLARE(Permute);
	TIMER_DECLARE(ElimBelow);

	typename Matrix1::RowIterator i_A, j_A;

	size_t col = 0;
	int k = start_row;
	Element a, x, xinv, negxinv, negaxinv;

	ctx.F.assign (a, ctx.F.zero ());
	ctx.F.assign (x, ctx.F.zero ());

	typename Matrix2::RowIterator i_L, j_L;

	if (compute_L) {
		SetIdentity (L, start_row);
		i_L = L.rowBegin () + k;
	}

	P.clear ();
	rank = 0;
	ctx.F.init (det, 1);

	for (i_A = A.rowBegin () + k; i_A != A.rowEnd (); ++k, ++i_A) {
		TIMER_START(GetPivot);
		int pivot = GetPivot (A, k, col);
		TIMER_STOP(GetPivot);

		if (pivot == -1)
			break;

		TIMER_START(Permute);
		if (k != pivot) {
			// DEBUG
			// report << "Permuting " << k << " and " << pivot << " (first " << k - start_row << " columns)" << std::endl;
			// report << "L before permutation:" << std::endl;
			// BLAS3::write (ctx, report, L);

			Transposition t (k, pivot);
			P.push_back (t);
			BLAS3::permute_rows (ctx, &t, &t + 1, A);

			if (compute_L) {
				typename Matrix2::SubmatrixType Lp (L, 0, 0, L.rowdim (), k - start_row);
				BLAS3::permute_rows (ctx, &t, &t + 1, Lp);
			}

			// DEBUG
			// report << "L after permutation:" << std::endl;
			// BLAS3::write (ctx, report, L);
		}
		TIMER_STOP(Permute);

		// DEBUG
		// report << "Row " << k << ", pivot is " << pivot << std::endl;
		// report << "Current A:" << std::endl;
		// BLAS3::write (ctx, report, A);
		// report << "Current L:" << std::endl;
		// BLAS3::write (ctx, report, L);

		BLAS1::head (ctx, x, *i_A);

		ctx.F.mulin (det, x);

		ctx.F.inv (xinv, x);
		ctx.F.neg (negxinv, xinv);

		TIMER_START(ElimBelow);

		if (compute_L) {
			j_L = i_L;
			++j_L;
		}

		for (j_A = i_A; ++j_A != A.rowEnd ();) {
			if (BLAS1::head (ctx, a, *j_A) == (int) col) {
				// DEBUG
				// report << "Eliminating row " << j_A - A.rowBegin () << " from row " << k << std::endl;

				ctx.F.mul (negaxinv, negxinv, a);
				BLAS1::axpy (ctx, negaxinv, *i_A, *j_A);

				if (compute_L)
					BLAS1::axpy (ctx, negaxinv, *i_L, *j_L);
			}

			if (compute_L)
				++j_L;
		}
		TIMER_STOP(ElimBelow);

		++rank;

		if ((i_A - A.rowBegin ()) % PROGRESS_STEP == PROGRESS_STEP - 1)
			commentator.progress ();

		if (compute_L)
			++i_L;
	}

	if (reduced)
		ReduceRowEchelon (A, L, compute_L, rank + start_row, start_row);

	TIMER_REPORT(GetPivot);
	TIMER_REPORT(Permute);
	TIMER_REPORT(ElimBelow);

	commentator.stop (MSG_DONE);
}

template <class Ring, class Modules>
template <class Matrix1, class Matrix2>
Matrix1 &Elimination<Ring, Modules>::ReduceRowEchelon (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row) const
{
	commentator.start ("Reducing row-echelon form", "Elimination::ReduceRowEchelon", A.rowdim () / PROGRESS_STEP);

	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	if (rank == 0)
		return A;

	typename Matrix1::RowIterator i_A, j_A;
	typename Matrix2::RowIterator i_L;

	if (compute_L)
		i_L = L.rowBegin () + (rank - 1);

	int current_row = rank - 1, elim_row;

	typename Ring::Element a, x;

	i_A = A.rowBegin () + current_row;

	do {
		// DEBUG
		// report << "Row " << current_row << ", current A:" << std::endl;
		// BLAS3::write (ctx, report, A);

		for (elim_row = rank - 1, j_A = A.rowBegin () + elim_row; elim_row > std::max (current_row, (int) start_row - 1); --elim_row, --j_A) {
			size_t col = BLAS1::head (ctx, a, *j_A);

			if (VectorUtils::getEntry (*i_A, x, col) && !ctx.F.isZero (x)) {
				// DEBUG
				// report << "Eliminating " << current_row << " from " << elim_row << std::endl;

				BLAS1::axpy (ctx, ctx.F.one (), *j_A, *i_A);

				if (compute_L)
					BLAS1::axpy (ctx, ctx.F.one (), *(L.rowBegin () + elim_row), *i_L);
			}
		}

		if (compute_L)
			--i_L;

		if ((rank - current_row) % PROGRESS_STEP == 0)
			commentator.progress ();

		--current_row;
	} while (i_A-- != A.rowBegin ());

	commentator.stop (MSG_DONE, NULL, "Elimination::ReduceRowEchelon");

	return A;
}

} // namespace LinBox

#endif // __LINBOX_ALGORITHMS_ELIMINATION_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
