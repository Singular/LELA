/* lela/algorithms/elimination.h
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

#ifndef __LELA_ALGORITHMS_ELIMINATION_H
#define __LELA_ALGORITHMS_ELIMINATION_H

#include "lela/vector/traits.h"
#include "lela/vector/bit-iterator.h"
#include "lela/ring/gf2.h"
#include "lela/algorithms/pivot-strategy.h"

namespace LELA
{

/** Classical Gaussian elimiation
 *
 * This implements the naive Gaussian elimination algorithm. No
 * attempt is made to take advantage of fast matrix-multiplication,
 * cache-optimisation, or other things.
 *
 * \ingroup algorithms
 */

template <class Ring, class Modules = AllModules<Ring> >
class Elimination
{
public:
	typedef typename Ring::Element Element;
	typedef std::pair<uint32, uint32> Transposition;
	typedef std::vector<Transposition> Permutation;

private:
	Context<Ring, Modules> &ctx;

public:
	/**
	 * \brief Constructor
	 *
	 * @param _ctx Context-object for computations
	 */
	Elimination (Context<Ring, Modules> &_ctx)
		: ctx (_ctx) {}

	/**
	 * \brief Compute the (non-reduced)
	 * row-echelon form of a matrix
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in (non-reduced)
	 * row-echelon form, L is lower triangular, and P is a
	 * permutation.
	 *
	 * If compute_L is set to true, then L is computed and stored
	 * in A in place of the part under the main diagonal. The
	 * diagonal-entries of L are of course omitted and may assumed
	 * to be one.
	 *
	 * @param A The matrix whose row-echelon form is to be
	 * computed. Will be replaced by its row-echelon form during
	 * computation.
	 *
	 * @param P The permutation into which to store the
	 * permutation P as defined above.
	 *
	 * @param rank An integer into which to store the
	 * computed rank of A.
	 *
	 * @param det A ring-element into which to store the
	 * computed determinant of the submatrix of A formed
	 * by taking pivot-rows and -columns.
	 *
	 * @param compute_L True if the matrix L should be
	 * computed. If false, then L is ignored.
	 */
	template <class Matrix>
	Matrix &echelonize (Matrix        &A,
			    Permutation   &P,
			    size_t        &rank,
			    Element       &det,
			    bool           compute_L = true) const
		{ return echelonize (A, P, rank, det, typename DefaultPivotStrategy<Ring, Modules, typename Matrix::Row>::Strategy (ctx), compute_L); }

	/** Compute the non-reduced row-echelon form of a matrix using
	 * the pivot-strategy provided
	 *
	 * Identical to echelonize above, but uses the given
	 * pivot-strategy.
	 */
	template <class Matrix, class PivotStrategy>
	Matrix &echelonize (Matrix        &A,
			    Permutation   &P,
			    size_t        &rank,
			    Element       &det,
			    PivotStrategy  PS,
			    bool           compute_L) const;

	/**
	 * \brief Compute the reduced row-echelon form of a matrix
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in reduced
	 * row-echelon form, L is the transform-matrix, and P is a
	 * permutation.
	 *
	 * @param A The matrix whose reduced row-echelon form is to be
	 * computed. Will be replaced by its reduced row-echelon form
	 * during computation.
	 *
	 * @param L The matrix in which to store the
	 * transform-matrix. Should be square of dimension equal to
	 * the row-dimension of A.
	 *
	 * @param P The permutation into which to store the
	 * permutation P as defined above.
	 *
	 * @param rank An integer into which to store the
	 * computed rank of A.
	 *
	 * @param det A ring-element into which to store the
	 * computed determinant of the submatrix of A formed
	 * by taking pivot-rows and -columns.
	 *
	 * @param compute_L True if the matrix L should be
	 * computed. If false, then L is ignored.
	 *
	 * @returns Reference to A
	 */
	template <class Matrix1, class Matrix2>
	Matrix1 &echelonize_reduced (Matrix1       &A,
				     Matrix2       &L,
				     Permutation   &P,
				     size_t        &rank,
				     Element       &det,
				     bool           compute_L = false) const
		{ return echelonize_reduced (A, L, P, rank, det, typename DefaultPivotStrategy<Ring, Modules, typename Matrix1::Row>::Strategy (ctx), compute_L); }

	/** Compute the reduced row-echelon form of a matrix using the
	 * pivot-strategy provided
	 *
	 * Identical to echelonize_reduced above, but uses the given
	 * pivot-strategy.
	 */
	template <class Matrix1, class Matrix2, class PivotStrategy>
	Matrix1 &echelonize_reduced (Matrix1       &A,
				     Matrix2       &L,
				     Permutation   &P,
				     size_t        &rank,
				     Element       &det,
				     PivotStrategy  PS,
				     bool           compute_L) const;

	/**
	 * \brief Compute the PLUQ-decomposition of a matrix
	 *
	 * At conclusion, the parameters will have the property that
	 * A=PLUQ, where L is unit lower triangular, U is upper
	 * triangular, and P and Q are permutations.
	 *
	 * The matrices L and U are stored in place in A, with L
	 * occupying the part below the main diagonal and U occupying
	 * the part above. The entries on the diagonal of L are
	 * omitted and are taken to be one.
	 *
	 * The matrix A must support row-iterators. This method can
	 * not be used on sparse matrices with rows in the hybrid 0-1
	 * format. Its use on dense 0-1 matrices is not currently
	 * recommended.
	 *
	 * @param A The sparse matrix whose reduced row-echelon form
	 * to compute. Will be replaced by the matrices L and U during
	 * computation.
	 *
	 * @param P The permutation into which to store the
	 * permutation P as defined above.
	 *
	 * @param rank An integer into which to store the
	 * computed rank of A.
	 *
	 * @param det A ring-element into which to store the
	 * computed determinant of the submatrix of A formed
	 * by taking pivot-rows and -columns.
	 */
	template <class Matrix>
	Matrix &pluq (Matrix        &A,
		      Permutation   &P,
		      Permutation   &Q,
		      size_t        &rank,
		      Element       &det) const
		{ return pluq (A, P, Q, rank, det, typename DefaultPivotStrategy<Ring, Modules, typename Matrix::Row>::Strategy (ctx)); }

	template <class Matrix, class PivotStrategy>
	Matrix &pluq (Matrix        &A,
		      Permutation   &P,
		      Permutation   &Q,
		      size_t        &rank,
		      Element       &det,
		      PivotStrategy  PS) const;

	/** Move the lower-triangular part of A to L and reset the
	 * lower-triangular part of A to 0
	 *
	 * This is a convenience-function. It is not written for
	 * efficiency and may perform quite poorly on 0-1 matrices.
	 *
	 * Assumes L is preset to 0
	 */
	template <class Matrix1, class Matrix2>
	void move_L (Matrix1 &L, Matrix2 &A) const;
};

} // namespace LELA

#include "lela/algorithms/elimination.tcc"

#endif // __LELA_ALGORITHMS_ELIMINATION_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
