/* lela/algorithms/gauss-jordan.h
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Asymptotically fast Gauss-Jordan elimination
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ALGORITHMS_GAUSS_JORDAN_H
#define __LELA_ALGORITHMS_GAUSS_JORDAN_H

#include "lela/util/commentator.h"
#include "lela/util/timer.h"
#include "lela/blas/context.h"
#include "lela/ring/gf2.h"
#include "lela/matrix/dense.h"
#include "lela/matrix/sparse.h"
#include "lela/algorithms/pivot-strategy.h"

namespace LELA
{

/** Implementation of asymptotically fast Gauss-Jordan elimination
 *
 * \ingroup algorithms
 */
template <class Ring, class Modules = AllModules<Ring> >
class GaussJordan
{
public:
	typedef typename Ring::Element Element;
	typedef std::pair<uint32, uint32> Transposition;
	typedef std::vector<Transposition> Permutation;

private:
	Context<Ring, Modules> &ctx;

	struct CompareSecond
	{
		bool operator () (const Transposition &t1, const Transposition &t2) const { return t1.second < t2.second; }
	};

	// Internal recursive procedure for the kth indexed
	// Gauss-Jordan transform. Uses a divide and conquer
	// method to maximise use of fast
	// matrix-multiplication. See Chapter 2 of "Algorithms
	// for Matrix Canonical Forms", Ph.D thesis by Arne
	// Storjohann.
	template <class Matrix1, class Matrix2, class PivotStrategy>
	void GaussJordanTransform (Matrix1      &A,
				   int           k,
				   Element       d_0,
				   Matrix2      &U,
				   Permutation  &P,
				   size_t       &r,
				   int          &h,
				   Element      &d,
				   DenseMatrix<typename Ring::Element> &S,
				   DenseMatrix<typename Ring::Element> &T,
				   PivotStrategy PS) const;

	// Internal recursive procedure for the Gauss transform. Uses
	// a divide and conquer method to maximise use of fast
	// matrix-multiplication. See Chapter 2 of "Algorithms for
	// Matrix Canonical Forms", Ph.D thesis by Arne Storjohann.
	template <class Matrix1, class Matrix2, class PivotStrategy>
	void GaussTransform (Matrix1     &A,
			     Matrix2     &L,
			     Element      d_0,
			     Permutation &P,
			     size_t      &r,
			     int         &h,
			     Element     &d,
			     PivotStrategy PS) const;

public:
	/**
	 * \brief Constructor
	 *
	 * @param _F Ring over which operations take place
	 */
	GaussJordan (Context<Ring, Modules> &_ctx)
		: ctx (_ctx)
		{}

	/**
	 * \brief Convert the matrix A into (non-reduced) row-echelon
	 * form.
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in (non-reduced)
	 * row-echelon form, L is lower triangular, and P is a
	 * permutation.
	 *
	 * L is computed and placed in the lower triangular part of
	 * A. Its diagonal-entries are omitted and taken in any case
	 * to be one.
	 *
	 * @param A Matrix A. Will be replaced by its reduced
	 * row-echelon form and L
	 *
	 * @param P Permutation in which to store the
	 * row-permutations of A made by the choice of pivots.
	 *
	 * @param rank Integer into which to store the
	 * computed rank
	 *
	 * @param det Ring-element into which to store the
	 * computed determinant
	 *
	 * @returns Reference to A
	 */
	template <class Matrix>
	Matrix &echelonize (Matrix      &A,
			    Permutation &P,
			    size_t      &rank,
			    Element     &det)
		{ return echelonize (A, P, rank, det, typename DefaultPivotStrategy<Ring, Modules, typename Matrix::Row>::Strategy (ctx)); }

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
			    PivotStrategy  PS);

	/**
	 * \brief Convert the matrix A into reduced row-echelon
	 * form.
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in reduced
	 * row-echelon form and P is a permutation.
	 *
	 * @param A Matrix A. Will be replaced by its reduced
	 * row-echelon form
	 *
	 * @param L Matrix in which to store L
	 *
	 * @param P Permutation in which to store the
	 * row-permutations of A made by the choice of pivots.
	 *
	 * @param rank Integer into which to store the
	 * computed rank
	 *
	 * @param det Ring-element into which to store the
	 * computed determinant
	 */
	template <class Matrix1, class Matrix2>
	Matrix1 &echelonize_reduced (Matrix1     &A,
				     Matrix2     &L,
				     Permutation &P,
				     size_t      &rank,
				     Element     &det)
		{ return echelonize_reduced (A, L, P, rank, det, typename DefaultPivotStrategy<Ring, Modules, typename Matrix1::Row>::Strategy (ctx)); }

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
				     PivotStrategy  PS);
};

} // namespace LELA

#include "lela/algorithms/gauss-jordan.tcc"

#endif // __LELA_ALGORITHMS_GAUSS_JORDAN_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
