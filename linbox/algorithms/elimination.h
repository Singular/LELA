/* linbox/algorithms/elimination.h
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

#ifndef __LINBOX_ALGORITHMS_ELIMINATION_H
#define __LINBOX_ALGORITHMS_ELIMINATION_H

#include "linbox/vector/traits.h"
#include "linbox/vector/bit-iterator.h"
#include "linbox/ring/gf2.h"

namespace LinBox
{

/** Classical Gaussian elimiation
 *
 * This implements the naive Gaussian elimination algorithm. No
 * attempt is made to take advantage of fast matrix-multiplication,
 * cache-optimisation, or other things.
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

	// Find a suitable pivot from the matrix A starting at
	// start_row. If no pivot can be found (i.e. all rows
	// from start_row onwards are already 0) return
	// -1. Otherwise fill in col with the pivot-column.
	template <class Matrix>
	int GetPivot (const Matrix &A, int start_row, size_t &col) const
		{ return GetPivotSpecialised (A, start_row, col, typename VectorTraits<Ring, typename Matrix::Row>::RepresentationType ()); }

	// Find the first nonzero element in the given column
	// starting at the row of the same index. Return -1 if
	// none found.
	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorRepresentationTypes::Dense) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorRepresentationTypes::Dense01) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorRepresentationTypes::Sparse) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorRepresentationTypes::Sparse01) const;

	template <class Matrix>
	int GetPivotSpecialised (const Matrix &A, int start_row, size_t &col, VectorRepresentationTypes::Hybrid01) const;

public:
	/**
	 * \brief Constructor
	 *
	 * @param _ctx Context-object for computations
	 */
	Elimination (Context<Ring, Modules> &_ctx)
		: ctx (_ctx) {}

	/**
	 * \brief Compute the (reduced or non-reduced)
	 * row-echelon form of a matrix
	 *
	 * At conclusion, the parameters will have the property that
	 * A_out=LPA_in, where A_out is the matrix A at output and
	 * A_in is the matrix A at input. A_out is in reduced
	 * row-echelon form, L is lower triangular, and P is a
	 * permutation.
	 *
	 * In comparison with @see RowEchelonForm, this
	 * version does not take advantage of fast
	 * matrix-multiplication and does not use a
	 * divide-and-conquer method. It also modifies the
	 * input-matrix A.
	 *
	 * The pivot-strategy is to find the row with the
	 * fewest elements. This seems to be the only sensible
	 * approach, given that we are not allowed to permute
	 * columns.
	 *
	 * @param A The sparse matrix whose reduced
	 * row-echelon form to compute. Will be replaced by
	 * its reduced row-echelon form during computation.
	 *
	 * @param L The dense matrix into which to store the
	 * matrix L as defined above. Should be of size n x n,
	 * with n the row-dimension of A.
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
	 * @param reduced True if the routine should compute
	 * the reduced row-echelon form and false if it should
	 * only compute the (non-reduced) row-echelon form.
	 *
	 * @param compute_L True if the matrix L should be
	 * computed. If false, then L is ignored.
	 *
	 * @param start_row Start at this row. Intended for
	 * internal use.
	 */
	template <class Matrix1, class Matrix2>
	void RowEchelonForm (Matrix1       &A,
			     Matrix2       &L,
			     Permutation   &P,
			     size_t        &rank,
			     Element       &det,
			     bool           reduced = false,
			     bool           compute_L = true,
			     size_t         start_row = 0) const;

	/** \brief Take a matrix of known rank in row-echelon
	 * form and convert to reduced row-echelon form
	 *
	 * @param A Input matrix A in row-echelon form;
	 * replaced by its reduced row-echelon form
	 *
	 * @param L Dense matrix L into which to store
	 * conversion-matrix
	 *
	 * @param compute_L bool, true if L should be
	 * computed; if false then L is left unchanged
	 *
	 * @param rank Rank of A; must be known a priori
	 * (though can be easily computed by scanning rows of
	 * A)
	 *
	 * @param start_row Pivot-rows begin at this
	 * row. Intended for internal use only.
	 *
	 * @return Reference to A
	 */
	template <class Matrix1, class Matrix2>
	Matrix1 &ReduceRowEchelon (Matrix1 &A, Matrix2 &L, bool compute_L, size_t rank, size_t start_row = 0) const;

	/// Set the given matrix to the identity-matrix
	template <class Matrix>
	void SetIdentity (Matrix &U, size_t start_row = 0) const;
};

} // namespace LinBox

#include "linbox/algorithms/elimination.tcc"

#endif // __LINBOX_ALGORITHMS_ELIMINATION_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
