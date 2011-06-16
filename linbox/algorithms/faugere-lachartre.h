/* linbox/algorithms/faugere-lachartre.h
 * Copyright 2010 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 *
 * Implementation of algorithm for computing reduced row-echelon form
 * of a matrix coming from the F4-algorithm
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ALGORITHMS_FAUGERE_LACHARTRE_H
#define __LINBOX_ALGORITHMS_FAUGERE_LACHARTRE_H

#include "linbox/blas/context.h"
#include "linbox/algorithms/gauss-jordan.h"
#include "linbox/solutions/echelon-form.h"
#include "linbox/solutions/echelon-form-gf2.h"
#include "linbox/util/splicer.h"

namespace LinBox
{

/**
 * \brief Implementation of algorithm for computing reduced row-echelon form
 * of a matrix coming from the F4-algorithm
 *
 * The algorithm used is based loosely on the algorithm
 * described in J.C. Faugère's paper "Parallel Gaussian
 * Elimination for Gröbner bases computations in finite
 * fields", PASCO 2010.
 */
template <class Ring, class Modules = AllModules<Ring> >
class FaugereLachartre {
public:
	typedef typename Ring::Element Element;
	typedef typename GaussJordan<Ring>::SparseMatrix SparseMatrix;
	typedef typename GaussJordan<Ring>::DenseMatrix DenseMatrix;

private:
	Context<Ring, Modules> &ctx;
	EchelonForm<Ring, Modules> EF;

	template <class Matrix>
	void setup_splicer (Splicer &splicer, Splicer &reconst_splicer, const Matrix &A, size_t &num_pivot_rows, typename Ring::Element &det) const;

public:
	/**
	 * \brief Construct a new FaugereLachartre
	 *
	 * @param _ctx Context-object for matrix-calculations
	 */
	FaugereLachartre (Context<Ring, Modules> &_ctx)
		: ctx (_ctx), EF (_ctx) {}

	/** 
	 * \brief Convert the matrix A into reduced
	 * row-echelon form
	 *
	 * @param R Matrix into which to store the reduced
	 * row-echelon form. Should have the same dimensions
	 * as A. May be the same matrix as A, in which case A
	 * is replaced by its reduced row-echelon form.
	 *
	 * @param A Matrix to be converted to row-echelon
	 * form. Not altered, unless R is the same matrix.
	 *
	 * @param rank Integer-reference into which to store
	 * computed rank
	 *
	 * @param det Ring-element-reference into which to
	 * store computed determinant of pivot-submatrix
	 */
	void RowEchelonForm (SparseMatrix &R, const SparseMatrix &X, size_t &rank, Element &det);
};

} // namespace LinBox

#include "linbox/algorithms/faugere-lachartre.tcc"

#endif // __LINBOX_ALGORITHMS_FAUGERE_LACHARTRE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
