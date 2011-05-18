/* linbox/solutions/echelon-form.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Compute the echelon-form of a matrix
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_SOLUTIONS_ECHELON_FORM_H
#define __LINBOX_SOLUTIONS_ECHELON_FORM_H

#include "linbox/blas/context.h"
#include "linbox/algorithms/gauss-jordan.h"
#include "linbox/matrix/dense.h"
#include "linbox/util/error.h"

namespace LinBox
{

/** Solution for computing the echelon-form of a matrix */
template <class Field, class Modules = AllModules<Field> >
class EchelonForm
{
	GaussJordan<Field, Modules> _GJ;

	DenseMatrix<typename Field::Element> _L;
	typename MatrixDomain<Field>::Permutation _P;

	// Map pointers to matrices to computed ranks
	std::map<const void *, size_t> _rank_table;

public:
	enum Method { METHOD_STANDARD_GJ, METHOD_ASYMPTOTICALLY_FAST_GJ };

	/** Constructor
	 *
	 * @param F Field over which to compute
	 */
	EchelonForm (Context<Field, Modules> &ctx) : _GJ (ctx) {}

	/** Compute the (possibly reduced) row-echelon form of a matrix
	 *
	 * @param A Input matrix, to be replaced by its row-echelon form
	 * @param reduced true if reduced row-echelon form should be computed, false if not
	 * @param method Method to be used. Must be METHOD_STANDARD_GJ if the matrix is not dense.
	 * @returns Reference to A
	 */
	template <class Matrix>
	Matrix &RowEchelonForm (Matrix &A, bool reduced = false, Method method = METHOD_STANDARD_GJ)
	{
		commentator.start ("Row-echelon form", __FUNCTION__);

		size_t rank;
		typename Field::Element d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_GJ.StandardRowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		_rank_table[&A] = rank;

		commentator.stop (MSG_DONE);

		return A;
	}

	// Specialisation for dense matrices
	DenseMatrix<typename Field::Element> &RowEchelonForm (DenseMatrix<typename Field::Element> &A, bool reduced = false, Method method = METHOD_ASYMPTOTICALLY_FAST_GJ)
	{
		commentator.start ("Row-echelon form", __FUNCTION__);

		size_t rank;
		typename Field::Element d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_GJ.StandardRowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		case METHOD_ASYMPTOTICALLY_FAST_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_GJ.DenseRowEchelonForm (A, _L, _P, rank, d, reduced);
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		_rank_table[&A] = rank;

		commentator.stop (MSG_DONE);

		return A;
	}

	/** Determine the rank of the given matrix
	 *
	 * @param A Input matrix. Must already have been an argument to RowEchelonForm.
	 * @returns rank
	 */
	template <class Matrix>
	size_t rank (const Matrix &A) const
		{ return _rank_table[&A]; }
};

} // namespace LinBox

#endif // __LINBOX_SOLUTIONS_ECHELON_FORM_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
