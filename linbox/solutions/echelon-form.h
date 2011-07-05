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
#include "linbox/algorithms/elimination.h"
#include "linbox/algorithms/gauss-jordan.h"
#include "linbox/algorithms/faugere-lachartre.h"
#include "linbox/matrix/dense.h"
#include "linbox/util/error.h"

namespace LinBox
{

/** Solution for computing the echelon-form of a matrix */
template <class Ring, class Modules = AllModules<Ring> >
class EchelonForm
{
	Context<Ring, Modules> _ctx;
	Elimination<Ring, Modules> _elim;
	GaussJordan<Ring, Modules> _GJ;

	DenseMatrix<typename Ring::Element> _L;
	typename GaussJordan<Ring, Modules>::Permutation _P;

	// Map pointers to matrices to computed ranks
	std::map<const void *, size_t> _rank_table;

public:
	enum Method { METHOD_UNKNOWN, METHOD_STANDARD_GJ, METHOD_ASYMPTOTICALLY_FAST_GJ, METHOD_FAUGERE_LACHARTRE };

	/** Constructor
	 *
	 * @param F Ring over which to compute
	 */
	EchelonForm (Context<Ring, Modules> &ctx) : _ctx (ctx), _elim (ctx), _GJ (ctx) {}

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
		static const char *method_names[] = { "unknown", "standard", "recursive", "M4RI", "Faugère-Lachartre" };

		std::ostringstream str;
		str << "Row-echelon form (method: " << method_names[method] << ")" << std::ends;

		commentator.start (str.str ().c_str (), __FUNCTION__);

		size_t rank;
		typename Ring::Element d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_elim.RowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		case METHOD_FAUGERE_LACHARTRE:
			if (reduced) {
				// Must do it this way to avoid an infinite loop of inclusion...
				FaugereLachartre<Ring, Modules> FL (_ctx);
				FL.RowEchelonForm (A, A, rank, d);
			} else
				throw LinboxError ("Only reduced row-echelon form is available with Faugère-Lachartre");
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		_rank_table[&A] = rank;

		commentator.stop (MSG_DONE);

		return A;
	}

	// Specialisation for dense matrices
	DenseMatrix<typename Ring::Element> &RowEchelonForm (DenseMatrix<typename Ring::Element> &A, bool reduced = false, Method method = METHOD_ASYMPTOTICALLY_FAST_GJ)
	{
		static const char *method_names[] = { "unknown", "standard", "recursive", "M4RI", "Faugère-Lachartre" };

		std::ostringstream str;
		str << "Row-echelon form (method: " << method_names[method] << ")" << std::ends;

		commentator.start (str.str ().c_str (), __FUNCTION__);

		size_t rank;
		typename Ring::Element d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_elim.RowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		case METHOD_ASYMPTOTICALLY_FAST_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_GJ.RowEchelonForm (A, _L, _P, rank, d, reduced);
			break;

		case METHOD_FAUGERE_LACHARTRE:
			if (reduced) {
				// Must do it this way to avoid an infinite loop of inclusion...
				FaugereLachartre<Ring, Modules> FL (_ctx);
				FL.RowEchelonForm (A, A, rank, d);
			} else
				throw LinboxError ("Only reduced row-echelon form is available with Faugère-Lachartre");
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
