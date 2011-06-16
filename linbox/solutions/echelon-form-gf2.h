/* linbox/solutions/echelon-form-gf2.h
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

#ifndef __LINBOX_SOLUTIONS_ECHELON_FORM_GF2_H
#define __LINBOX_SOLUTIONS_ECHELON_FORM_GF2_H

// No point specialising if we aren't using libm4ri
#ifdef __LINBOX_HAVE_M4RI

#include <m4ri/m4ri.h>

#include "linbox/util/commentator.h"
#include "linbox/ring/gf2.h"
#include "linbox/solutions/echelon-form.h"
#include "linbox/matrix/m4ri-matrix.h"
#include "linbox/blas/level1.h"
#include "linbox/blas/level3.h"

namespace LinBox
{

// Specialisation of EchelonForm to GF2 to take advantage of M4RI-routines
template <>
class EchelonForm<GF2, AllModules<GF2> >
{
	Context<GF2, AllModules<GF2> > &_ctx;
	GaussJordan<GF2, AllModules<GF2> > _GJ;

	DenseMatrix<bool> _L;
	GaussJordan<GF2, AllModules<GF2> >::Permutation _P;

	// Map pointers to matrices to computed ranks
	std::map<const void *, size_t> _rank_table;

public:
	enum Method { METHOD_STANDARD_GJ, METHOD_ASYMPTOTICALLY_FAST_GJ, METHOD_M4RI };

	EchelonForm (Context<GF2, AllModules<GF2> > &ctx) : _ctx (ctx), _GJ (ctx) {}

	template <class Matrix>
	Matrix &RowEchelonForm (Matrix &A, bool reduced = false, Method method = METHOD_STANDARD_GJ)
	{
		static const char *method_names[] = { "standard", "recursive", "M4RI" };

		std::ostringstream str;
		str << "Row-echelon form (method: " << method_names[method] << ")" << std::ends;

		commentator.start (str.str ().c_str (), __FUNCTION__);

		size_t rank;
		bool d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_GJ.StandardRowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		_rank_table[&A] = rank;

		commentator.stop (MSG_DONE);

		return A;
	}

	// Specialisation for M4RI-matrices
	DenseMatrix<bool> &RowEchelonForm (DenseMatrix<bool> &A, bool reduced = false, Method method = METHOD_M4RI)
	{
		static const char *method_names[] = { "standard", "recursive", "M4RI" };

		std::ostringstream str;
		str << "Row-echelon form (method: " << method_names[method] << ")" << std::ends;

		commentator.start (str.str ().c_str (), __FUNCTION__);

		size_t rank;
		bool d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_GJ.StandardRowEchelonForm (A, _L, _P, rank, d, reduced, false);
			_rank_table[&A] = rank;
			break;

		case METHOD_ASYMPTOTICALLY_FAST_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_GJ.DenseRowEchelonForm (A, _L, _P, rank, d, reduced);
			_rank_table[&A] = rank;
			break;

		case METHOD_M4RI:
			mzd_echelonize_pluq (A._rep, reduced ? 1 : 0);
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		commentator.stop (MSG_DONE);

		return A;
	}

	template <class Matrix>
	size_t rank (const Matrix &A) const
		{ return _rank_table[&A]; }

	size_t rank (const DenseMatrix<bool> &A)
	{
		if (_rank_table.find (&A) != _rank_table.end ())
			return _rank_table[&A];
		else {
			// Not computed yet, so we compute it ourselves. Assume matrix is in row-echelon form.
			if (BLAS3::is_zero (_ctx, A)) {
				_rank_table[&A] = 0;
				return 0;
			}

			DenseMatrix<bool>::ConstRowIterator i = A.rowBegin () + (A.rowdim () - 1);

			size_t r = A.rowdim ();

			while (BLAS1::is_zero (_ctx, *i)) {
				--i; --r;
			}

			_rank_table[&A] = r;

			return r;
		}
	}
};

} // namespace LinBox

#endif // __LINBOX_HAVE_M4RI

#endif // __LINBOX_SOLUTIONS_ECHELON_FORM_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
