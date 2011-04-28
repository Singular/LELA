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

#include "linbox/solutions/echelon-form.h"
#include "linbox/field/gf2.h"

// No point specialising if we aren't using libm4ri
#ifdef __LINBOX_HAVE_M4RI

#include <m4ri/m4ri.h>

#include "linbox/matrix/m4ri-matrix.h"

namespace LinBox
{

template <>
class EchelonForm<GF2>
{
	GaussJordan<GF2> _GJ;

	DenseMatrix<bool> _L;
	MatrixDomain<GF2>::Permutation _P;

public:
	enum Method { METHOD_STANDARD_GJ, METHOD_ASYMPTOTICALLY_FAST_GJ, METHOD_M4RI };

	EchelonForm (const GF2 &F) : _GJ (F) {}

	template <class Matrix>
	Matrix &RowEchelonForm (Matrix &A, bool reduced = false, Method method = METHOD_STANDARD_GJ)
	{
		size_t rank;
		bool d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_GJ.StandardRowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		return A;
	}

	// Specialisation for M4RI-matrices
	M4RIMatrix &RowEchelonForm (M4RIMatrix &A, bool reduced = false, Method method = METHOD_M4RI)
	{
		size_t rank;
		bool d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			_GJ.StandardRowEchelonForm (A, _L, _P, rank, d, reduced, false);
			break;

		case METHOD_ASYMPTOTICALLY_FAST_GJ:
			_L.resize (A.rowdim (), A.rowdim ());
			_GJ.DenseRowEchelonForm (A, _L, _P, rank, d, reduced);
			break;

		case METHOD_M4RI:
			mzd_echelonize_pluq (A._rep, reduced ? 1 : 0);
			break;

		default:
			throw LinboxError ("Invalid method for choice of matrix");
		}

		return A;
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
