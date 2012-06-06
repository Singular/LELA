/* lela/solutions/echelon-form-gf2.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Compute the echelon-form of a matrix
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_SOLUTIONS_ECHELON_FORM_GF2_H
#define __LELA_SOLUTIONS_ECHELON_FORM_GF2_H

// No point specialising if we aren't using libm4ri
#ifdef __LELA_HAVE_M4RI

#include <m4ri/m4ri.h>

#include "lela/util/commentator.h"
#include "lela/ring/gf2.h"
#include "lela/solutions/echelon-form.h"
#include "lela/algorithms/elimination.h"
#include "lela/algorithms/gauss-jordan.h"
#include "lela/algorithms/faugere-lachartre.h"
#include "lela/matrix/m4ri-matrix.h"
#include "lela/blas/level1.h"
#include "lela/blas/level3.h"

namespace LELA
{

// Specialisation of EchelonForm to GF2 to take advantage of M4RI-routines
template <>
class EchelonForm<GF2, AllModules<GF2> >
{
	Context<GF2, AllModules<GF2> > &_ctx;
	Elimination<GF2, AllModules<GF2> > _elim;
	GaussJordan<GF2, AllModules<GF2> > _GJ;

	DenseMatrix<bool> _L;
	GaussJordan<GF2, AllModules<GF2> >::Permutation _P;

	// Map pointers to matrices to computed ranks
	std::map<const void *, size_t> _rank_table;

public:
	enum Method { METHOD_UNKNOWN, METHOD_STANDARD_GJ, METHOD_ASYMPTOTICALLY_FAST_GJ, METHOD_M4RI, METHOD_FAUGERE_LACHARTRE };

	EchelonForm (Context<GF2, AllModules<GF2> > &ctx) : _ctx (ctx), _elim (ctx), _GJ (ctx) {}

	template <class Matrix>
	Matrix &echelonize (Matrix &A, bool reduced = false, Method method = METHOD_STANDARD_GJ)
	{
		static const char *method_names[] = { "unknown", "standard", "recursive", "M4RI", "Faugère-Lachartre" };

		std::ostringstream str;
		str << "Row-echelon form (method: " << method_names[method] << ")" << std::ends;

		commentator.start (str.str ().c_str (), __FUNCTION__);

		size_t rank;
		bool d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			if (reduced)
				_elim.echelonize_reduced (A, _L, _P, rank, d, false);
			else
				_elim.echelonize (A, _P, rank, d, false);

			break;

		case METHOD_FAUGERE_LACHARTRE:
			{
				// Must do it this way to avoid an infinite loop of inclusion...
				FaugereLachartre<GF2, AllModules<GF2> > FL (_ctx);
				FL.echelonize (A, A, rank, d, reduced);
			}

			break;

		default:
			throw LELAError ("Invalid method for choice of matrix");
		}

		_rank_table[&A] = rank;

		commentator.stop (MSG_DONE);

		return A;
	}

	// Specialisation for M4RI-matrices
	DenseMatrix<bool> &echelonize (DenseMatrix<bool> &A, bool reduced = false, Method method = METHOD_M4RI)
	{
		static const char *method_names[] = { "unknown", "standard", "recursive", "M4RI", "Faugère-Lachartre" };

		std::ostringstream str;
		str << "Row-echelon form (method: " << method_names[method] << ")" << std::ends;

		commentator.start (str.str ().c_str (), __FUNCTION__);

		size_t rank;
		bool d;

		switch (method) {
		case METHOD_STANDARD_GJ:
			if (reduced)
				_elim.echelonize_reduced (A, _L, _P, rank, d, false);
			else
				_elim.echelonize (A, _P, rank, d, false);

			_rank_table[&A] = rank;
			break;

		case METHOD_ASYMPTOTICALLY_FAST_GJ:
			_L.resize (A.rowdim (), A.rowdim ());

			if (reduced)
				_GJ.echelonize_reduced (A, _L, _P, rank, d);
			else {
				_GJ.echelonize (A, _P, rank, d);
				_elim.move_L (A, A);
			}

			_rank_table[&A] = rank;
			break;

		case METHOD_M4RI:
			mzd_echelonize (A._rep, reduced ? 1 : 0);
			break;

		case METHOD_FAUGERE_LACHARTRE:
			{
				// Must do it this way to avoid an infinite loop of inclusion...
				FaugereLachartre<GF2, AllModules<GF2> > FL (_ctx);
				FL.echelonize (A, A, rank, d, reduced);
			}

			break;

		default:
			throw LELAError ("Invalid method for choice of matrix");
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

} // namespace LELA

#endif // __LELA_HAVE_M4RI

#endif // __LELA_SOLUTIONS_ECHELON_FORM_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
