/* lela/algorithms/gauss-jordan.tcc
 * Copyright 2010, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Gauss-Jordan elimination
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ALGORITHMS_GAUSS_JORDAN_TCC
#define __LELA_ALGORITHMS_GAUSS_JORDAN_TCC

#include "lela/algorithms/gauss-jordan.h"

#include "lela/blas/level1.h"
#include "lela/blas/level3.h"

namespace LELA
{

template <class Ring, class Modules>
template <class Matrix1, class Matrix2, class PivotStrategy>
void GaussJordan<Ring, Modules>::GaussJordanTransform (Matrix1      &A,
						       int           k,
						       Element       d_0,
						       Matrix2      &U,
						       Permutation  &P,
						       size_t       &r,
						       int          &h,
						       Element      &d,
						       DenseMatrix<typename Ring::Element> &S,
						       DenseMatrix<typename Ring::Element> &T,
						       PivotStrategy PS) const
{
	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	// DEBUG
	// report << "enter" << std::endl;
	// report << "A =" << std::endl;
	// BLAS3::write (ctx, report, A);
	// report << "U =" << std::endl;
	// BLAS3::write (ctx, report, U);
	// report << "k = " << k << ", d_0 = " << d_0 << std::endl;

	typename Matrix1::SubmatrixType Aw (A, k, 0, A.rowdim () - k, A.coldim ());

	if (BLAS3::is_zero (ctx, Aw)) {
		// DEBUG
		// report << "A is 0" << std::endl;

		r = 0;
		h = A.rowdim () - k;
		ctx.F.copy (d, d_0);
	}
	else if (A.coldim () == 1) {
		// DEBUG
		// report << "A.coldim == 1" << std::endl;

		size_t l = P.size (), j, row = k, col = 0;

		typename Ring::Element aii, negaiiinv, aij;

		// Find the pivot
		if (!PS.getPivot (A, aii, row, col))
			throw LELAError ("Could not find pivot even though matrix is reported as nonzero");

		lela_check (col == 0);

		// Set output variable d
		ctx.F.mulin (d, aii);

		// Prepare pivot for elimination
		if (!ctx.F.inv (negaiiinv, aii))
			throw LELAError ("Pivot not invertible in the ring");

		ctx.F.negin (negaiiinv);

		if (row != (size_t) k) {
			// Set the permutation
			P.push_back (Transposition (k, row));

			// Swap the pivot-entry and the entry at (0,0)
			if (A.getEntry (aij, k, 0)) {
				A.setEntry (row, 0, aij);
			} else {
				A.setEntry (row, 0, ctx.F.zero ());
				A.eraseEntry (row, 0);
			}

			A.setEntry (k, 0, aii);
		}

		// Perform elimination-steps to complete L
		for (j = 0; j < A.rowdim (); ++j) {
			if (j != (size_t) k && A.getEntry (aij, j, 0)) {
				A.setEntry (j, 0, ctx.F.zero ());
				A.eraseEntry (j, 0);
				ctx.F.mulin (aij, negaiiinv);
				U.setEntry (j, k, aij);
			}
		}

		// Set output-variable r
		r = 1;

		// DEBUG
		// report << "A after elimination:" << std::endl;
		// BLAS3::write (ctx, report, A);

		if (!P.empty ())
			h = std::min_element (P.begin () + l, P.end (), CompareSecond ())->second;
		else
			h = k;
	}
	else {
		// DEBUG
		// report << "A.coldim () > 1" << std::endl;

		int m_1 = A.coldim () / 2, m_2 = A.coldim () - m_1;
		typename Matrix1::SubmatrixType A_1 (A, 0, 0,   A.rowdim (), m_1);
		typename Matrix1::SubmatrixType A_2 (A, 0, m_1, A.rowdim (), m_2);

		size_t r_1, r_2;
		int h_1, h_2;
		Element d_1;

		ctx.F.copy (d_1, ctx.F.one ());

		GaussJordanTransform (A_1, k, d_0, U, P, r_1, h_1, d_1, S, T, PS);

		// DEBUG
		// report << "Status after first recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "P = ";
		// BLAS3::write_permutation (report, P) << std::endl;
		// report << "R =" << std::endl;
		// BLAS3::write (ctx, report, R);
		// report << "r_1 = " << r_1 << std::endl;
		// report << "d_1 = " << d_1 << std::endl;

		typename Matrix2::SubmatrixType U_2  (U, 0,       k, U.rowdim (),           r_1);
		typename Matrix2::SubmatrixType U_21 (U, 0,       k, k,                     r_1);
		typename Matrix2::SubmatrixType U_22 (U, k,       k, r_1,                   r_1);
		typename Matrix2::SubmatrixType U_23 (U, r_1 + k, k, U.rowdim () - r_1 - k, r_1);

		typename Matrix1::SubmatrixType A_21 (A, 0,       m_1, k,                     m_2);
		typename Matrix1::SubmatrixType A_22 (A, k,       m_1, r_1,                   m_2);
		typename Matrix1::SubmatrixType A_23 (A, k + r_1, m_1, A.rowdim () - r_1 - k, m_2);

		BLAS3::permute_rows (ctx, P.begin (), P.end (), A_2);

		T.resize (r_1, m_2);
		BLAS3::copy (ctx, A_22, T);
		BLAS3::gemm (ctx, ctx.F.one (), U_21, T, ctx.F.one (),  A_21);
		BLAS3::gemm (ctx, ctx.F.one (), U_22, T, ctx.F.zero (), A_22);
		BLAS3::gemm (ctx, ctx.F.one (), U_23, T, ctx.F.one (),  A_23);

		Permutation P_2;

		// DEBUG
		// report << "Status before second recursive call:" << std::endl;
		// report << "B =" << std::endl;
		// BLAS3::write (ctx, report, R_2);

		GaussJordanTransform (A_2, k + r_1, d_1, U, P_2, r_2, h_2, d, S, T, PS);

		// DEBUG
		// report << "Status after second recursive call:" << std::endl;
		// report << "U =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "R =" << std::endl;
		// BLAS3::write (ctx, report, R);
		// report << "r_2 = " << r_2 << std::endl;
		// report << "d = " << d << std::endl;

		typename Matrix2::SubmatrixType P_2U_2  (U,  0,             k,       U.rowdim (),                  r_1);
		typename Matrix2::SubmatrixType U_212   (U,  0,             k,       k + r_1,                      r_1);
		typename Matrix2::SubmatrixType P_2U_23 (U,  k + r_1,       k,       r_2,                          r_1);
		typename Matrix2::SubmatrixType P_2U_24 (U,  k + r_1 + r_2, k,       U.rowdim () - r_1 - r_2 - k,  r_1);
		typename Matrix2::SubmatrixType U_312   (U,  0,             k + r_1, k + r_1,                      r_2);
		typename Matrix2::SubmatrixType U_33    (U,  k + r_1,       k + r_1, r_2,                          r_2);
		typename Matrix2::SubmatrixType U_34    (U,  k + r_1 + r_2, k + r_1, U.rowdim () - r_1 - r_2 - k,  r_2);

		// DEBUG
		// report << "U_2 =" << std::endl;
		// BLAS3::write (ctx, report, U);
		// report << "P_2 = ";
		// BLAS3::write_permutation (report, P_2) << std::endl;

		BLAS3::permute_rows (ctx, P_2.begin (), P_2.end (), P_2U_2);
				
		// DEBUG
		// report << "P_2U_2 =" << std::endl;
		// BLAS3::write (ctx, report, P_2U_2);
		// report << "P_2U_23 =" << std::endl;
		// BLAS3::write (ctx, report, P_2U_23);
		// report << "U_3 =" << std::endl;;
		// BLAS3::write (ctx, report, U_3);

		S.resize (r_2, r_1);
		BLAS3::copy (ctx, P_2U_23, S);

		BLAS3::gemm (ctx, ctx.F.one (), U_312, S, ctx.F.one (), U_212);
		BLAS3::gemm (ctx, ctx.F.one (), U_33, S, ctx.F.zero (), P_2U_23);
		BLAS3::gemm (ctx, ctx.F.one (), U_34, S, ctx.F.one (), P_2U_24);

		// DEBUG
		// report << "U_3P_2U_23 =" << std::endl;
		// BLAS3::write (ctx, report, U_3P_2U_23);

		P.insert (P.end (), P_2.begin (), P_2.end ());
		r = r_1 + r_2;
		h = std::min (h_1, h_2);
	}

	// DEBUG
	// report << "Status at end:" << std::endl;
	// report << "U =" << std::endl;
	// BLAS3::write (ctx, report, U);
	// report << "P = ";
	// BLAS3::write_permutation (report, P) << std::endl;
	// report << "R =" << std::endl;
	// BLAS3::write (ctx, report, R);
	// report << "k = " << k << ", r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
}

template <class Ring, class Modules>
template <class Matrix1, class Matrix2, class PivotStrategy>
void GaussJordan<Ring, Modules>::GaussTransform (Matrix1     &A,
						 Matrix2     &L,
						 Element      d_0,
						 Permutation &P,
						 size_t      &r,
						 int         &h,
						 Element     &d,
						 PivotStrategy PS) const
{
	// std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	// DEBUG
	// report << "enter" << std::endl;
	// report << "A =" << std::endl;
	// BLAS3::write (ctx, report, A);

	if (BLAS3::is_zero (ctx, A)) {
		// DEBUG
		// report << "A is 0" << std::endl;

		r = 0;
		h = A.rowdim ();
		ctx.F.copy (d, d_0);
	}
	else if (A.coldim () == 1) {
		// DEBUG
		// report << "A.coldim == 1" << std::endl;

		size_t l = P.size (), j, row = 0, col = 0;

		typename Ring::Element aii, negaiiinv, aij;

		// Find the pivot
		if (!PS.getPivot (A, aii, row, col))
			throw LELAError ("Could not find pivot even though matrix is reported as nonzero");

		lela_check (col == 0);

		// Set output variable d
		ctx.F.mulin (d, aii);

		// Prepare pivot for elimination
		if (!ctx.F.inv (negaiiinv, aii))
			throw LELAError ("Pivot not invertible in the ring");

		ctx.F.negin (negaiiinv);

		if (row != 0) {
			// Set the permutation
			P.push_back (Transposition (0, row));

			// Swap the pivot-entry and the entry at (0,0)
			if (A.getEntry (aij, 0, 0)) {
				A.setEntry (row, 0, aij);
			} else {
				A.setEntry (row, 0, ctx.F.zero ());
				A.eraseEntry (row, 0);
			}
			    
			A.setEntry (0, 0, aii);
		}

		// Perform elimination-steps to complete L
		for (j = 1; j < A.rowdim (); ++j) {
			if (A.getEntry (aij, j, 0)) {
				A.setEntry (j, 0, ctx.F.zero ());
				A.eraseEntry (j, 0);
				ctx.F.mulin (aij, negaiiinv);
				L.setEntry (j, 0, aij);
			}
		}

		// Set output-variable r
		r = 1;

		// DEBUG
		// report << "A after elimination:" << std::endl;
		// BLAS3::write (ctx, report, A);

		// Set output-variable h
		if (!P.empty ())
			h = std::min_element (P.begin () + l, P.end (), CompareSecond ())->second;
		else
			h = 0;
	}
	else {
		// DEBUG
		// report << "A.coldim () > 1" << std::endl;

		int m_1 = A.coldim () / 2, m_2 = A.coldim () - m_1;
		typename Matrix1::SubmatrixType A_1 (A, 0, 0, A.rowdim (), m_1);

		size_t r_1, r_2;
		int h_1, h_2;
		Element d_1;

		ctx.F.copy (d_1, ctx.F.one ());

		GaussTransform (A_1, L, d_0, P, r_1, h_1, d_1, PS);

		// DEBUG
		// report << "Status after first recursive call:" << std::endl;
		// report << "A =" << std::endl;
		// BLAS3::write (ctx, report, A);
		// report << "P = ";
		// BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;
		// report << "r_1 = " << r_1 << std::endl;
		// report << "d_1 = " << d_1 << std::endl;

		typename Matrix2::SubmatrixType L_11 (L, 0,   0,   r_1,               r_1);
		typename Matrix2::SubmatrixType L_12 (L, r_1, 0,   L.rowdim () - r_1, r_1);
		typename Matrix2::SubmatrixType L_22 (L, r_1, r_1, L.rowdim () - r_1, L.coldim () - r_1);

		typename Matrix1::SubmatrixType A_2  (A, 0,   m_1, A.rowdim (),       m_2);
		typename Matrix1::SubmatrixType A_21 (A, 0,   m_1, r_1,               m_2);
		typename Matrix1::SubmatrixType A_22 (A, r_1, m_1, A.rowdim () - r_1, m_2);

		BLAS3::permute_rows (ctx, P.begin (), P.end (), A_2);
		BLAS3::gemm (ctx, ctx.F.one (), L_12, A_21, ctx.F.one (), A_22);
		BLAS3::trmm (ctx, ctx.F.one (), L_11, A_21, LowerTriangular, true);

		Permutation P_2;

		// DEBUG
		// report << "Status before second recursive call:" << std::endl;
		// report << "A_22 =" << std::endl;
		// BLAS3::write (ctx, report, A_22);

		GaussTransform (A_22, L_22, d_1, P_2, r_2, h_2, d, PS);

		// DEBUG
		// report << "Status after second recursive call:" << std::endl;
		// report << "A =" << std::endl;
		// BLAS3::write (ctx, report, A);
		// report << "r_2 = " << r_2 << std::endl;
		// report << "d = " << d << std::endl;

		// DEBUG
		// report << "P_2 = ";
		// BLAS1::write_permutation (report, P_2.begin (), P_2.end ()) << std::endl;

		BLAS3::permute_rows (ctx, P_2.begin (), P_2.end (), L_12);

		typename Matrix2::SubmatrixType L_12p (L, r_1, 0, r_2, r_1);
		typename Matrix2::SubmatrixType L_13p (L, r_1 + r_2, 0, L.rowdim () - (r_1 + r_2), r_1);
		typename Matrix2::SubmatrixType L_22p (L, r_1, r_1, r_2, r_2);
		typename Matrix2::SubmatrixType L_23p (L, r_1 + r_2, r_1, L.rowdim () - (r_1 + r_2), r_2);

		BLAS3::gemm (ctx, ctx.F.one (), L_23p, L_12p, ctx.F.one (), L_13p);
		BLAS3::trmm (ctx, ctx.F.one (), L_22p, L_12p, LowerTriangular, true);

		// DEBUG
		// report << "P_2L_12 =" << std::endl;
		// BLAS3::write (ctx, report, P_2L_12);

		size_t P_len = P.size ();
		P.insert (P.end (), P_2.begin (), P_2.end ());

		// Update indices in permutation to refer to whole matrix
		typename Permutation::iterator i;

		for (i = P.begin () + P_len; i != P.end (); ++i) {
			i->first += r_1;
			i->second += r_1;
		}

		r = r_1 + r_2;
		h = std::min (h_1, h_2);
	}

	// DEBUG
	// report << "Status at end:" << std::endl;
	// report << "P = ";
	// BLAS1::write_permutation (report, P.begin (), P.end ()) << std::endl;
	// report << "A =" << std::endl;
	// BLAS3::write (ctx, report, A);
	// report << "r = " << r << ", d_0 = " << d_0 << ", d = " << d << std::endl;
}

template <class Ring, class Modules>
template <class Matrix, class PivotStrategy>
Matrix &GaussJordan<Ring, Modules>::echelonize (Matrix      &A,
						Permutation &P,
						size_t      &rank,
						Element     &det,
						PivotStrategy PS)
{
	commentator.start ("Asymptotically fast row-echelon form", __FUNCTION__);

	int h;

	ctx.F.copy (det, ctx.F.one ());

	GaussTransform (A, A, ctx.F.one (), P, rank, h, det, PS);

	commentator.stop (MSG_DONE);

	return A;
}

template <class Ring, class Modules>
template <class Matrix1, class Matrix2, class PivotStrategy>
Matrix1 &GaussJordan<Ring, Modules>::echelonize_reduced (Matrix1       &A,
							 Matrix2       &L,
							 Permutation   &P,
							 size_t        &rank,
							 Element       &det,
							 PivotStrategy  PS)
{
	lela_check (L.rowdim () == A.rowdim ());
	lela_check (L.coldim () == A.rowdim ());

	commentator.start ("Asymptotically fast reduced row-echelon form", __FUNCTION__);

	int h;

	ctx.F.copy (det, ctx.F.one ());

	DenseMatrix<typename Ring::Element> S, T;

	// Set L first to the identity-matrix
	StandardBasisStream<Ring, typename Matrix2::Row> s (ctx.F, L.rowdim ());
	typename Matrix2::RowIterator i_L;

	for (i_L = L.rowBegin (); i_L != L.rowEnd (); ++i_L)
		s >> *i_L;

	GaussJordanTransform (A, 0, ctx.F.one (), L, P, rank, h, det, S, T, PS);

	commentator.stop (MSG_DONE);

	return A;
}

} // namespace LELA

#endif // __LELA_ALGORITHMS_GAUSS_JORDAN_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
