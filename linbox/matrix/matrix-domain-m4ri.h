/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/matrix-domain-m4ri.h
 * Copyright 2010 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __MATRIX_MATRIX_DOMAIN_M4RI_H
#define __MATRIX_MATRIX_DOMAIN_M4RI_H

#include "linbox/linbox-config.h"

#ifndef __LINBOX_HAVE_M4RI
#  error "This header file requires that LinBox be configured with libm4ri enabled. Please ensure that libm4ri is properly installed and re-run configure."
#endif

#include <m4ri/m4ri.h>

#include <iostream>
#include <time.h>

#include "linbox/field/gf2.h"
#include "linbox/matrix/matrix-domain-gf2.h"
#include "linbox/matrix/m4ri-matrix.h"

namespace LinBox 
{

/** Wrapper for libm4ri
 *
 * Provides the parts of the interface of \ref MatrixDomainSupport
 * which are supported by libm4ri.
 */
class MatrixDomainM4RI : public MatrixDomainSupportGF2
{
    public:
	MatrixDomainM4RI (const GF2 &F) : MatrixDomainSupportGF2 (F) {}

	template <class Matrix1, class Matrix2>
	inline Matrix1 &copy (Matrix1 &B, const Matrix2 &A) const
		{ return MatrixDomainSupportGF2::copy (B, A); }

	inline M4RIMatrix &copy (M4RIMatrix &B, const M4RIMatrix &A) const
		{ mzd_copy (B._rep, A._rep); return B; }

	inline Submatrix<M4RIMatrix> &copy (Submatrix<M4RIMatrix> &B, const M4RIMatrix &A) const
		{ mzd_copy (B._rep, A._rep); return B; }

	inline M4RIMatrix &copy (M4RIMatrix &B, const Submatrix<M4RIMatrix> &A) const
		{ mzd_copy (B._rep, A._rep); return B; }

	inline Submatrix<M4RIMatrix> &copy (Submatrix<M4RIMatrix> &B, const Submatrix<M4RIMatrix> &A) const
		{ mzd_copy (B._rep, A._rep); return B; }

	template <class Matrix1, class Matrix2>
	inline bool areEqual (const Matrix1 &A, const Matrix2 &B) const
		{ return MatrixDomainSupportGF2::areEqual (A, B); }

	inline bool areEqual (const M4RIMatrix &A, const M4RIMatrix &B) const
		{ return mzd_equal (A._rep, B._rep); }

 	inline bool areEqual (const Submatrix<M4RIMatrix> &A, const M4RIMatrix &B) const
		{ return mzd_equal (A._rep, B._rep); }

 	inline bool areEqual (const M4RIMatrix &A, const Submatrix<M4RIMatrix> &B) const
		{ return mzd_equal (A._rep, B._rep); }

 	inline bool areEqual (const Submatrix<M4RIMatrix> &A, const Submatrix<M4RIMatrix> &B) const
		{ return mzd_equal (A._rep, B._rep); }

	template <class Matrix>
	inline bool isZero (const Matrix &A) const
		{ return MatrixDomainSupportGF2::isZero (A); }

	inline bool isZero (const M4RIMatrix &A) const
		{ return mzd_is_zero (A._rep); }

	inline bool isZero (const Submatrix<M4RIMatrix> &A) const
		{ return mzd_is_zero (A._rep); }

	template <class Matrix>
	inline Matrix &scal (Matrix &A, const bool &a) const
		{ return MatrixDomainSupportGF2::scal (A, a); }

	inline M4RIMatrix &scal (M4RIMatrix &A, const bool &a) const
	{
		size_t i;

		if (!a)
			for (i = 0; i < A.rowdim (); ++i)
				mzd_row_clear_offset (A._rep, i, 0);
		return A;
	}

	inline Submatrix<M4RIMatrix> &scal (Submatrix<M4RIMatrix> &A, const bool &a) const
		{ scal ((M4RIMatrix &) A, a); return A; }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &axpy (const bool &a, const Matrix1 &A, Matrix2 &B) const
		{ return MatrixDomainSupportGF2::axpy (a, A, B); }

	inline M4RIMatrix &axpy (const bool &a, const M4RIMatrix &A, M4RIMatrix &B) const
	{
		if (a)
			mzd_add (B._rep, B._rep, A._rep);

		return B;
	}

	inline M4RIMatrix &axpy (const bool &a, const Submatrix<M4RIMatrix> &A, M4RIMatrix &B) const
		{ axpy (a, (M4RIMatrix &) A, B); return B; }

	inline M4RIMatrix &axpy (const bool &a, const M4RIMatrix &A, Submatrix<M4RIMatrix> &B) const
		{ axpy (a, A, (M4RIMatrix &) B); return B; }

	inline Submatrix<M4RIMatrix> &axpy (const bool &a, const Submatrix<M4RIMatrix> &A, Submatrix<M4RIMatrix> &B) const
		{ axpy (a, (M4RIMatrix &) A, (M4RIMatrix &) B); return B; }

	template <class Matrix1, class Matrix2, class Matrix3>
	inline Matrix3 &gemm (const bool &a, const Matrix1 &A, const Matrix2 &B, const bool &b, Matrix3 &C) const
		{ return MatrixDomainSupportGF2::gemm (a, A, B, b, C); }

	inline M4RIMatrix &gemm (const bool &a, const M4RIMatrix &A, const M4RIMatrix &B, const bool &b, M4RIMatrix &C) const
	{
		if (a) {
			if (b)
				mzd_addmul (C._rep, A._rep, B._rep, STRASSEN_MUL_CUTOFF);
			else
				mzd_mul (C._rep, A._rep, B._rep, STRASSEN_MUL_CUTOFF);
		}
		else if (!b)
			scal (C, false);

		return C;
	}

	inline M4RIMatrix &gemm (const bool &a, const Submatrix<M4RIMatrix> &A, const M4RIMatrix &B, const bool &b, M4RIMatrix &C) const
		{ gemm (a, (const M4RIMatrix &) A, B, b, C); return C; }

	inline M4RIMatrix &gemm (const bool &a, const M4RIMatrix &A, const Submatrix<M4RIMatrix> &B, const bool &b, M4RIMatrix &C) const
		{ gemm (a, A, (const M4RIMatrix &) B, b, C); return C; }

	inline M4RIMatrix &gemm (const bool &a, const Submatrix<M4RIMatrix> &A, const Submatrix<M4RIMatrix> &B, const bool &b, M4RIMatrix &C) const
		{ gemm (a, (const M4RIMatrix &) A, (const M4RIMatrix &) B, b, C); return C; }

	inline Submatrix<M4RIMatrix> &gemm (const bool &a, const M4RIMatrix &A, const M4RIMatrix &B, const bool &b, Submatrix<M4RIMatrix> &C) const
		{ gemm (a, A, B, b, (M4RIMatrix &) C); return C; }

	inline Submatrix<M4RIMatrix> &gemm (const bool &a, const Submatrix<M4RIMatrix> &A, const M4RIMatrix &B, const bool &b, Submatrix<M4RIMatrix> &C) const
		{ gemm (a, (const M4RIMatrix &) A, B, b, (M4RIMatrix &) C); return C; }

	inline Submatrix<M4RIMatrix> &gemm (const bool &a, const M4RIMatrix &A, const Submatrix<M4RIMatrix> &B, const bool &b, Submatrix<M4RIMatrix> &C) const
		{ gemm (a, A, (const M4RIMatrix &) B, b, (M4RIMatrix &) C); return C; }

	inline Submatrix<M4RIMatrix> &gemm (const bool &a, const Submatrix<M4RIMatrix> &A, const Submatrix<M4RIMatrix> &B, const bool &b, Submatrix<M4RIMatrix> &C) const
		{ gemm (a, (const M4RIMatrix &) A, (const M4RIMatrix &) B, b, (M4RIMatrix &) C); return C; }

	template <class Matrix1, class Matrix2>
	inline Matrix2 &trsm (const bool &a, const Matrix1 &A, Matrix2 &B) const
		{ return MatrixDomainSupportGF2::trsm (a, A, B); }

	inline M4RIMatrix &trsm (const bool &a, const M4RIMatrix &A, M4RIMatrix &B) const
	{
		if (a)
			mzd_trsm_upper_right (A._rep, B._rep, STRASSEN_MUL_CUTOFF);
		else
			scal (B, false);

		return B;
	}

	inline M4RIMatrix &trsm (const bool &a, const Submatrix<M4RIMatrix> &A, M4RIMatrix &B) const
		{ trsm (a, (const M4RIMatrix &) A, B); return B; }

	inline Submatrix<M4RIMatrix> &trsm (const bool &a, const M4RIMatrix &A, Submatrix<M4RIMatrix> &B) const
		{ trsm (a, A, (M4RIMatrix &) B); return B; }

	inline Submatrix<M4RIMatrix> &trsm (const bool &a, const Submatrix<M4RIMatrix> &A, Submatrix<M4RIMatrix> &B) const
		{ trsm (a, (const M4RIMatrix &) A, (M4RIMatrix &) B); return B; }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRows (Matrix   &A,
				    Iterator  P_start,
				    Iterator  P_end) const
		{ return MatrixDomainSupportGF2::permuteRows (A, P_start, P_end); }

#if 0 // Not working, disabled
	template <class Iterator>
	inline M4RIMatrix &permuteRows (M4RIMatrix &A,
					Iterator    P_start,
					Iterator    P_end) const
	{
		mzp_t *P = make_m4ri_permutation (P_start, P_end, A.rowdim ());
		mzd_apply_p_left (A._rep, P);
		mzp_free (P);
		return A;
	}

	template <class Iterator>
	inline Submatrix<M4RIMatrix> &permuteRows (Submatrix<M4RIMatrix> &A,
						   Iterator    P_start,
						   Iterator    P_end) const
		{ permuteRows ((M4RIMatrix &) A, P_start, P_end); return A; }
#endif // Not working, disabled

	template <class Matrix, class Iterator>
	inline Matrix &permuteColumns (Matrix   &A,
				       Iterator  P_start,
				       Iterator  P_end) const
		{ return MatrixDomainSupportGF2::permuteColumns (A, P_start, P_end); }

#if 0 // Not working, disabled
	template <class Iterator>
	inline M4RIMatrix &permuteColumns (M4RIMatrix &A,
					   Iterator    P_start,
					   Iterator    P_end) const
	{
		mzp_t *P = make_m4ri_permutation (P_start, P_end, A.coldim ());
		mzd_apply_p_right (A._rep, P);
		mzp_free (P);
		return A;
	}

	template <class Iterator>
	inline Submatrix<M4RIMatrix> &permuteColumns (Submatrix<M4RIMatrix> &A,
						      Iterator    P_start,
						      Iterator    P_end) const
		{ permuteColumns ((M4RIMatrix &) A, P_start, P_end); return A; }
#endif // Not working, disabled

    private:
	// Convert a partition in MatrixDomain-format to LAPACK-format for use in libm4ri
	template <class Iterator>
	mzp_t *make_m4ri_permutation (Iterator P_start, Iterator P_end, size_t len) const {
		mzp_t *P;
		Iterator i_P;

		P = mzp_init (len);

		for (i_P = P_start; i_P != P_end; ++i_P) {
			linbox_check (i_P->first < len);
			linbox_check (i_P->second < len);
			P->values[i_P->first] = P->values[i_P->second];
		}

		return P;
	}
};

} // namespace LinBox

#endif // __MATRIX_MATRIX_DOMAIN_M4RI_H
