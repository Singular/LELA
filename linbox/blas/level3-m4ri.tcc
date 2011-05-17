/* linbox/blas/level3-m4ri.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_M4RI_TCC
#define __BLAS_LEVEL3_M4RI_TCC

#include "linbox/blas/level3-m4ri.h"

namespace LinBox
{

namespace BLAS3
{

DenseMatrix<bool> &scal_impl (const GF2 &F, M4RIModule &M, bool a, DenseMatrix<bool> &A, MatrixCategories::RowMatrixTag)
{
	size_t i;

	if (A.rowdim () == 0 || A.coldim () == 0)
		return A;

	if (!a)
		for (i = 0; i < A.rowdim (); ++i)
			mzd_row_clear_offset (A._rep, i, 0);
	return A;
}


DenseMatrix<bool> &axpy_impl (const GF2 &F, M4RIModule &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	if (A.rowdim () == 0 || A.coldim () == 0)
		return B;

	if (a)
		mzd_add (B._rep, B._rep, A._rep);

	return B;
}

DenseMatrix<bool> &gemm_impl (const GF2 &F, M4RIModule &M,
			      bool a, const DenseMatrix<bool> &A, const DenseMatrix<bool> &B, bool b, DenseMatrix<bool> &C,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	if (A.rowdim () == 0 || B.rowdim () == 0 || B.coldim () == 0)
		return C;

	if (a) {
		if (b)
			mzd_addmul (C._rep, A._rep, B._rep, STRASSEN_MUL_CUTOFF);
		else
			mzd_mul (C._rep, A._rep, B._rep, STRASSEN_MUL_CUTOFF);
	}
	else if (!b)
		_scal (F, M, false, C);

	return C;
}

DenseMatrix<bool> &trmm_impl (const GF2 &F, M4RIModule &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B, TriangularMatrixType type, bool diagIsOne,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.rowdim () == A.coldim ());

	// Our method doesn't work if the diagonal isn't one, so just fall back to the generic version in that case
	if (!diagIsOne) {
		GenericModule M1 (M);
		return trmm_impl (F, M1, a, A, B, type, false, MatrixCategories::RowMatrixTag (), MatrixCategories::RowMatrixTag ());
	}

	if (A.rowdim () == 0 || A.coldim () == 0 || B.coldim () == 0)
		return B;

	if (a) {
		// Since M4RI doesn't seem to have in-place
		// trmm, we'll do it by two calls to
		// trsm. First we construct an
		// identity-matrix, then we use trsm to get
		// the inverse of A, then we use trsm again to
		// construct AB.
		M._tmp.resize (A.rowdim (), A.coldim ());

		DenseMatrix<bool>::RowIterator i;
		StandardBasisStream<GF2, DenseMatrix<bool>::Row> stream (F, A.rowdim ());

		for (i = M._tmp.rowBegin (); i != M._tmp.rowEnd (); ++i)
			stream >> *i;

		if (type == UpperTriangular) {
			mzd_trsm_upper_left (A._rep, M._tmp._rep, STRASSEN_MUL_CUTOFF);
			mzd_trsm_upper_left (M._tmp._rep, B._rep, STRASSEN_MUL_CUTOFF);
		}
		else if (type == LowerTriangular) {
			mzd_trsm_lower_left (A._rep, M._tmp._rep, STRASSEN_MUL_CUTOFF);
			mzd_trsm_lower_left (M._tmp._rep, B._rep, STRASSEN_MUL_CUTOFF);
		}
	} else
		_scal (F, M, false, B);

	return B;
}

DenseMatrix<bool> &trsm_impl (const GF2 &F, M4RIModule &M, bool a, const DenseMatrix<bool> &A, DenseMatrix<bool> &B, TriangularMatrixType type, bool diagIsOne,
			      MatrixCategories::RowMatrixTag, MatrixCategories::RowMatrixTag)
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.rowdim () == A.coldim ());

	if (A.rowdim () == 0 || A.coldim () == 0 || B.coldim () == 0)
		return B;

	if (a) {
		if (type == UpperTriangular)
			mzd_trsm_upper_left (A._rep, B._rep, STRASSEN_MUL_CUTOFF);
		else if (type == LowerTriangular)
			mzd_trsm_lower_left (A._rep, B._rep, STRASSEN_MUL_CUTOFF);
		// FIXME: Should throw error at this point -- invalid type
	} else
		_scal (F, M, false, B);

	return B;
}

template <class Iterator>
mzp_t *make_m4ri_permutation (Iterator P_start, Iterator P_end, size_t len)
{
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

template <class Iterator>
DenseMatrix<bool> &permute_rows_impl (const GF2 &F, M4RIModule &M, Iterator P_begin, Iterator P_end, DenseMatrix<bool> &A, MatrixCategories::RowMatrixTag)
{
	mzp_t *P = make_m4ri_permutation (P_begin, P_end, A.rowdim ());
	mzd_apply_p_left (A._rep, P);
	mzp_free (P);
	return A;
}

template <class Iterator>
DenseMatrix<bool> &permute_cols_impl (const GF2 &F, M4RIModule &M, Iterator P_begin, Iterator P_end, DenseMatrix<bool> &A, MatrixCategories::RowMatrixTag)
{
	mzp_t *P = make_m4ri_permutation (P_begin, P_end, A.coldim ());
	mzd_apply_p_right (A._rep, P);
	mzp_free (P);
	return A;
}

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_M4RI_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
