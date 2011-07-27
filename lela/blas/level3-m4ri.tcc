/* lela/blas/level3-m4ri.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 3 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_M4RI_TCC
#define __BLAS_LEVEL3_M4RI_TCC

#include "lela/blas/level3-m4ri.h"
#include "lela/blas/level3-ll.h"

// The macro STRASSEN_MUL_CUTOFF was renamed in M4RI 20110601
#ifndef __LELA_HAVE_M4RI_GE_20110601
#  define __M4RI_STRASSEN_MUL_CUTOFF STRASSEN_MUL_CUTOFF
#endif // !__LELA_HAVE_M4RI_GE_20110601

namespace LELA
{

namespace BLAS3
{

template <class Modules>
M4RIMatrixBase &_scal<GF2, M4RIModule::Tag>::op (const GF2 &F, Modules &M, bool a, M4RIMatrixBase &A)
{
	size_t i;

	if (A.rowdim () == 0 || A.coldim () == 0)
		return A;

	if (!a)
		for (i = 0; i < A.rowdim (); ++i)
			mzd_row_clear_offset (A._rep, i, 0);
	return A;
}

template <class Modules>
M4RIMatrixBase &_axpy<GF2, M4RIModule::Tag>::op (const GF2 &F, Modules &M, bool a, const M4RIMatrixBase &A, M4RIMatrixBase &B)
{
	lela_check (A.rowdim () == B.rowdim ());
	lela_check (A.coldim () == B.coldim ());

	if (A.rowdim () == 0 || A.coldim () == 0)
		return B;

	if (a)
		mzd_add (B._rep, B._rep, A._rep);

	return B;
}

template <class Modules>
M4RIMatrixBase &_gemm<GF2, M4RIModule::Tag>::op (const GF2 &F, Modules &M,
						 bool a, const M4RIMatrixBase &A, const M4RIMatrix &B, bool b, M4RIMatrixBase &C)
{
	lela_check (A.rowdim () == C.rowdim ());
	lela_check (A.coldim () == B.rowdim ());
	lela_check (B.coldim () == C.coldim ());

	if (A.rowdim () == 0 || B.rowdim () == 0 || B.coldim () == 0)
		return C;

	if (a) {
		if (b)
			mzd_addmul (C._rep, A._rep, B._rep, __M4RI_STRASSEN_MUL_CUTOFF);
		else
			mzd_mul (C._rep, A._rep, B._rep, __M4RI_STRASSEN_MUL_CUTOFF);
	}
	else if (!b)
		_scal<GF2, typename Modules::Tag>::op (F, M, false, C);

	return C;
}

template <class Modules>
M4RIMatrixBase &_trsm<GF2, M4RIModule::Tag>::op (const GF2 &F, Modules &M, bool a, const M4RIMatrixBase &A, M4RIMatrixBase &B, TriangularMatrixType type, bool diagIsOne)
{
	lela_check (A.rowdim () == B.rowdim ());
	lela_check (A.rowdim () == A.coldim ());

	if (A.rowdim () == 0 || A.coldim () == 0 || B.coldim () == 0)
		return B;

	if (a) {
		if (type == UpperTriangular)
			mzd_trsm_upper_left (A._rep, B._rep, __M4RI_STRASSEN_MUL_CUTOFF);
		else if (type == LowerTriangular)
			mzd_trsm_lower_left (A._rep, B._rep, __M4RI_STRASSEN_MUL_CUTOFF);
		// FIXME: Should throw error at this point -- invalid type
	} else
		_scal<GF2, typename Modules::Tag>::op (F, M, false, B);

	return B;
}

template <class Iterator>
mzp_t *make_m4ri_permutation (Iterator P_start, Iterator P_end, size_t len)
{
	mzp_t *P;
	Iterator i_P;

	P = mzp_init (len);

	for (i_P = P_start; i_P != P_end; ++i_P) {
		lela_check (i_P->first < len);
		lela_check (i_P->second < len);
		P->values[i_P->first] = P->values[i_P->second];
	}

	return P;
}

template <class Modules, class Iterator>
M4RIMatrixBase &_permute_rows<GF2, M4RIModule::Tag>::op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, M4RIMatrixBase &A)
{
	mzp_t *P = make_m4ri_permutation (P_begin, P_end, A.rowdim ());
	mzd_apply_p_left (A._rep, P);
	mzp_free (P);
	return A;
}

template <class Modules, class Iterator>
M4RIMatrixBase &_permute_cols<GF2, M4RIModule::Tag>::op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, M4RIMatrixBase &A)
{
	mzp_t *P = make_m4ri_permutation (P_begin, P_end, A.coldim ());
	mzd_apply_p_right (A._rep, P);
	mzp_free (P);
	return A;
}

} // namespace BLAS3

} // namespace LELA

#endif // __BLAS_LEVEL3_M4RI_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
