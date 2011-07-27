/* lela/blas/level3-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 3 BLAS interface for Z/p
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_MODULAR_TCC
#define __BLAS_LEVEL3_MODULAR_TCC

// Disable this entire module if we don't have BLAS
#ifdef __LELA_BLAS_AVAILABLE

#include "lela/blas/level3-modular.h"

namespace LELA
{

namespace BLAS3
{

template <class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm<Modular<float>, typename ZpModule<float>::Tag>::gemm_impl
	(const Modular<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
	 MatrixStorageTypes::Dense, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
{
	lela_check (A.coldim () == B.rowdim ());
	lela_check (A.rowdim () == C.rowdim ());
	lela_check (B.coldim () == C.coldim ());

	float ainvb;

	if (M.block_size == 1 || !F.div (ainvb, b, a))
		return _gemm<Modular<float>, typename ZpModule<float>::Tag::Parent>::op (F, M, a, A, B, b, C);

	if (F.isZero (a) || A.coldim () == 0)
		return _scal<Modular<float>, typename ZpModule<float>::Tag>::op (F, M, b, C);

	if (C.coldim () == 0)
		return C;

	TypeWrapperRing<float> Rp;

	size_t k, first_block_end = A.coldim () % (M.block_size - 1);

	if (first_block_end == 0)
		first_block_end = M.block_size - 1;

	typename Matrix1::ConstSubmatrixType A_sub_1 (A, 0, 0, A.rowdim (), first_block_end);
	typename Matrix2::ConstSubmatrixType B_sub_1 (B, 0, 0, first_block_end, B.coldim ());

	_gemm<TypeWrapperRing<float>, typename ZpModule<float>::Tag::TWParent>::op (Rp, M.TWM, F.one (), A_sub_1, B_sub_1, ainvb, C);

	typename Matrix3::RawIterator i_C;

	for (i_C = C.rawBegin (); i_C != C.rawEnd (); ++i_C)
		ModularTraits<float>::reduce (*i_C, *i_C, F._modulus);

	for (k = first_block_end; k < A.coldim (); k += M.block_size - 1) {
		typename Matrix1::ConstSubmatrixType A_sub (A, 0, k, A.rowdim (), M.block_size - 1);
		typename Matrix2::ConstSubmatrixType B_sub (B, k, 0, M.block_size - 1, B.coldim ());

		_gemm<TypeWrapperRing<float>, typename ZpModule<float>::Tag::TWParent>::op (Rp, M.TWM, F.one (), A_sub, B_sub, F.one (), C);

		for (i_C = C.rawBegin (); i_C != C.rawEnd (); ++i_C)
			ModularTraits<float>::reduce (*i_C, *i_C, F._modulus);
	}

	return _scal<Modular<float>, typename ZpModule<float>::Tag>::op (F, M, a, C);
}

template <class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm<Modular<double>, typename ZpModule<double>::Tag>::gemm_impl
	(const Modular<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
	 MatrixStorageTypes::Dense, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense)
{
	lela_check (A.coldim () == B.rowdim ());
	lela_check (A.rowdim () == C.rowdim ());
	lela_check (B.coldim () == C.coldim ());

	double ainvb;

	if (M.block_size == 1 || !F.div (ainvb, b, a))
		return _gemm<Modular<double>, typename ZpModule<double>::Tag::Parent>::op (F, M, a, A, B, b, C);

	if (F.isZero (a) || A.coldim () == 0)
		return _scal<Modular<double>, typename ZpModule<double>::Tag>::op (F, M, b, C);

	if (C.coldim () == 0)
		return C;

	TypeWrapperRing<double> Rp;

	size_t k, first_block_end = A.coldim () % (M.block_size - 1);

	if (first_block_end == 0)
		first_block_end = M.block_size - 1;

	typename Matrix1::ConstSubmatrixType A_sub_1 (A, 0, 0, A.rowdim (), first_block_end);
	typename Matrix2::ConstSubmatrixType B_sub_1 (B, 0, 0, first_block_end, B.coldim ());

	_gemm<TypeWrapperRing<double>, typename ZpModule<double>::Tag::TWParent>::op (Rp, M.TWM, F.one (), A_sub_1, B_sub_1, ainvb, C);

	typename Matrix3::RawIterator i_C;

	for (i_C = C.rawBegin (); i_C != C.rawEnd (); ++i_C)
		ModularTraits<double>::reduce (*i_C, *i_C, F._modulus);

	for (k = first_block_end; k < A.coldim (); k += M.block_size - 1) {
		typename Matrix1::ConstSubmatrixType A_sub (A, 0, k, A.rowdim (), M.block_size - 1);
		typename Matrix2::ConstSubmatrixType B_sub (B, k, 0, M.block_size - 1, B.coldim ());

		_gemm<TypeWrapperRing<double>, typename ZpModule<double>::Tag::TWParent>::op (Rp, M.TWM, F.one (), A_sub, B_sub, F.one (), C);

		for (i_C = C.rawBegin (); i_C != C.rawEnd (); ++i_C)
			ModularTraits<double>::reduce (*i_C, *i_C, F._modulus);
	}

	return _scal<Modular<double>, typename ZpModule<double>::Tag>::op (F, M, a, C);
}

} // namespace BLAS3

} // namespace LELA

#endif // __LELA_BLAS_AVAILABLE

#endif // __BLAS_LEVEL3_MODULAR_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
