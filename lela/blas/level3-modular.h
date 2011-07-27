/* lela/blas/level3-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 3 BLAS interface for Z/p
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_MODULAR_H
#define __BLAS_LEVEL3_MODULAR_H

// Disable this entire module if we don't have BLAS
#ifdef __LELA_BLAS_AVAILABLE

#include "lela/ring/modular.h"
#include "lela/blas/context.h"
#include "lela/matrix/traits.h"
#include "lela/blas/level3-ll.h"

namespace LELA
{

namespace BLAS3
{

template <>
class _gemm<Modular<float>, ZpModule<float>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Modular<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _gemm<Modular<float>, typename ZpModule<float>::Tag::Parent>::op (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Modular<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense);

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const Modular<float> &F, Modules &M, float a, const Matrix1 &A, const Matrix2 &B, float b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C,
				    typename Matrix1::StorageType (),
				    typename Matrix2::StorageType (),
				    typename Matrix3::StorageType ()); }
};

template <>
class _gemm<Modular<double>, ZpModule<double>::Tag>
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Modular<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::Generic, MatrixStorageTypes::Generic, MatrixStorageTypes::Generic)
		{ return _gemm<Modular<double>, typename ZpModule<double>::Tag::Parent>::op (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Modular<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C,
				   MatrixStorageTypes::Dense, MatrixStorageTypes::Dense, MatrixStorageTypes::Dense);

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const Modular<double> &F, Modules &M, double a, const Matrix1 &A, const Matrix2 &B, double b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C,
				    typename Matrix1::StorageType (),
				    typename Matrix2::StorageType (),
				    typename Matrix3::StorageType ()); }
};

} // namespace BLAS3

} // namespace LELA

#endif // __LELA_BLAS_AVAILABLE

#endif // __BLAS_LEVEL3_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
