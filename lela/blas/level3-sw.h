/* lela/blas/level3-sw.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * BLAS-interface to Strassen-Winograd implementation
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL3_SW_H
#define __BLAS_LEVEL3_SW_H

#include "lela/blas/context.h"
#include "lela/blas/level3-ll.h"
#include "lela/algorithms/strassen-winograd.h"

namespace LELA
{

/** This namespace contains the level 3 BLAS interface */
namespace BLAS3
{

template <class Ring, class ParentModule>
class _gemm<Ring, StrassenModuleTag<Ring, ParentModule> >
{
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   VectorRepresentationTypes::Generic)
		{ return _gemm<Ring, typename ParentModule::Tag>::op (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   VectorRepresentationTypes::Dense)
		{ return ((StrassenModule<Ring, ParentModule> &) M).sw.gemm (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   VectorRepresentationTypes::Dense01)
		{ return ((StrassenModule<Ring, ParentModule> &) M).sw.gemm (F, M, a, A, B, b, C); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::Row)
		{ return gemm_impl (F, M, a, A, B, b, C, typename VectorTraits<Ring, typename Matrix3::Row>::RepresentationType ()); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::Col)
		{ return gemm_impl (F, M, a, A, B, b, C, typename VectorTraits<Ring, typename Matrix3::Col>::RepresentationType ()); }

	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &gemm_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C,
				   MatrixIteratorTypes::RowCol)
		{ return gemm_impl (F, M, a, A, B, b, C, typename VectorTraits<Ring, typename Matrix3::Row>::RepresentationType ()); }

public:
	template <class Modules, class Matrix1, class Matrix2, class Matrix3>
	static Matrix3 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Ring::Element &b, Matrix3 &C)
		{ return gemm_impl (F, M, a, A, B, b, C, typename Matrix3::IteratorType ()); }
};

} // namespace BLAS3

} // namespace LELA

#endif // __BLAS_LEVEL3_SW_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
