/* linbox/blas/level2-ll.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_LL_TCC
#define __BLAS_LEVEL2_LL_TCC

#include "linbox/blas/context.h"
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

/** This namespace contains the level 2 BLAS interface */
namespace BLAS2
{

template <class Field, class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv (const Field                   &F,
		Modules                       &M,
		const typename Field::Element &a,
		const Matrix                  &A,
		const Vector1                 &x,
		const typename Field::Element &b,
		Vector2                       &y,
		size_t                         start_idx,
		size_t                         end_idx)
	{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Modules, class Matrix, class Vector>
Vector &_trmv (const Field &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return trmv_impl (F, M, A, x, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Matrix, class Vector>
Vector &_trsv (const Field &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return trsv_impl (F, M, A, x, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger (const Field &F, Modules &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
	{ return ger_impl (F, M, a, x, y, A,
			   typename VectorTraits<Field, Vector1>::VectorCategory (),
			   typename VectorTraits<Field, Vector2>::VectorCategory (),
			   typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_LL_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
