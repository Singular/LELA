/* linbox/blas/level2-ll.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-2 interface, to be used by implementation-functions
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_LL_H
#define __BLAS_LEVEL2_LL_H

#include "linbox/blas/context.h"

namespace LinBox
{

/** This namespace contains the level 2 BLAS interface */
namespace BLAS2
{

template <class Ring, class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv (const Ring                   &F,
		Modules                       &M,
		const typename Ring::Element &a,
		const Matrix                  &A,
		const Vector1                 &x,
		const typename Ring::Element &b,
		Vector2                       &y,
		size_t                         start_idx = 0,
		size_t                         end_idx = (size_t) -1);

template <class Ring, class Modules, class Matrix, class Vector>
Vector &_trmv (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne);

template <class Ring, class Modules, class Matrix, class Vector>
Vector &_trsv (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne);

template <class Ring, class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A);

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_LL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
