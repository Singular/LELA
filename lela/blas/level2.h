/* lela/blas/level2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * BLAS Level 2 public interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_H
#define __BLAS_LEVEL2_H

#include "lela/blas/context.h"
#include "lela/blas/level2-ll.h"

namespace LELA
{

/** This namespace contains the level 2 BLAS interface. This includes
 * arithmetic involving both vectors and matrices.
 *
 * \ingroup blas
 */
namespace BLAS2
{

/** General matrix-vector multiply, y <- a Ax + by
 *
 * Currently, specialisations in which the the two vectors x and y are
 * of different representations (e.g. x is sparse and y is dense) are
 * not implemented.
 *
 * The optional parameters start_idx and end_idx permit the
 * computation of the product of only part of vector x by the
 * respective part of the matrix A.
 *
 * @param a Ring::Element scalar a
 * @param A Matrix A
 * @param x Vector x
 * @param b Ring::Element scalar b
 * @param y Vector y, to be replaced by result of calculation
 * @returns Reference to y
 */

template <class Ring, class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &gemv (Context<Ring, Modules>       &ctx,
	       const typename Ring::Element &a,
	       const Matrix                 &A,
	       const Vector1                &x,
	       const typename Ring::Element &b,
	       Vector2                      &y)
	{ return _gemv<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, a, A, x, b, y); }

/** Triangular matrix-vector multiply, x <- Ax, where A is triangular
 *
 * A must be square.
 *
 * x must have a dense or dense 0-1 representation. This function is
 * not available for sparse or hybrid vectors.
 *
 * @param A Matrix A
 * @param x Vector x
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */

template <class Ring, class Modules, class Matrix, class Vector>
Vector &trmv (Context<Ring, Modules> &ctx, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return _trmv<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, A, x, type, diagIsOne); }

/** Triangular matrix-vector solve, x <- A^-1 x, where A is triangular
 *
 * A must be square. If diagIsOne is false and there are zeros on the
 * diagonal, an exception is thrown.
 *
 * x must have a dense or dense 0-1 representation. This function is
 * not available for sparse or hybrid vectors.
 *
 * @param A Matrix A
 * @param x Vector x
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */

template <class Ring, class Modules, class Matrix, class Vector>
Vector &trsv (Context<Ring, Modules> &ctx, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return _trsv<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, A, x, type, diagIsOne); }

/** General rank-1 update, A <- a x y^T + A
 *
 * Currently, specialisations in which the the two vectors x and y are
 * of different representations (e.g. x is sparse and y is dense) are
 * not implemented.
 *
 * @param a Ring::Element scalar a
 * @param x Vector x
 * @param y Vector y
 * @param A Matrix A, to be replaced by result of calculation
 * @returns Reference to A
 */

template <class Ring, class Modules, class Vector1, class Vector2, class Matrix>
Matrix &ger (Context<Ring, Modules> &ctx, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
	{ return _ger<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, a, x, y, A); }

} // namespace BLAS2

} // namespace LELA

#endif // __BLAS_LEVEL2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
