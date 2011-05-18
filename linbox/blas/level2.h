/* linbox/blas/level2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_H
#define __BLAS_LEVEL2_H

#include "linbox/blas/context.h"
#include "linbox/blas/level2-generic.h"

namespace LinBox
{

/** This namespace contains the level 2 BLAS interface */
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
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param x Vector x
 * @param b Field::Element scalar b
 * @param y Vector y, to be replaced by result of calculation
 * @param start_idx Starting index in vector x of product
 * @param end_idx Ending index in vector x of product (-1 for the whole vector)
 * @returns Reference to y
 */

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &_gemv (const Field                   &F,
		Modules                       &M,
		const typename Field::Element &a,
		const Matrix                  &A,
		const Vector1                 &x,
		const typename Field::Element &b,
		Vector2                       &y,
		size_t                         start_idx = 0,
		size_t                         end_idx = (size_t) -1)
	{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Modules, class Matrix, class Vector1, class Vector2>
Vector2 &gemv (Context<Field, Modules>       &ctx,
	       const typename Field::Element &a,
	       const Matrix                  &A,
	       const Vector1                 &x,
	       const typename Field::Element &b,
	       Vector2                       &y,
	       size_t                         start_idx = 0,
	       size_t                         end_idx = (size_t) -1)
	{ return _gemv (ctx.F, ctx.M, a, A, x, b, y, start_idx, end_idx); }

/** Triangular matrix-vector multiply, x <- Ax, where A is triangular
 *
 * A must be square
 *
 * @param A Matrix A
 * @param x Vector x
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix, class Vector>
Vector &_trmv (const Field &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return trmv_impl (F, M, A, x, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Matrix, class Vector>
Vector &trmv (Context<Field, Modules> &ctx, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return _trmv (ctx.F, ctx.M, A, x, type, diagIsOne); }

/** Triangular matrix-vector solve, x <- A^-1 x, where A is triangular
 *
 * A must be square. If diagIsOne is false and there are zeros on the
 * diagonal, an exception is thrown.
 *
 * @param A Matrix A
 * @param x Vector x
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix, class Vector>
Vector &_trsv (const Field &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return trsv_impl (F, M, A, x, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Matrix, class Vector>
Vector &trsv (Context<Field, Modules> &ctx, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
	{ return _trsv (ctx.F, ctx.M, A, x, type, diagIsOne); }

/** General rank-1 update, A <- a x y^T + A
 *
 * Currently, specialisations in which the the two vectors x and y are
 * of different representations (e.g. x is sparse and y is dense) are
 * not implemented.
 *
 * @param a Field::Element scalar a
 * @param x Vector x
 * @param y Vector y
 * @param A Matrix A, to be replaced by result of calculation
 * @returns Reference to A
 */

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector1, class Vector2, class Matrix>
Matrix &_ger (const Field &F, Modules &M, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
	{ return ger_impl (F, M, a, x, y, A,
			   typename VectorTraits<Field, Vector1>::VectorCategory (),
			   typename VectorTraits<Field, Vector2>::VectorCategory (),
			   typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

template <class Field, class Modules, class Vector1, class Vector2, class Matrix>
Matrix &ger (Context<Field, Modules> &ctx, const typename Field::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
	{ return _ger (ctx.F, ctx.M, a, x, y, A); }

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
