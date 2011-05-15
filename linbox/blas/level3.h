/* linbox/blas/level3.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL3_H
#define __BLAS_LEVEL3_H

#include "linbox/blas/context.h"

namespace LinBox
{

/** This namespace contains the level-2 BLAS interface */
namespace BLAS3
{

/** Copy A into B
 *
 * A and B must both be defined over the same field.
 *
 * @param ctx @ref Context object for calculation
 * @param A Origin matrix
 * @param B Destination matrix
 * @returns Reference to B
 */
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &copy (Context<Field, Modules> &ctx, const Matrix1 &A, Matrix2 &B)
	{ return _copy (ctx.F, ctx.M, A, B); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_copy (const Field &F, Modules &M, const Matrix1 &A, Matrix2 &B)
	{ return copy_impl (F, M, A, B,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory ()); }

/** Scale the matrix A by a scalar, A <- a A
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element scalar
 * @param A Matrix, to be replaced by result of computation
 * @returns Reference to A
 */
template <class Field, class Modules, class Matrix>
bool scal (Context<Field, Modules> &ctx, const typename Field::Element &a, Matrix &A)
	{ return _scal (ctx.F, ctx.M, a, A); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix>
bool _scal (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A)
	{ return scal_impl (F, M, a, A, typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

/** General matrix-matrix multiply, C <- a AB + b C
 *
 * Not all combinations of iterator-types are available.
 *
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param B Matrix B
 * @param b Field::Element scalar b
 * @param C Matrix C, to be replaced by result of calculation
 * @returns Reference to y
 */
template <class Field, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &gemm (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C)
	{ return _gemm (ctx.F, ctx.M, a, A, B, b, C); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2, class Matrix3>
Matrix3 &_gemm (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, const Matrix2 &B, const typename Field::Element &b, Matrix3 &C)
	{ return gemm_impl (F, M, a, A, B, b, C,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix3>::MatrixCategory>::MatrixCategory ()); }

/** Triangular matrix-matrix multiply, B <- a AB, where A is triangular
 *
 * A must be square
 *
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param B Matrix B
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &trmm (Context<Field, Modules> &ctx, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return _trmm (ctx.F, ctx.M, A, B, type, diagIsOne); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trmm (const Field &F, Modules &M, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return trmm_impl (F, M, A, B, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

/** Triangular matrix-matrix solve, B <- a A^-1 B, where A is triangular
 *
 * A must be square. If diagIsOne is false and there are zeros on the
 * diagonal, an exception is thrown.
 *
 * @param a Field::Element scalar a
 * @param A Matrix A
 * @param B Matrix B
 * @param type Whether A is upper or lower triangular
 * @param diagIsOne Whether to assume that the entries on the diagonal of A are one
 * @returns Reference to y
 */
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &trsm (Context<Field, Modules> &ctx, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return _trsm (ctx.F, ctx.M, A, B, type, diagIsOne); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trsm (const Field &F, Modules &M, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return trsm_impl (F, M, A, B, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

/** Test whether A and B are equal
 *
 * @param ctx @ref Context object for calculation
 * @param A First matrix
 * @param B Second matrix
 * @returns true if equal, false otherwise
 */
template <class Field, class Modules, class Matrix1, class Matrix2>
bool equal (Context<Field, Modules> &ctx, const Matrix1 &A, const Matrix2 &B)
	{ return _equal (ctx.F, ctx.M, A, B); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
bool _equal (const Field &F, Modules &M, const Matrix1 &A, const Matrix2 &B)
	{ return equal_impl (F, M, A, B,
			     typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			     typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

/** Test whether A is the zero matrix
 *
 * @param ctx @ref Context object for calculation
 * @param A Matrix
 * @returns true if A is zero, false otherwise
 */
template <class Field, class Modules, class Matrix>
bool is_zero (Context<Field, Modules> &ctx, const Matrix &A)
	{ return _is_zero (ctx.F, ctx.M, A); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix>
bool _is_zero (const Field &F, Modules &M, const Matrix &A)
	{ return is_zero_impl (F, M, A, typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

} // namespace BLAS3

} // namespace LinBox

#endif // __BLAS_LEVEL3_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
