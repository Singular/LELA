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
#include "linbox/matrix/matrix-traits.h"

namespace LinBox
{

/** This namespace contains the level 3 BLAS interface */
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
Matrix &scal (Context<Field, Modules> &ctx, const typename Field::Element &a, Matrix &A)
	{ return _scal (ctx.F, ctx.M, a, A); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix>
Matrix &_scal (const Field &F, Modules &M, const typename Field::Element &a, Matrix &A)
	{ return scal_impl (F, M, a, A, typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

/** Matrix-scalar axpy, B <- a A + B
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element scalar
 * @param A Matrix, to be replaced by result of computation
 * @returns Reference to A
 */
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &axpy (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B)
	{ return _axpy (ctx.F, ctx.M, a, A, B); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_axpy (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B)
	{ return axpy_impl (F, M, a, A, B,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

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
Matrix2 &trmm (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return _trmm (ctx.F, ctx.M, a, A, B, type, diagIsOne); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trmm (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return trmm_impl (F, M, a, A, B, type, diagIsOne,
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
Matrix2 &trsm (Context<Field, Modules> &ctx, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return _trsm (ctx.F, ctx.M, a, A, B, type, diagIsOne); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Matrix1, class Matrix2>
Matrix2 &_trsm (const Field &F, Modules &M, const typename Field::Element &a, const Matrix1 &A, Matrix2 &B, TriangularMatrixType type, bool diagIsOne)
	{ return trsm_impl (F, M, a, A, B, type, diagIsOne,
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix1>::MatrixCategory>::MatrixCategory (),
			    typename MatrixIteratorTypes<typename MatrixTraits<Matrix2>::MatrixCategory>::MatrixCategory ()); }

/** Permute rows of A, A <- PA, where P is a permutation
 *
 * @param P_begin beginning iterator of permutation
 * @param P_end ending iterator of permutation
 * @param A Matrix A, to be replaced by result of computation
 * @returns Reference to A
 */
template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_rows (Context<Field, Modules> &ctx, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return _permute_rows (ctx.F, ctx.M, P_begin, P_end, A); }

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &_permute_rows (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return permute_rows_impl (F, M, P_begin, P_end, A,
				    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

/** Permute columns of A, A <- AP, where P is a permutation
 *
 * @param P_begin beginning iterator of permutation
 * @param P_end ending iterator of permutation
 * @param A Matrix A, to be replaced by result of computation
 * @returns Reference to A
 */
template <class Field, class Modules, class Iterator, class Matrix>
Matrix &permute_cols (Context<Field, Modules> &ctx, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return _permute_cols (ctx.F, ctx.M, P_begin, P_end, A); }

template <class Field, class Modules, class Iterator, class Matrix>
Matrix &_permute_cols (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Matrix &A)
	{ return permute_cols_impl (F, M, P_begin, P_end, A,
				     typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

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
	{ return is_zero_impl (F, M, A, typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }

} // namespace BLAS3

} // namespace LinBox

#include "linbox/blas/level3-generic.h"

#endif // __BLAS_LEVEL3_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
