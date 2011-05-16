/* linbox/blas/level1.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_H
#define __BLAS_LEVEL1_H

#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

/** This namespace contains the level 1 BLAS interface */
namespace BLAS1 
{

/// @name Operations on vectors
//@{

/** Compute the dot-product of the two vectors
 *
 * The optional parameters start_idx and end_idx allow the computation
 * of the dot-product of just part of the vectors.
 *
 * The vectors x and y may be of different types which different
 * representations, but must be defined over the same field.
 *
 * @param ctx @ref Context object for calculation
 * @param res Field-element in which to place output
 * @param x First vector
 * @param y Second vector
 * @param start_idx Starting index in vector
 * @param end_idx Ending index in vector, or (size_t) -1 for the whole vector
 * @returns Reference to res
 */
template <class Field, class Modules, class Vector1, class Vector2>
typename Field::Element &dot (Context<Field, Modules> &ctx, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
			      size_t start_idx = 0, size_t end_idx = (size_t) -1)
	{ return _dot (ctx.F, ctx.M, res, x, y, start_idx, end_idx); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector1, class Vector2>
typename Field::Element &_dot (const Field &F, Modules &M, typename Field::Element &res, const Vector1 &x, const Vector2 &y,
			      size_t start_idx = 0, size_t end_idx = (size_t) -1)
	{ return dot_impl (F, M, res, x, y, start_idx, end_idx,
			   typename VectorTraits<Field, Vector1>::VectorCategory (),
			   typename VectorTraits<Field, Vector2>::VectorCategory ()); }

/** Swap two vectors
 *
 * The two vectors must be of the same type
 *
 * @param ctx @ref Context object for calculation
 * @param x First vector
 * @param y Second vector
 */
template <class Field, class Modules, class Vector>
void swap (Context<Field, Modules> &ctx, Vector &x, Vector &y)
	{ _swap (ctx.F, ctx.M, x, y); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector>
void _swap (const Field &F, Modules &M, Vector &x, Vector &y)
	{ swap_impl (F, M, x, y, typename VectorTraits<Field, Vector>::VectorCategory ()); }

/** Copy x into y
 *
 * x and y may be of different types, but must be defined over the same field
 *
 * @param ctx @ref Context object for calculation
 * @param x Origin vector
 * @param y Destination vector
 * @returns Reference to y
 */
template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &copy (Context<Field, Modules> &ctx, const Vector1 &x, Vector2 &y)
	{ return _copy (ctx.F, ctx.M, x, y); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &_copy (const Field &F, Modules &M, const Vector1 &x, Vector2 &y)
	{ return copy_impl (F, M, x, y,
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

/** y -> ax + y
 *
 * x and y may be of different types, but must be defined over the same field
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element scalar a
 * @param x Vector x
 * @param y Vector y, to be replaced by result of calculation
 * @returns Reference to y
 */
template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &axpy (Context<Field, Modules> &ctx, const typename Field::Element &a, const Vector1 &x, Vector2 &y)
	{ return _axpy (ctx.F, ctx.M, a, x, y); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &_axpy (const Field &F, Modules &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y)
	{ return axpy_impl (F, M, a, x, y,
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

/** x -> ax
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element scalar a
 * @param x Vector x, to be replaced by result of calculation
 */
template <class Field, class Modules, class Vector>
Vector &scal (Context<Field, Modules> &ctx, const typename Field::Element &a, Vector &x)
	{ return _scal (ctx.F, ctx.M, a, x); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector>
Vector &_scal (const Field &F, Modules &M, const typename Field::Element &a, Vector &x)
	{ return scal_impl (F, M, a, x,
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

/** Permute entries of v, v <- Pv, where P is a permutation
 *
 * @param P_begin beginning iterator of permutation
 * @param P_end ending iterator of permutation
 * @param v Vector v, to be replaced by result of computation
 * @returns Reference to v
 */
template <class Field, class Modules, class Iterator, class Vector>
Vector &permute (Context<Field, Modules> &ctx, Iterator P_begin, Iterator P_end, Vector &v)
	{ return _permute (ctx.F, ctx.M, P_begin, P_end, v); }

template <class Field, class Modules, class Iterator, class Vector>
Vector &_permute (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v)
	{ return permute_impl (F, M, P_begin, P_end, v,
			       typename VectorTraits<Field, Vector>::VectorCategory ()); }

//@} Operations on vectors

/// @name Queries on vectors
//@{

/** Test whether x and y are equal
 *
 * x and y may be of different types, but must be defined over the same field
 *
 * @param ctx @ref Context object for calculation
 * @param x First vector
 * @param y Second vector
 * @returns true if equal, false otherwise
 */
template <class Field, class Modules, class Vector1, class Vector2>
bool equal (Context<Field, Modules> &ctx, const Vector1 &x, const Vector2 &y)
	{ return _equal (ctx.F, ctx.M, x, y); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector1, class Vector2>
bool _equal (const Field &F, Modules &M, const Vector1 &x, const Vector2 &y)
	{ return equal_impl (F, M, x, y,
			     typename VectorTraits<Field, Vector1>::VectorCategory (),
			     typename VectorTraits<Field, Vector2>::VectorCategory ()); }

/** Test whether x is the zero vector
 *
 * @param ctx @ref Context object for calculation
 * @param x Vector
 * @returns true if x is zero, false otherwise
 */
template <class Field, class Modules, class Vector>
bool is_zero (Context<Field, Modules> &ctx, const Vector &x)
	{ return _is_zero (ctx.F, ctx.M, x); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector>
bool _is_zero (const Field &F, Modules &M, const Vector &x)
	{ return is_zero_impl (F, M, x, typename VectorTraits<Field, Vector>::VectorCategory ()); }

/** Find the first nonzero element of the vector and return its index, or -1 if the vector is zero
 *
 * @param ctx @ref Context object for calculation
 * @param a Field::Element into which to store first nonzero element
 * @param x Vector
 * @returns Index of first nonzero element, or -1 if the vector is zero
 */
template <class Field, class Modules, class Vector>
int head (Context<Field, Modules> &ctx, typename Field::Element &a, const Vector &x)
	{ return _head (ctx.F, ctx.M, a, x); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector>
int _head (const Field &F, Modules &M, typename Field::Element &a, const Vector &x)
	{ return head_impl (F, M, a, x, typename VectorTraits<Field, Vector>::VectorCategory ()); }

//@} Queries on vectors

/// @name I/O of vectors
//@{

/** Read the given vector from the given stream
 *
 * @param ctx Context-object
 * @param is Input-stream
 * @param v Vector into which to place result
 * @returns Reference to is
 */
template <class Field, class Modules, class Vector>
std::istream &read (Context<Field, Modules> &ctx, std::istream &is, Vector &v)
	{ return _read (ctx.F, ctx.M, is, v); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector>
std::istream &_read (const Field &F, Modules &M, std::istream &is, Vector &v)
	{ return read_impl (F, M, is, v, typename VectorTraits<Field, Vector>::VectorCategory ()); }

/** Write the given vector to the given stream
 *
 * @param ctx Context-object
 * @param os Output-stream
 * @param v Vector to be written
 * @returns Reference to os
 */
template <class Field, class Modules, class Vector>
std::ostream &write (Context<Field, Modules> &ctx, std::ostream &os, const Vector &v)
	{ return _write (ctx.F, ctx.M, os, v); }

/// Version specifying the field and module directly rather than in a Context object
template <class Field, class Modules, class Vector>
std::ostream &_write (const Field &F, Modules &M, std::ostream &os, const Vector &v)
	{ return write_impl (F, M, os, v, typename VectorTraits<Field, Vector>::VectorCategory ()); }

//@} I/O of vectors

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_H

#include "linbox/blas/level1-generic.h"

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
