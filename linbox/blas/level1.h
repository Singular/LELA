/* linbox/blas/level1.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * BLAS Level 1 public interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_H
#define __BLAS_LEVEL1_H

#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/blas/level1-ll.h"

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
 * representations, but must be defined over the same ring.
 *
 * @param ctx @ref Context object for calculation
 * @param res Ring-element in which to place output
 * @param x First vector
 * @param y Second vector
 * @param start_idx Starting index in vector
 * @param end_idx Ending index in vector, or (size_t) -1 for the whole vector
 * @returns Reference to res
 */

template <class reference, class Ring, class Modules, class Vector1, class Vector2>
reference &dot (Context<Ring, Modules> &ctx, reference &res, const Vector1 &x, const Vector2 &y,
		size_t start_idx = 0, size_t end_idx = (size_t) -1)
	{ return _dot<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, res, x, y, start_idx, end_idx); }

/** Swap two vectors
 *
 * The two vectors must be of the same type
 *
 * @param ctx @ref Context object for calculation
 * @param x First vector
 * @param y Second vector
 */

template <class Ring, class Modules, class Vector>
void swap (Context<Ring, Modules> &ctx, Vector &x, Vector &y)
	{ _swap<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, x, y); }

/** Copy x into y
 *
 * x and y may be of different types, but must be defined over the same ring
 *
 * @param ctx @ref Context object for calculation
 * @param x Origin vector
 * @param y Destination vector
 * @returns Reference to y
 */

template <class Ring, class Modules, class Vector1, class Vector2>
Vector2 &copy (Context<Ring, Modules> &ctx, const Vector1 &x, Vector2 &y)
	{ return _copy<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, x, y); }

/** y -> ax + y
 *
 * x and y may be of different types, but must be defined over the same ring
 *
 * @param ctx @ref Context object for calculation
 * @param a Ring::Element scalar a
 * @param x Vector x
 * @param y Vector y, to be replaced by result of calculation
 * @returns Reference to y
 */

template <class Ring, class Modules, class Vector1, class Vector2>
Vector2 &axpy (Context<Ring, Modules> &ctx, const typename Ring::Element &a, const Vector1 &x, Vector2 &y)
	{ return _axpy<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, a, x, y); }

/** x -> ax
 *
 * @param ctx @ref Context object for calculation
 * @param a Ring::Element scalar a
 * @param x Vector x, to be replaced by result of calculation
 */

template <class Ring, class Modules, class Vector>
Vector &scal (Context<Ring, Modules> &ctx, const typename Ring::Element &a, Vector &x)
	{ return _scal<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, a, x); }

/** Permute entries of v, v <- Pv, where P is a permutation
 *
 * @param P_begin beginning iterator of permutation
 * @param P_end ending iterator of permutation
 * @param v Vector v, to be replaced by result of computation
 * @returns Reference to v
 */

template <class Ring, class Modules, class Iterator, class Vector>
Vector &permute (Context<Ring, Modules> &ctx, Iterator P_begin, Iterator P_end, Vector &v)
	{ return _permute<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, P_begin, P_end, v); }

//@} Operations on vectors

/// @name Queries on vectors
//@{

/** Test whether x and y are equal
 *
 * x and y may be of different types, but must be defined over the same ring
 *
 * @param ctx @ref Context object for calculation
 * @param x First vector
 * @param y Second vector
 * @returns true if equal, false otherwise
 */

template <class Ring, class Modules, class Vector1, class Vector2>
bool equal (Context<Ring, Modules> &ctx, const Vector1 &x, const Vector2 &y)
	{ return _equal<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, x, y); }

/** Test whether x is the zero vector
 *
 * @param ctx @ref Context object for calculation
 * @param x Vector
 * @returns true if x is zero, false otherwise
 */

template <class Ring, class Modules, class Vector>
bool is_zero (Context<Ring, Modules> &ctx, const Vector &x)
	{ return _is_zero<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, x); }

/** Find the first nonzero element of the vector and return its index, or -1 if the vector is zero
 *
 * @param ctx @ref Context object for calculation
 * @param a Ring::Element into which to store first nonzero element
 * @param x Vector
 * @returns Index of first nonzero element, or -1 if the vector is zero
 */

template <class Ring, class Modules, class Vector>
int head (Context<Ring, Modules> &ctx, typename Ring::Element &a, const Vector &x)
	{ return _head<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, a, x); }

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

template <class Ring, class Modules, class Vector>
std::istream &read (Context<Ring, Modules> &ctx, std::istream &is, Vector &v)
	{ return _read<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, is, v); }

/** Write the given vector to the given stream
 *
 * @param ctx Context-object
 * @param os Output-stream
 * @param v Vector to be written
 * @returns Reference to os
 */

template <class Ring, class Modules, class Vector>
std::ostream &write (Context<Ring, Modules> &ctx, std::ostream &os, const Vector &v)
	{ return _write<Ring, typename Modules::Tag>::op (ctx.F, ctx.M, os, v); }

/** Output the given permutation to the given stream
 *
 * @param os Output-stream
 * @param P_begin Beginning of permutation
 * @param P_end End of permutation
 * @returns Reference to os
 */
template <class Iterator>
std::ostream &write_permutation (std::ostream &os, Iterator P_begin, Iterator P_end);

//@} I/O of vectors

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
