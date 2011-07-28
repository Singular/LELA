/* Copyright 1994-1997 by Givaro Team
 * Copyright 2000-2002 by LELA Team 
 *
 * Created by M. Samama, T. Gautier
 *
 * Modified Jean-Guillaume.Dumas <Jean-Guillaume.Dumas@imag.fr>
 *          B. David Saunders <saunders@cis.udel.edu>,
 *          Bradford Hovinen <hovinen@gmail.com>
 *          Gilles Villard <Gilles.Villard@ens-lyon.fr>
 *                        JGD Random functions back.                          
 *                        (2002/02/12 16:05:24) 
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_integer_H
#define __LELA_integer_H

#include "lela/lela-config.h"
#include "lela/util/double-word.h"

#include <gmpxx.h>

namespace LELA
{

/** \brief This is a representation of arbitrary integers.  
 *
 * It is a wrapper of GMP integers.  Arithmetic operations are via C++
 * infix operator forms (eg. a*b). It is for ``casual'' uses such as
 * characteristics and cardinalities and when initializing field
 * elements. The integers are also represented as a LELA ring for use
 * in integer matrix computation, see integers.h
 *
 * \ingroup lela
 */
typedef mpz_class integer;

typedef signed __LELA_INT8 int8;
typedef signed __LELA_INT16 int16;

/** \memo This is a representation of 32 bit ints, usually equivalent to `int'.
 *
 * The use of `int32' ensures you are working with 
 * 32 bit signed ints, [-2^31..2^31).  Similarly, int8, int16, and int64 are defined.
 *
 * \ingroup lela
 */
typedef signed __LELA_INT32 int32;

typedef signed __LELA_INT64 int64;

typedef unsigned __LELA_INT8 uint8;
typedef unsigned __LELA_INT16 uint16;

/** This is a representation of 32 bit unsigned ints, usually equivalent to `unsigned int'.
 *
 * The use of `uint32' ensures you are working with 
 * 32 bit unsigned ints, [0..2^32).  Similarly, uint8, uint16, and uint64 are defined.
 *
 * \ingroup lela
 */
typedef unsigned __LELA_INT32 uint32;

typedef unsigned __LELA_INT64 uint64;

#ifdef __LELA_UINT128
typedef __LELA_UINT128 uint128;
#else
typedef DoubleWord<unsigned __LELA_INT64> uint128;
#endif // __LELA_UINT128

#ifdef __LELA_UINT256
typedef __LELA_UINT256 uint256;
#endif // __LELA_UINT256

template <class T>
T abs (const T &a) { return a <= 0 ? a * -1 : a; }

/// Compile-time computation of the gcd of two integers
///
/// \ingroup lela
template <int n, int m> struct const_gcd { static const int val = const_gcd<m, n % m>::val; };
template <int n> struct const_gcd<n, 0> { static const int val = n; };

/// Compile-time computation of the lcm of two integers
///
/// \ingroup lela
template <int n, int m> struct const_lcm { static const int val = n * m / const_gcd<n, m>::val; };

} // namespace LELA

#endif // __LELA_integer_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

