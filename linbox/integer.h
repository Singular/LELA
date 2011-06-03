/* Copyright 1994-1997 by Givaro Team
 * Copyright 2000-2002 by LinBox Team 
 * see the COPYING file for license information.
 * Created by M. Samama, T. Gautier
 *
 * Modified Jean-Guillaume.Dumas <Jean-Guillaume.Dumas@imag.fr>
 *          B. David Saunders <saunders@cis.udel.edu>,
 *          Bradford Hovinen <hovinen@cis.udel.edu>
 *          Gilles Villard <Gilles.Villard@ens-lyon.fr>
 *                        JGD Random functions back.                          
 *                        (2002/02/12 16:05:24) 
 *
 */

#ifndef __LINBOX_integer_H
#define __LINBOX_integer_H

#include "linbox/linbox-config.h"

#include <gmpxx.h>

namespace LinBox
{

/** \brief This is a representation of arbitrary integers.  
 *
 * \ingroup linbox
 *
 * It is a wrapper of GMP integers.  Arithmetic operations are via C++
 * infix operator forms (eg. a*b). It is for ``casual'' uses such as
 * characteristics and cardinalities and when initializing field
 * elements.  The integers are also represented as a LinBox ring for use
 * in integer matrix computation, see pid-integers.h or ntl-ZZ.h.
 */
typedef mpz_class integer;

typedef signed __LINBOX_INT8 int8;
typedef signed __LINBOX_INT16 int16;

/** \memo This is a representation of 32 bit ints, usually equivalent to `int'.
 *
 * The use of `int32' ensures you are working with 
 * 32 bit signed ints, [-2^31..2^31).  Similarly, int8, int16, and int64 are defined.
 */
typedef signed __LINBOX_INT32 int32;

typedef signed __LINBOX_INT64 int64;

typedef unsigned __LINBOX_INT8 uint8;
typedef unsigned __LINBOX_INT16 uint16;

/** This is a representation of 32 bit unsigned ints, usually equivalent to `unsigned int'.
 *
 * The use of `uint32' ensures you are working with 
 * 32 bit unsigned ints, [0..2^32).  Similarly, uint8, uint16, and uint64 are defined.
 */
typedef unsigned __LINBOX_INT32 uint32;

typedef unsigned __LINBOX_INT64 uint64;

#ifdef __LINBOX_UINT128
typedef __LINBOX_UINT128 uint128;
#endif // __LINBOX_UINT128

#ifdef __LINBOX_UINT256
typedef __LINBOX_UINT256 uint256;
#endif // __LINBOX_UINT256

template <class T>
T abs (const T &a) { return a <= 0 ? a * -1 : a; }

} // namespace LinBox

#endif // __LINBOX_integer_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

