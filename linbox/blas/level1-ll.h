/* linbox/blas/level1-ll.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-1 interface, to be used by implementation-functions
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_LL_H
#define __BLAS_LEVEL1_LL_H

#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

/** This namespace contains the level 1 BLAS interface */
namespace BLAS1 
{

template <class reference, class Field, class Modules, class Vector1, class Vector2>
reference &_dot (const Field &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
		 size_t start_idx = 0, size_t end_idx = (size_t) -1);

template <class Field, class Modules, class Vector>
void _swap (const Field &F, Modules &M, Vector &x, Vector &y);

template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &_copy (const Field &F, Modules &M, const Vector1 &x, Vector2 &y);

template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &_axpy (const Field &F, Modules &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y);

template <class Field, class Modules, class Vector>
Vector &_scal (const Field &F, Modules &M, const typename Field::Element &a, Vector &x);

template <class Field, class Modules, class Iterator, class Vector>
Vector &_permute (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v);

template <class Field, class Modules, class Vector1, class Vector2>
bool _equal (const Field &F, Modules &M, const Vector1 &x, const Vector2 &y);

template <class Field, class Modules, class Vector>
bool _is_zero (const Field &F, Modules &M, const Vector &x);

template <class Field, class Modules, class Vector>
int _head (const Field &F, Modules &M, typename Field::Element &a, const Vector &x);

template <class Field, class Modules, class Vector>
std::istream &_read (const Field &F, Modules &M, std::istream &is, Vector &v);

template <class Field, class Modules, class Vector>
std::ostream &_write (const Field &F, Modules &M, std::ostream &os, const Vector &v);

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_LL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
