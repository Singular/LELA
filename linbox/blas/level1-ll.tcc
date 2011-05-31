/* linbox/blas/level1-ll.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-1 interface, to be used by implementation-functions
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_LL_TCC
#define __BLAS_LEVEL1_LL_TCC

#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/blas/level1-ll.h"

namespace LinBox
{

/** This namespace contains the level 1 BLAS interface */
namespace BLAS1 
{

template <class reference, class Field, class Modules, class Vector1, class Vector2>
reference &_dot (const Field &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
		 size_t start_idx, size_t end_idx)
	{ return dot_impl (F, M, res, x, y, start_idx, end_idx,
			   typename VectorTraits<Field, Vector1>::VectorCategory (),
			   typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Modules, class Vector>
void _swap (const Field &F, Modules &M, Vector &x, Vector &y)
	{ swap_impl (F, M, x, y,
		     typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &_copy (const Field &F, Modules &M, const Vector1 &x, Vector2 &y)
	{ return copy_impl (F, M, x, y,
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Modules, class Vector1, class Vector2>
Vector2 &_axpy (const Field &F, Modules &M, const typename Field::Element &a, const Vector1 &x, Vector2 &y)
	{ return axpy_impl (F, M, a, x, y,
			    typename VectorTraits<Field, Vector1>::VectorCategory (),
			    typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Modules, class Vector>
Vector &_scal (const Field &F, Modules &M, const typename Field::Element &a, Vector &x)
	{ return scal_impl (F, M, a, x,
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Iterator, class Vector>
Vector &_permute (const Field &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v)
	{ return permute_impl (F, M, P_begin, P_end, v,
			       typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Vector1, class Vector2>
bool _equal (const Field &F, Modules &M, const Vector1 &x, const Vector2 &y)
	{ return equal_impl (F, M, x, y,
			     typename VectorTraits<Field, Vector1>::VectorCategory (),
			     typename VectorTraits<Field, Vector2>::VectorCategory ()); }

template <class Field, class Modules, class Vector>
bool _is_zero (const Field &F, Modules &M, const Vector &x)
	{ return is_zero_impl (F, M, x,
			       typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Vector>
int _head (const Field &F, Modules &M, typename Field::Element &a, const Vector &x)
	{ return head_impl (F, M, a, x,
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Vector>
std::istream &_read (const Field &F, Modules &M, std::istream &is, Vector &v)
	{ return read_impl (F, M, is, v,
			    typename VectorTraits<Field, Vector>::VectorCategory ()); }

template <class Field, class Modules, class Vector>
std::ostream &_write (const Field &F, Modules &M, std::ostream &os, const Vector &v)
	{ return write_impl (F, M, os, v,
			     typename VectorTraits<Field, Vector>::VectorCategory ()); }

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_LL_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
