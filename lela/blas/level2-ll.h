/* lela/blas/level2-ll.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-2 interface, to be used by implementation-functions
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_LL_H
#define __BLAS_LEVEL2_LL_H

#include "lela/blas/context.h"

namespace LELA
{

/** This namespace contains the level 2 BLAS interface */
namespace BLAS2
{

template <class Ring, class ModulesTag>
class _gemv
{
public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const Ring                   &F,
			    Modules                      &M,
			    const typename Ring::Element &a,
			    const Matrix                 &A,
			    const Vector1                &x,
			    const typename Ring::Element &b,
			    Vector2                      &y)
		{ return _gemv<Ring, typename ModulesTag::Parent>::op (F, M, a, A, x, b, y); }
};

template <class Ring, class ModulesTag>
class _trmv
{
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return _trmv<Ring, typename ModulesTag::Parent>::op (F, M, A, x, type, diagIsOne); }
};

template <class Ring, class ModulesTag>
class _trsv
{
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return _trsv<Ring, typename ModulesTag::Parent>::op (F, M, A, x, type, diagIsOne); }
};

template <class Ring, class ModulesTag>
class _ger
{
public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return _ger<Ring, typename ModulesTag::Parent>::op (F, M, a, x, y, A); }
};

} // namespace BLAS2

} // namespace LELA

#include "lela/blas/level2-generic.h"

#endif // __BLAS_LEVEL2_LL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
