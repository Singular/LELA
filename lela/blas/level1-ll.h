/* lela/blas/level1-ll.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Low-level BLAS-1 interface, to be used by implementation-functions
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_LL_H
#define __BLAS_LEVEL1_LL_H

#include <iostream>

#include "lela/blas/context.h"

namespace LELA
{

/** This namespace contains the level 1 BLAS interface */
namespace BLAS1 
{

template <class Ring, class ModulesTag>
class _dot
{
public:
	template <class Modules, class T, class Vector1, class Vector2>
	static T &op (const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y)
		{ return _dot<Ring, typename ModulesTag::Parent>::op (F, M, res, x, y); }

	template <class Modules, class Iterator, class Endianness, class Vector1, class Vector2>
	static BitVectorReference<Iterator, Endianness> &op (const Ring &F, Modules &M, BitVectorReference<Iterator, Endianness> &res,
							     const Vector1 &x, const Vector2 &y)
		{ return _dot<Ring, typename ModulesTag::Parent>::op (F, M, res, x, y); }
};

template <class Ring, class ModulesTag>
class _swap
{
public:
	template <class Modules, class Vector>
	static void op (const Ring &F, Modules &M, Vector &x, Vector &y)
		{ return _swap<Ring, typename ModulesTag::Parent>::op (F, M, x, y); }
};

template <class Ring, class ModulesTag>
class _copy
{
public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return _copy<Ring, typename ModulesTag::Parent>::op (F, M, x, y); }
};

template <class Ring, class ModulesTag>
class _axpy
{
public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y)
		{ return _axpy<Ring, typename ModulesTag::Parent>::op (F, M, a, x, y); }
};

template <class Ring, class ModulesTag>
class _scal
{
public:
	template <class Modules, class Vector>
	static Vector &op (const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x)
		{ return _scal<Ring, typename ModulesTag::Parent>::op (F, M, a, x); }
};

template <class Ring, class ModulesTag>
class _permute
{
public:
	template <class Modules, class Iterator, class Vector>
	static Vector &op (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v)
		{ return _permute<Ring, typename ModulesTag::Parent>::op (F, M, P_begin, P_end, v); }
};

template <class Ring, class ModulesTag>
class _equal
{
public:
	template <class Modules, class Vector1, class Vector2>
	static bool op (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y)
		{ return _equal<Ring, typename ModulesTag::Parent>::op (F, M, x, y); }
};

template <class Ring, class ModulesTag>
class _is_zero
{
public:
	template <class Modules, class Vector>
	static bool op (const Ring &F, Modules &M, const Vector &x)
		{ return _is_zero<Ring, typename ModulesTag::Parent>::op (F, M, x); }
};

template <class Ring, class ModulesTag>
class _head
{
public:
	template <class Modules, class Vector>
	static int op (const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x)
		{ return _head<Ring, typename ModulesTag::Parent>::op (F, M, a, x); }
};

template <class Ring, class ModulesTag>
class _read
{
public:
	template <class Modules, class Vector>
	static std::istream &op (const Ring &F, Modules &M, std::istream &is, Vector &v)
		{ return _read<Ring, typename ModulesTag::Parent>::op (F, M, is, v); }
};

template <class Ring, class ModulesTag>
class _write
{
public:
	template <class Modules, class Vector>
	static std::ostream &op (const Ring &F, Modules &M, std::ostream &os, const Vector &v)
		{ return _write<Ring, typename ModulesTag::Parent>::op (F, M, os, v); }
};

} // namespace BLAS1

} // namespace LELA

#include "lela/blas/level1-generic.h"

#endif // __BLAS_LEVEL1_LL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
