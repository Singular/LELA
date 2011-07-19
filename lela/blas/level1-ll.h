/* lela/blas/level1-ll.h
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

#include "lela/blas/context.h"

namespace LELA
{

/** This namespace contains the level 1 BLAS interface */
namespace BLAS1 
{

/** Notes on structure of low-level interface
 *
 * Each BLAS-routine has its own class whose name is the name of the
 * routine preceeded by an underscore. Within the routine is a single
 * static method called op. The classes are parametrised by ring and
 * an empty structure called ModulesTag, which can be retrieved from a
 * Modules-class with Modules::Tag. The method is then parametrised
 * normally by Modules (whose tag need not be ModulesTag, in case
 * ModulesTag refers to a less specialised implementation) and the
 * vector- and matrix-types.
 *
 * To create a specialised implementation or a higher level algorithm,
 * just create the corresponding Modules class with an empty Tag and
 * partially specialise the corresponding class to that tag. If the
 * implementation only works with a given ring, then the class may
 * also be specialised to the corresponding ring.
 *
 * The default implementation of the method op simply retrieves the
 * parent of the tag (accessible via ModulesTag::Parent) and invokes
 * the same method with the same operations on the corresponding class
 * instanted with that tag. This permits that a given module,
 * identified by a certain ModulesTag, need not implement all
 * BLAS-routines -- when it fails to implement a routine, a call to
 * that routine is just passed to the parent.
 *
 * The reason to put these methods into classes, as opposed to using
 * naked functions, is to ensure that the correct specialisation is
 * always invoked. Otherwise the compiler may in some cases (which are
 * notoriously tricky and hard to predict) actually prefer a more
 * generic parametrised method to a more specialised one, perhaps even
 * invoking the most generic method in a never ending recursion.
 *
 * When a method must be further specified by, say, vector-type or
 * matrix-iterator-type, these specialisations should appear as
 * private static methods invoked by the method op.
 */

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
