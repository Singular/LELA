/* lela/blas/level1-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 1 BLAS interface
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_GENERIC_H
#define __BLAS_LEVEL1_GENERIC_H

#include <algorithm>
#include <iostream>

#include "lela/blas/context.h"
#include "lela/vector/traits.h"
#include "lela/blas/level1-ll.h"

// These are needed for the specialisations of std::swap to that the correct overload is used
#include "lela/vector/bit-subvector-word-aligned.h"
#include "lela/vector/sparse.h"
#include "lela/vector/sparse-subvector.h"

namespace LELA
{

namespace BLAS1 
{

template <class Ring>
class _dot<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class T, class Vector1, class Vector2>
	static T &dot_impl (const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
			    VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Modules, class T, class Vector1, class Vector2>
	static T &dot_impl (const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
			    VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Modules, class T, class Vector1, class Vector2>
	static T &dot_impl (const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
			    VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Modules, class T, class Vector1, class Vector2>
	static T &dot_impl (const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y,
			    VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class T, class Vector1, class Vector2>
	static T &op (const Ring &F, Modules &M, T &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Ring, Vector1>::RepresentationType (),
				   typename VectorTraits<Ring, Vector2>::RepresentationType ()); }
};

template <class Ring>
class _swap<Ring, typename GenericModule<Ring>::Tag>
{
public:
	template <class Modules, class Vector>
	static void op (const Ring &F, Modules &M, Vector &x, Vector &y)
		{ std::swap (x, y); }
};

template <class Ring>
class _copy<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const Ring &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return copy_impl (F, M, x, y,
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }
};

template <class Ring>
class _axpy<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, Vector2 &y)
		{ return axpy_impl (F, M, a, x, y,
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }
};

template <class Ring>
class _scal<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector>
	static Vector &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector>
	static Vector &scal_impl (const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector>
	static Vector &op (const Ring &F, Modules &M, const typename Ring::Element &a, Vector &x)
		{ return scal_impl (F, M, a, x,
				    typename VectorTraits<Ring, Vector>::RepresentationType ()); }
};

template <class Ring>
class _permute<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Iterator, class Vector>
	static Vector &permute_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Dense);

	template <class Modules, class Iterator, class Vector>
	static Vector &permute_impl (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Iterator, class Vector>
	static Vector &op (const Ring &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v)
		{ return permute_impl (F, M, P_begin, P_end, v,
				       typename VectorTraits<Ring, Vector>::RepresentationType ()); }
};

template <class Ring>
class _equal<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return _equal<Ring, typename Modules::Tag>::op (F, M, y, x); }

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector1, class Vector2>
	static bool op (const Ring &F, Modules &M, const Vector1 &x, const Vector2 &y)
		{ return equal_impl (F, M, x, y,
				     typename VectorTraits<Ring, Vector1>::RepresentationType (),
				     typename VectorTraits<Ring, Vector2>::RepresentationType ()); }
};

template <class Ring>
class _is_zero<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector>
	static bool is_zero_impl (const Ring &F, Modules &M, const Vector &x, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector>
	static bool is_zero_impl (const Ring &F, Modules &M, const Vector &x, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector>
	static bool op (const Ring &F, Modules &M, const Vector &x)
		{ return is_zero_impl (F, M, x,
				       typename VectorTraits<Ring, Vector>::RepresentationType ()); }
};

template <class Ring>
class _head<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector>
	static int head_impl (const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector>
	static int head_impl (const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector>
	static int op (const Ring &F, Modules &M, typename Ring::Element &a, const Vector &x)
		{ return head_impl (F, M, a, x,
				    typename VectorTraits<Ring, Vector>::RepresentationType ()); }
};

template <class Ring>
class _read<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector>
	static std::istream &read_impl (const Ring &F, Modules &M, std::istream &is, Vector &v, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector>
	static std::istream &read_impl (const Ring &F, Modules &M, std::istream &is, Vector &v, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector>
	static std::istream &op (const Ring &F, Modules &M, std::istream &is, Vector &v)
		{ return read_impl (F, M, is, v,
				    typename VectorTraits<Ring, Vector>::RepresentationType ()); }
};

template <class Ring>
class _write<Ring, typename GenericModule<Ring>::Tag>
{
	template <class Modules, class Vector>
	static std::ostream &write_impl (const Ring &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Dense);

	template <class Modules, class Vector>
	static std::ostream &write_impl (const Ring &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class Vector>
	static std::ostream &op (const Ring &F, Modules &M, std::ostream &os, const Vector &v)
		{ return write_impl (F, M, os, v,
				     typename VectorTraits<Ring, Vector>::RepresentationType ()); }
};

} // namespace BLAS1

} // namespace LELA

#include "lela/blas/level1-generic.tcc"

#endif // __BLAS_LEVEL1_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
