/* lela/blas/level1-gf2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for GF2
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_GF2_H
#define __BLAS_LEVEL1_GF2_H

#include <algorithm>
#include <iostream>

#include "lela/blas/context.h"
#include "lela/ring/gf2.h"
#include "lela/vector/traits.h"
#include "lela/blas/level1-ll.h"

namespace LELA
{

namespace BLAS1 
{

template <>
class _dot<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01);

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01);

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01)
		{ return op (F, M, res, y, x); }

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01)
		{ return op (F, M, res, y, x); }

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01)
		{ return op (F, M, res, y, x); }

	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &dot_impl (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y,
				    VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01);

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const GF2 &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<GF2, Vector1>::RepresentationType (),
				   typename VectorTraits<GF2, Vector2>::RepresentationType ()); }
};

template <>
class _copy<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01)
	{
		lela_check (x.size () == y.size ());

		std::copy (x.word_begin (), x.word_end (), y.word_begin ());

		if (x.size () > 0)
			y.back_word () = x.back_word ();

		return y;
	}

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01)
		{ y.assign (x.begin (), x.end ()); return y; }

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &copy_impl (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01)
		{ y.assign (x.begin (), x.end ()); return y; }

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const GF2 &F, Modules &M, const Vector1 &x, Vector2 &y)
		{ return copy_impl (F, M, x, y,
				    typename VectorTraits<GF2, Vector1>::RepresentationType (),
				    typename VectorTraits<GF2, Vector2>::RepresentationType ()); }
};

template <>
class _axpy<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Generic);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static Vector2 &axpy_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y,
				   VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01);

public:
	template <class Modules, class Vector1, class Vector2>
	static Vector2 &op (const GF2 &F, Modules &M, bool a, const Vector1 &x, Vector2 &y)
		{ return axpy_impl (F, M, a, x, y,
				    typename VectorTraits<GF2, Vector1>::RepresentationType (),
				    typename VectorTraits<GF2, Vector2>::RepresentationType ()); }
};

template <>
class _scal<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector>
	static Vector &scal_impl (const GF2 &F, Modules &M, bool a, Vector &x, VectorRepresentationTypes::Dense01)
		{ if (!a) { std::fill (x.word_begin (), x.word_end (), 0); x.back_word () = 0; } return x; }

	template <class Modules, class Vector>
	static Vector &scal_impl (const GF2 &F, Modules &M, bool a, Vector &x, VectorRepresentationTypes::Sparse01)
		{ if (!a) x.clear (); return x; }

	template <class Modules, class Vector>
	static Vector &scal_impl (const GF2 &F, Modules &M, bool a, Vector &x, VectorRepresentationTypes::Hybrid01)
		{ if (!a) x.clear (); return x; }

public:
	template <class Modules, class Vector>
	static Vector &op (const GF2 &F, Modules &M, bool a, Vector &x)
		{ return scal_impl (F, M, a, x,
				    typename VectorTraits<GF2, Vector>::RepresentationType ()); }
};

template <>
class _permute<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Iterator, class Vector>
	static Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Dense01);

	template <class Modules, class Iterator, class Vector>
	static Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Iterator, class Vector>
	static Vector &permute_impl (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v, VectorRepresentationTypes::Hybrid01);

public:
	template <class Modules, class Iterator, class Vector>
	static Vector &op (const GF2 &F, Modules &M, Iterator P_begin, Iterator P_end, Vector &v)
		{ return permute_impl (F, M, P_begin, P_end, v,
				       typename VectorTraits<GF2, Vector>::RepresentationType ()); }
};

template <>
class _equal<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Dense01)
		{ return op (F, M, y, x); }

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense01, VectorRepresentationTypes::Hybrid01);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Dense01)
		{ return op (F, M, y, x); }

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse01, VectorRepresentationTypes::Hybrid01)
		{ return op (F, M, y, x); }

	template <class Modules, class Vector1, class Vector2>
	static bool equal_impl (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Hybrid01, VectorRepresentationTypes::Hybrid01);

public:
	template <class Modules, class Vector1, class Vector2>
	static bool op (const GF2 &F, Modules &M, const Vector1 &x, const Vector2 &y)
		{ return equal_impl (F, M, x, y,
				     typename VectorTraits<GF2, Vector1>::RepresentationType (),
				     typename VectorTraits<GF2, Vector2>::RepresentationType ()); }
};

template <>
class _is_zero<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector>
	static bool is_zero_impl (const GF2 &F, Modules &M, const Vector &x, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector>
	static bool is_zero_impl (const GF2 &F, Modules &M, const Vector &x, VectorRepresentationTypes::Sparse01)
		{ return x.empty (); }

	template <class Modules, class Vector>
	static bool is_zero_impl (const GF2 &F, Modules &M, const Vector &x, VectorRepresentationTypes::Hybrid01)
		{ return x.empty (); }

public:
	template <class Modules, class Vector>
	static bool op (const GF2 &F, Modules &M, const Vector &x)
		{ return is_zero_impl (F, M, x,
				       typename VectorTraits<GF2, Vector>::RepresentationType ()); }
};

template <>
class _head<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class reference, class Vector>
	static int head_impl (const GF2 &F, Modules &M, reference &a, const Vector &x, VectorRepresentationTypes::Dense01);

	template <class Modules, class reference, class Vector>
	static int head_impl (const GF2 &F, Modules &M, reference &a, const Vector &x, VectorRepresentationTypes::Sparse01)
		{ if (x.empty ()) return -1; a = true; return x.front (); }

	template <class Modules, class reference, class Vector>
	static int head_impl (const GF2 &F, Modules &M, reference &a, const Vector &x, VectorRepresentationTypes::Hybrid01);

public:
	template <class Modules, class reference, class Vector>
	static int op (const GF2 &F, Modules &M, reference &a, const Vector &x)
		{ return head_impl (F, M, a, x,
				    typename VectorTraits<GF2, Vector>::RepresentationType ()); }
};

template <>
class _read<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector>
	static std::istream &read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &v, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector>
	static std::istream &read_impl (const GF2 &F, Modules &M, std::istream &is, const Vector &v, VectorRepresentationTypes::Sparse01);

public:
	template <class Modules, class Vector>
	static std::istream &op (const GF2 &F, Modules &M, std::istream &is, Vector &v)
		{ return read_impl (F, M, is, v,
				    typename VectorTraits<GF2, Vector>::RepresentationType ()); }
};

template <>
class _write<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector>
	static std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Dense01);

	template <class Modules, class Vector>
	static std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Sparse01);

	template <class Modules, class Vector>
	static std::ostream &write_impl (const GF2 &F, Modules &M, std::ostream &os, const Vector &v, VectorRepresentationTypes::Hybrid01);

public:
	template <class Modules, class Vector>
	static std::ostream &op (const GF2 &F, Modules &M, std::ostream &os, const Vector &v)
		{ return write_impl (F, M, os, v,
				     typename VectorTraits<GF2, Vector>::RepresentationType ()); }
};

} // namespace BLAS1

} // namespace LELA

#endif // __BLAS_LEVEL1_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
