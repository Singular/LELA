/* lela/blas/level1-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for Z/p
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_MODULAR_H
#define __BLAS_LEVEL1_MODULAR_H

#include "lela/ring/modular.h"
#include "lela/blas/context.h"
#include "lela/vector/traits.h"
#include "lela/blas/level1-ll.h"

namespace LELA
{

namespace BLAS1 
{

template <class Element>
class _dot<Modular<Element>, typename ZpModule<Element>::Tag>
{
	template <class Vector1, class Vector2>
	static Element &dot_impl (const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
				  VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static Element &dot_impl (const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
				  VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static Element &dot_impl (const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
				  VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Vector1, class Vector2>
	static Element &dot_impl (const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
				  VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const Modular<Element> &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Modular<Element>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<Element>, Vector2>::RepresentationType ()); }

	template <class Modules, class Iterator, class Accessor, class Vector1, class Vector2>
	static Element &op (const Modular<Element> &F, Modules &M, Property<Iterator, Accessor> res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res.ref (), x, y,
				   typename VectorTraits<Modular<Element>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<Element>, Vector2>::RepresentationType ()); }
};

template <>
class _dot<Modular<uint32>, ZpModule<uint32>::Tag>
{
	template <class Vector1, class Vector2>
	static uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Vector1, class Vector2>
	static uint32 &dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const Modular<uint32> &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Modular<uint32>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<uint32>, Vector2>::RepresentationType ()); }
};

} // namespace BLAS1

} // namespace LELA

#endif // __BLAS_LEVEL1_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
