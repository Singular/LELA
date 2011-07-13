/* lela/blas/level1-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for Z/p
 * ------------------------------------
 *
 * See COPYING for license information.
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

template <>
class _dot<Modular<uint8>, ZpModule<uint8>::Tag>
{
	template <class Vector1, class Vector2>
	static uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Vector1, class Vector2>
	static uint8 &dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const Modular<uint8> &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Modular<uint8>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<uint8>, Vector2>::RepresentationType ()); }
};

template <>
class _dot<Modular<uint16>, ZpModule<uint16>::Tag>
{
	template <class Vector1, class Vector2>
	static uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Vector1, class Vector2>
	static uint16 &dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse);

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const Modular<uint16> &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Modular<uint16>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<uint16>, Vector2>::RepresentationType ()); }
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

template <>
class _dot<Modular<float>, ZpModule<float>::Tag>
{
	template <class Vector1, class Vector2>
	static float &dot_impl (const Modular<float> &F, ZpModule<float> &M, float &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static float &dot_impl (const Modular<float> &F, ZpModule<float> &M, float &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static float &dot_impl (const Modular<float> &F, ZpModule<float> &M, float &res, const Vector1 &x, const Vector2 &y,
				VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Vector1, class Vector2>
	static double &dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
		{ return _dot<Modular<double>, ZpModule<float>::Tag::Parent>::op (F, M, res, x, y); }

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const Modular<float> &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Modular<float>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<float>, Vector2>::RepresentationType ()); }
};

template <>
class _dot<Modular<double>, ZpModule<double>::Tag>
{
	template <class Vector1, class Vector2>
	static double &dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense);
	
	template <class Vector1, class Vector2>
	static double &dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense);

	template <class Vector1, class Vector2>
	static double &dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Sparse)
		{ return op (F, M, res, y, x); }

	template <class Vector1, class Vector2>
	static double &dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
				 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
		{ return _dot<Modular<double>, ZpModule<double>::Tag::Parent>::op (F, M, res, x, y); }

public:
	template <class Modules, class reference, class Vector1, class Vector2>
	static reference &op (const Modular<double> &F, Modules &M, reference &res, const Vector1 &x, const Vector2 &y)
		{ return dot_impl (F, M, res, x, y,
				   typename VectorTraits<Modular<double>, Vector1>::RepresentationType (),
				   typename VectorTraits<Modular<double>, Vector2>::RepresentationType ()); }
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
