/* lela/blas/level2-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for Z/p
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_MODULAR_H
#define __BLAS_LEVEL2_MODULAR_H

#include "lela/ring/modular.h"
#include "lela/blas/context.h"
#include "lela/vector/traits.h"
#include "lela/matrix/traits.h"
#include "lela/blas/level2-ll.h"

namespace LELA
{

namespace BLAS2
{

template <class Element>
class _gemv<Modular<Element>, typename ZpModule<Element>::Tag>
{
	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<Element> &F, ZpModule<Element> &M,
				   Element a, const Matrix &A, const Vector1 &x, Element b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Dense);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<Element> &F, ZpModule<Element> &M,
				   Element a, const Matrix &A, const Vector1 &x, Element b, Vector2 &y,
				   MatrixIteratorTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Generic)
		{ return _gemv<Modular<Element>, typename ZpModule<Element>::Tag::Parent>::op (F, M, a, A, x, b, y); }

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const Modular<Element> &F,
			    Modules                &M,
			    Element                 a,
			    const Matrix           &A,
			    const Vector1          &x,
			    Element                 b,
			    Vector2                &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::IteratorType (),
				    typename VectorTraits<Modular<Element>, Vector1>::RepresentationType (),
				    typename VectorTraits<Modular<Element>, Vector2>::RepresentationType ()); }
};

template <>
class _gemv<Modular<uint32>, ZpModule<uint32>::Tag>
{
	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
					uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
					VectorRepresentationTypes::Dense,
					VectorRepresentationTypes::Dense);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
					uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
					VectorRepresentationTypes::Sparse,
					VectorRepresentationTypes::Dense);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
					uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
					VectorRepresentationTypes::Generic,
					VectorRepresentationTypes::Generic)
		{ return _gemv<Modular<uint32>, ZpModule<uint32>::Tag::Parent>::op (F, M, a, A, x, b, y); }

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint32> &F, ZpModule<uint32> &M,
				   uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Generic)
		{ return gemv_col_dense (F, M, a, A, x, b, y,
					 typename VectorTraits<Modular<uint32>, typename Matrix::Column>::RepresentationType (),
					 typename VectorTraits<Modular<uint32>, Vector2>::RepresentationType ()); }

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint32> &F, ZpModule<uint32> &M,
				   uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
				   MatrixIteratorTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Generic)
		{ return _gemv<Modular<uint32>, ZpModule<uint32>::Tag::Parent>::op (F, M, a, A, x, b, y); }

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const Modular<uint32> &F,
			    Modules               &M,
			    uint32                 a,
			    const Matrix          &A,
			    const Vector1         &x,
			    uint32                 b,
			    Vector2               &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::IteratorType (),
				    typename VectorTraits<Modular<uint32>, Vector1>::RepresentationType (),
				    typename VectorTraits<Modular<uint32>, Vector2>::RepresentationType ()); }
};

} // namespace BLAS2

} // namespace LELA

#endif // __BLAS_LEVEL2_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
