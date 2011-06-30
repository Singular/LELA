/* linbox/blas/level2-modular.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for Z/p
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_MODULAR_H
#define __BLAS_LEVEL2_MODULAR_H

#include "linbox/ring/modular.h"
#include "linbox/blas/context.h"
#include "linbox/vector/traits.h"
#include "linbox/matrix/traits.h"
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

namespace BLAS2
{

template <>
class _gemv<Modular<uint8>, ZpModule<uint8>::Tag>
{
	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
					uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
					VectorRepresentationTypes::Dense);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
					uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
					VectorRepresentationTypes::Sparse);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint8> &F, ZpModule<uint8> &M,
				   uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Generic)
		{ return gemv_col_dense (F, M, a, A, x, b, y, typename VectorTraits<Modular<uint8>, typename Matrix::Column>::RepresentationType ()); }

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint8> &F, ZpModule<uint8> &M,
				   uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
				   MatrixIteratorTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Generic)
		{ return _gemv<Modular<uint8>, ZpModule<uint8>::Tag::Parent>::op (F, M, a, A, x, b, y); }

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const Modular<uint8> &F,
			    Modules              &M,
			    uint8                 a,
			    const Matrix         &A,
			    const Vector1        &x,
			    uint8                 b,
			    Vector2              &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::IteratorType (),
				    typename VectorTraits<Modular<uint8>, Vector1>::RepresentationType (),
				    typename VectorTraits<Modular<uint8>, Vector2>::RepresentationType ()); }
};

template <>
class _gemv<Modular<uint16>, ZpModule<uint16>::Tag>
{
	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
					uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
					VectorRepresentationTypes::Dense);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
					uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
					VectorRepresentationTypes::Sparse);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint16> &F, ZpModule<uint16> &M,
				   uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Generic)
		{ return gemv_col_dense (F, M, a, A, x, b, y, typename VectorTraits<Modular<uint16>, typename Matrix::Column>::RepresentationType ()); }

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint16> &F, ZpModule<uint16> &M,
				   uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
				   MatrixIteratorTypes::Generic,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Generic)
		{ return _gemv<Modular<uint16>, ZpModule<uint16>::Tag::Parent>::op (F, M, a, A, x, b, y); }

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const Modular<uint16> &F,
			    Modules               &M,
			    uint16                 a,
			    const Matrix          &A,
			    const Vector1         &x,
			    uint16                 b,
			    Vector2               &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::IteratorType (),
				    typename VectorTraits<Modular<uint16>, Vector1>::RepresentationType (),
				    typename VectorTraits<Modular<uint16>, Vector2>::RepresentationType ()); }
};

template <>
class _gemv<Modular<uint32>, ZpModule<uint32>::Tag>
{
	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
					uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
					VectorRepresentationTypes::Dense);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
					uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
					VectorRepresentationTypes::Sparse);

	template <class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Modular<uint32> &F, ZpModule<uint32> &M,
				   uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Generic)
		{ return gemv_col_dense (F, M, a, A, x, b, y, typename VectorTraits<Modular<uint32>, typename Matrix::Column>::RepresentationType ()); }

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

} // namespace LinBox

#endif // __BLAS_LEVEL2_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
