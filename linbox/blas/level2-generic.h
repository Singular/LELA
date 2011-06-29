/* linbox/blas/level2-generic.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Generic implementations of level 2 BLAS interface
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_GENERIC_H
#define __BLAS_LEVEL2_GENERIC_H

#include <algorithm>

#include "linbox/blas/context.h"
#include "linbox/vector/traits.h"
#include "linbox/matrix/traits.h"
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

namespace BLAS2
{

template <class Ring>
class _gemv<Ring, GenericModule::Tag>
{
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::Row,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Dense);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::Row,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Sparse);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Generic);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Sparse,
				   VectorRepresentationTypes::Generic);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::RowCol,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Dense)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixIteratorTypes::Col (),
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::RowCol,
				   VectorRepresentationTypes::Dense,
				   VectorRepresentationTypes::Sparse)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixIteratorTypes::Col (),
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::RowCol,
				   VectorRepresentationTypes::Sparse,
				   VectorRepresentationTypes::Dense)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixIteratorTypes::Col (),
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixIteratorTypes::RowCol,
				   VectorRepresentationTypes::Sparse,
				   VectorRepresentationTypes::Sparse)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixIteratorTypes::Col (),
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }

public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const Ring                   &F,
			    Modules                       &M,
			    const typename Ring::Element &a,
			    const Matrix                  &A,
			    const Vector1                 &x,
			    const typename Ring::Element &b,
			    Vector2                       &y,
			    size_t                         start_idx = 0,
			    size_t                         end_idx = (size_t) -1)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    typename Matrix::IteratorType (),
				    typename VectorTraits<Ring, Vector1>::RepresentationType (),
				    typename VectorTraits<Ring, Vector2>::RepresentationType ()); }
};

template <class Ring>
class _trmv<Ring, GenericModule::Tag>
{
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne);
};

template <class Ring>
class _trsv<Ring, GenericModule::Tag>
{
public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne);
};

template <class Ring>
class _ger<Ring, GenericModule::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense,
				 VectorRepresentationTypes::Dense,
				 MatrixIteratorTypes::Row);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Sparse,
				 VectorRepresentationTypes::Sparse,
				 MatrixIteratorTypes::Row);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense,
				 VectorRepresentationTypes::Dense,
				 MatrixIteratorTypes::Col);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Sparse,
				 VectorRepresentationTypes::Sparse,
				 MatrixIteratorTypes::Col);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense,
				 VectorRepresentationTypes::Dense,
				 MatrixIteratorTypes::RowCol)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<Ring, Vector1>::RepresentationType (),
				   typename VectorTraits<Ring, Vector2>::RepresentationType (),
				   MatrixIteratorTypes::Row ()); }

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Sparse,
				 VectorRepresentationTypes::Sparse,
				 MatrixIteratorTypes::RowCol)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<Ring, Vector1>::RepresentationType (),
				   typename VectorTraits<Ring, Vector2>::RepresentationType (),
				   MatrixIteratorTypes::Row ()); }

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<Ring, Vector1>::RepresentationType (),
				   typename VectorTraits<Ring, Vector2>::RepresentationType (),
				   typename Matrix::IteratorType ()); }
};

} // namespace BLAS2

} // namespace LinBox

#include "linbox/blas/level2-generic.tcc"

#endif // __BLAS_LEVEL2_GENERIC_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
