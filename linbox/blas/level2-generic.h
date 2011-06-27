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
				   MatrixCategories::RowMatrixTag,
				   VectorCategories::GenericVectorTag,
				   VectorCategories::DenseVectorTag);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowMatrixTag,
				   VectorCategories::GenericVectorTag,
				   VectorCategories::SparseVectorTag);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::ColMatrixTag,
				   VectorCategories::DenseVectorTag,
				   VectorCategories::GenericVectorTag);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::ColMatrixTag,
				   VectorCategories::SparseVectorTag,
				   VectorCategories::GenericVectorTag);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowColMatrixTag,
				   VectorCategories::DenseVectorTag,
				   VectorCategories::DenseVectorTag)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixCategories::ColMatrixTag (),
				    typename VectorTraits<Ring, Vector1>::VectorCategory (),
				    typename VectorTraits<Ring, Vector2>::VectorCategory ()); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowColMatrixTag,
				   VectorCategories::DenseVectorTag,
				   VectorCategories::SparseVectorTag)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixCategories::ColMatrixTag (),
				    typename VectorTraits<Ring, Vector1>::VectorCategory (),
				    typename VectorTraits<Ring, Vector2>::VectorCategory ()); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowColMatrixTag,
				   VectorCategories::SparseVectorTag,
				   VectorCategories::DenseVectorTag)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixCategories::ColMatrixTag (),
				    typename VectorTraits<Ring, Vector1>::VectorCategory (),
				    typename VectorTraits<Ring, Vector2>::VectorCategory ()); }

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const Ring &F, Modules &M,
				   const typename Ring::Element &a, const Matrix &A, const Vector1 &x, const typename Ring::Element &b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowColMatrixTag,
				   VectorCategories::SparseVectorTag,
				   VectorCategories::SparseVectorTag)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    MatrixCategories::ColMatrixTag (),
				    typename VectorTraits<Ring, Vector1>::VectorCategory (),
				    typename VectorTraits<Ring, Vector2>::VectorCategory ()); }

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
				    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
				    typename VectorTraits<Ring, Vector1>::VectorCategory (),
				    typename VectorTraits<Ring, Vector2>::VectorCategory ()); }
};

template <class Ring>
class _trmv<Ring, GenericModule::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trmv_impl (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixCategories::RowMatrixTag,
				  VectorCategories::DenseVectorTag);

public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trmv_impl (F, M, A, x, type, diagIsOne,
				    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
				    typename VectorTraits<Ring, Vector>::VectorCategory ()); }
};

template <class Ring>
class _trsv<Ring, GenericModule::Tag>
{
	template <class Modules, class Matrix, class Vector>
	static Vector &trsv_impl (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne,
				  MatrixCategories::RowMatrixTag,
				  VectorCategories::DenseVectorTag);

public:
	template <class Modules, class Matrix, class Vector>
	static Vector &op (const Ring &F, Modules &M, const Matrix &A, Vector &x, TriangularMatrixType type, bool diagIsOne)
		{ return trsv_impl (F, M, A, x, type, diagIsOne,
				    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
				    typename VectorTraits<Ring, Vector>::VectorCategory ()); }
};

template <class Ring>
class _ger<Ring, GenericModule::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorCategories::DenseVectorTag,
				 VectorCategories::DenseVectorTag,
				 MatrixCategories::RowMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorCategories::SparseVectorTag,
				 VectorCategories::SparseVectorTag,
				 MatrixCategories::RowMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorCategories::DenseVectorTag,
				 VectorCategories::DenseVectorTag,
				 MatrixCategories::ColMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorCategories::SparseVectorTag,
				 VectorCategories::SparseVectorTag,
				 MatrixCategories::ColMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorCategories::DenseVectorTag,
				 VectorCategories::DenseVectorTag,
				 MatrixCategories::RowColMatrixTag)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<Ring, Vector1>::VectorCategory (),
				   typename VectorTraits<Ring, Vector2>::VectorCategory (),
				   MatrixCategories::RowMatrixTag ()); }

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorCategories::SparseVectorTag,
				 VectorCategories::SparseVectorTag,
				 MatrixCategories::RowColMatrixTag)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<Ring, Vector1>::VectorCategory (),
				   typename VectorTraits<Ring, Vector2>::VectorCategory (),
				   MatrixCategories::RowMatrixTag ()); }

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const Ring &F, Modules &M, const typename Ring::Element &a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<Ring, Vector1>::VectorCategory (),
				   typename VectorTraits<Ring, Vector2>::VectorCategory (),
				   typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory ()); }
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
