/* linbox/blas/level2-gf2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for GF2
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_GF2_H
#define __BLAS_LEVEL2_GF2_H

#include <algorithm>
#include <iostream>

#include "linbox/blas/context.h"
#include "linbox/ring/gf2.h"
#include "linbox/vector/traits.h"
#include "linbox/matrix/traits.h"
#include "linbox/blas/level2-ll.h"

namespace LinBox
{

namespace BLAS2
{

template <>
class _gemv<GF2, GenericModule::Tag>
{
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowMatrixTag,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Dense01);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowMatrixTag,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Sparse01);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::RowMatrixTag,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Hybrid01);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::ColMatrixTag,
				   VectorRepresentationTypes::Dense01,
				   VectorRepresentationTypes::Generic);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::ColMatrixTag,
				   VectorRepresentationTypes::Sparse01,
				   VectorRepresentationTypes::Generic);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   size_t start_idx, size_t end_idx,
				   MatrixCategories::ColMatrixTag,
				   VectorRepresentationTypes::Hybrid01,
				   VectorRepresentationTypes::Generic);
public:
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &op (const GF2     &F,
			    Modules       &M,
			    bool           a,
			    const Matrix  &A,
			    const Vector1 &x,
			    bool           b,
			    Vector2       &y,
			    size_t         start_idx = 0,
			    size_t         end_idx = (size_t) -1)
		{ return gemv_impl (F, M, a, A, x, b, y, start_idx, end_idx,
				    typename MatrixIteratorTypes<typename MatrixTraits<Matrix>::MatrixCategory>::MatrixCategory (),
				    typename VectorTraits<GF2, Vector1>::RepresentationType (),
				    typename VectorTraits<GF2, Vector2>::RepresentationType ()); }
};

template <>
class _ger<GF2, GenericModule::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense01,
				 VectorRepresentationTypes::Dense01,
				 MatrixCategories::RowMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Sparse01,
				 VectorRepresentationTypes::Sparse01,
				 MatrixCategories::RowMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Hybrid01,
				 VectorRepresentationTypes::Hybrid01,
				 MatrixCategories::RowMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense01,
				 VectorRepresentationTypes::Dense01,
				 MatrixCategories::ColMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Sparse01,
				 VectorRepresentationTypes::Sparse01,
				 MatrixCategories::ColMatrixTag);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Hybrid01,
				 VectorRepresentationTypes::Hybrid01,
				 MatrixCategories::ColMatrixTag);

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A);
};

} // namespace BLAS2

} // namespace LinBox

#endif // __BLAS_LEVEL2_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
