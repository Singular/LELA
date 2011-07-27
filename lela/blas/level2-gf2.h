/* lela/blas/level2-gf2.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for GF2
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_GF2_H
#define __BLAS_LEVEL2_GF2_H

#include <algorithm>
#include <iostream>

#include "lela/blas/context.h"
#include "lela/ring/gf2.h"
#include "lela/vector/traits.h"
#include "lela/matrix/traits.h"
#include "lela/blas/level2-ll.h"

namespace LELA
{

namespace BLAS2
{

template <>
class _gemv<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   MatrixIteratorTypes::Row,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Dense01);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   MatrixIteratorTypes::Row,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Sparse01);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   MatrixIteratorTypes::Row,
				   VectorRepresentationTypes::Generic,
				   VectorRepresentationTypes::Hybrid01);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Dense01,
				   VectorRepresentationTypes::Generic);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   MatrixIteratorTypes::Col,
				   VectorRepresentationTypes::Sparse01,
				   VectorRepresentationTypes::Generic);

	template <class Modules, class Matrix, class Vector1, class Vector2>
	static Vector2 &gemv_impl (const GF2 &F, Modules &M,
				   bool a, const Matrix &A, const Vector1 &x, bool b, Vector2 &y,
				   MatrixIteratorTypes::Col,
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
			    Vector2       &y)
		{ return gemv_impl (F, M, a, A, x, b, y,
				    typename Matrix::IteratorType (),
				    typename VectorTraits<GF2, Vector1>::RepresentationType (),
				    typename VectorTraits<GF2, Vector2>::RepresentationType ()); }
};

template <>
class _ger<GF2, GenericModule<GF2>::Tag>
{
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Dense01,
				 VectorRepresentationTypes::Generic,
				 MatrixIteratorTypes::Row);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Sparse01,
				 VectorRepresentationTypes::Generic,
				 MatrixIteratorTypes::Row);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Hybrid01,
				 VectorRepresentationTypes::Generic,
				 MatrixIteratorTypes::Row);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorRepresentationTypes::Dense01,
				 MatrixIteratorTypes::Col);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorRepresentationTypes::Sparse01,
				 MatrixIteratorTypes::Col);

	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &ger_impl (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A,
				 VectorRepresentationTypes::Generic,
				 VectorRepresentationTypes::Hybrid01,
				 MatrixIteratorTypes::Col);

public:
	template <class Modules, class Vector1, class Vector2, class Matrix>
	static Matrix &op (const GF2 &F, Modules &M, bool a, const Vector1 &x, const Vector2 &y, Matrix &A)
		{ return ger_impl (F, M, a, x, y, A,
				   typename VectorTraits<GF2, Vector1>::RepresentationType (),
				   typename VectorTraits<GF2, Vector2>::RepresentationType (),
				   typename Matrix::IteratorType ()); }
};

} // namespace BLAS2

} // namespace LELA

#endif // __BLAS_LEVEL2_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
