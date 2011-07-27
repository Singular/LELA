/* lela/blas/level2-modular.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 * Portions derived from LinBox, linbox/field/modular-int32.h
 *
 * Implementations of level 2 BLAS interface for Z/p
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL2_MODULAR_TCC
#define __BLAS_LEVEL2_MODULAR_TCC

#include <algorithm>

#include "lela/blas/level2-modular.h"
#include "lela/ring/type-wrapper.h"

namespace LELA
{

namespace BLAS2
{

template <class Element>
template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<Element>, typename ZpModule<Element>::Tag>::gemv_impl
	(const Modular<Element> &F, ZpModule<Element> &M,
	 Element a, const Matrix &A, const Vector1 &x, Element b, Vector2 &y,
	 MatrixIteratorTypes::Col,
	 VectorRepresentationTypes::Dense,
	 VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<Element> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<Element> > (y, A.rowdim ()));

	if (M.block_size == 1)
		return _gemv<Modular<Element>, typename ZpModule<Element>::Tag::Parent>::op (F, M, a, A, x, b, y);

	size_t block_size = (M.block_size == 0) ? x.size () : M.block_size;
	size_t first_col = x.size () % block_size;

	typename Vector1::const_iterator i;
	typename Vector1::const_iterator block_1_end_x = x.begin () + first_col;
	typename std::vector<typename ModularTraits<Element>::DoubleFatElement>::iterator l;

	M._tmp.resize (y.size ());
	std::fill (M._tmp.begin (), M._tmp.end (), 0);

	TypeWrapperRing<Element> Rp;

	Subvector<typename Vector1::const_iterator> x_sub_1 (x.begin (), block_1_end_x);
	typename Matrix::ConstSubmatrixType A_sub_1 (A, 0, 0, A.rowdim (), first_col);

	_gemv<TypeWrapperRing<Element>, typename ZpModule<Element>::Tag::TWParent>::op (Rp, M.TWM, Rp.one (), A_sub_1, x_sub_1, Rp.one (), M._tmp);

	for (l = M._tmp.begin (); l != M._tmp.end (); ++l)
		ModularTraits<Element>::reduce (*l, *l, F._modulus);

	for (i = block_1_end_x; i != x.end (); i += block_size, first_col += block_size) {
		Subvector<typename Vector1::const_iterator> x_sub (i, i + block_size);
		typename Matrix::ConstSubmatrixType A_sub (A, 0, first_col, A.rowdim (), block_size);

		_gemv<TypeWrapperRing<Element>, typename ZpModule<Element>::Tag::TWParent>::op (Rp, M.TWM, Rp.one (), A_sub, x_sub, Rp.one (), M._tmp);

		for (l = M._tmp.begin (); l != M._tmp.end (); ++l)
			ModularTraits<Element>::reduce (*l, *l, F._modulus);
	}

	typename Vector2::iterator j;

	for (j = y.begin (), l = M._tmp.begin (); j != y.end (); ++j, ++l)
		ModularTraits<Element>::reduce (*j, a * *l + b * *j, F._modulus);

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint32>, ZpModule<uint32>::Tag>::gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
									uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
									VectorRepresentationTypes::Dense,
									VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<uint32> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint32> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector1::const_iterator j;
	typename Matrix::ConstColumn::const_iterator k;
	std::vector<ModularTraits<uint32>::DoubleFatElement>::iterator l;

	uint64 t;

	M._tmp.resize (y.size ());
	std::fill (M._tmp.begin (), M._tmp.end (), 0);

	for (j = x.begin (); j != x.end (); ++j, ++i) {
		for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) *k) * ((uint64) *j);

			*l += t;

			if (*l < t)
				*l += M.two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l) {
		uint32 al, byj;

		F.mul (al, a, uint32 (*l % F._modulus));
		F.mul (byj, b, *y_j);
		F.add (*y_j, al, byj);
	}

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint32>, ZpModule<uint32>::Tag>::gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
									uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
									VectorRepresentationTypes::Sparse,
									VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<uint32> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint32> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<ModularTraits<uint32>::DoubleFatElement>::iterator l;

	uint64 t;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	for (j = x.begin (); j != x.end (); ++j, ++i) {
		for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) k->second) * ((uint64) *j);

			M._tmp[k->first] += t;

			if (M._tmp[k->first] < t)
				M._tmp[k->first] += M.two_64;
		}
	}

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l) {
		uint32 al, byj;

		F.mul (al, a, uint32 (*l % F._modulus));
		F.mul (byj, b, *y_j);
		F.add (*y_j, al, byj);
	}

	return y;
}

} // namespace BLAS2

} // namespace LELA

#endif // __BLAS_LEVEL2_MODULAR_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
