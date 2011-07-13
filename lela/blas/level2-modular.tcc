/* lela/blas/level2-modular.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 2 BLAS interface for Z/p
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL2_MODULAR_TCC
#define __BLAS_LEVEL2_MODULAR_TCC

#include <algorithm>

#include "lela/blas/level2-modular.h"

namespace LELA
{

namespace BLAS2
{

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint8>, ZpModule<uint8>::Tag>::gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
								      uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
								      VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<uint8> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint8> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::ConstColumn::const_iterator k;
	std::vector<ModularTraits<uint8>::DoubleFatElement>::iterator l, l_end;

	M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + std::min (x.size (), M.block_size);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += std::min (x.size () - (j_end - x.begin ()), M.block_size);

		for (l = M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint8>, ZpModule<uint8>::Tag>::gemv_col_dense (const Modular<uint8> &F, ZpModule<uint8> &M,
								      uint8 a, const Matrix &A, const Vector1 &x, uint8 b, Vector2 &y,
								      VectorRepresentationTypes::Sparse)
{
	lela_check (VectorUtils::hasDim<Modular<uint8> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint8> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::ConstColumn::const_iterator k;
	std::vector<ModularTraits<uint8>::DoubleFatElement>::iterator l, l_end;

	M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + std::min (x.size (), M.block_size);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				M._tmp[k->first] += k->second * *j;

		j_end += std::min (x.size () - (j_end - x.begin ()), M.block_size);

		for (l =M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint16>, ZpModule<uint16>::Tag>::gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
									uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
									VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<uint16> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint16> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = x.begin (), j_end;
	typename Matrix::ConstColumn::const_iterator k;
	std::vector<ModularTraits<uint16>::DoubleFatElement>::iterator l, l_end;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + std::min (x.size (), M.block_size);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += std::min (x.size () - (j_end - x.begin ()), M.block_size);

		for (l = M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint16>, ZpModule<uint16>::Tag>::gemv_col_dense (const Modular<uint16> &F, ZpModule<uint16> &M,
									uint16 a, const Matrix &A, const Vector1 &x, uint16 b, Vector2 &y,
									VectorRepresentationTypes::Sparse)
{
	lela_check (VectorUtils::hasDim<Modular<uint16> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint16> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::ConstColumn::const_iterator k;
	std::vector<ModularTraits<uint16>::DoubleFatElement>::iterator l, l_end;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

	l_end = M._tmp.begin () + y.size ();

	do {
		j = x.begin ();
		j_end = j + std::min (x.size (), M.block_size);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = M._tmp.begin (); k != i->end (); ++k, ++l)
				M._tmp[k->first] += k->second * *j;

		j_end += std::min (x.size () - (j_end - x.begin ()), M.block_size);

		for (l =M._tmp.begin (); l != l_end; ++l)
			*l %= F._modulus;

	} while (j_end != x.end ());

	typename Vector2::iterator y_j;

	for (y_j = y.begin (), l = M._tmp.begin (); y_j != y.end (); ++y_j, ++l)
		*y_j = (a * *l + b * *y_j) % F._modulus;

	return y;
}

template <class Matrix, class Vector1, class Vector2>
Vector2 &_gemv<Modular<uint32>, ZpModule<uint32>::Tag>::gemv_col_dense (const Modular<uint32> &F, ZpModule<uint32> &M,
									uint32 a, const Matrix &A, const Vector1 &x, uint32 b, Vector2 &y,
									VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<uint32> > (x, A.coldim ()));
	lela_check (VectorUtils::hasDim<Modular<uint32> > (y, A.rowdim ()));

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::ConstColumn::const_iterator k;
	std::vector<ModularTraits<uint32>::DoubleFatElement>::iterator l;

	uint64 t;

	if (M._tmp.size () < y.size ())
		M._tmp.resize (y.size ());

	std::fill (M._tmp.begin (), M._tmp.begin () + y.size (), 0);

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
									VectorRepresentationTypes::Sparse)
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
