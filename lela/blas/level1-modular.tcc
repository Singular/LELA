/* lela/blas/level1-modular.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 * Portions derived from LinBox, linbox/field/modular-int32.h
 *
 * Implementations of level 1 BLAS interface for Z/p
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __BLAS_LEVEL1_MODULAR_TCC
#define __BLAS_LEVEL1_MODULAR_TCC

#include "lela/blas/level1-modular.h"
#include "lela/ring/type-wrapper.h"

namespace LELA
{

namespace BLAS1 
{

template <class Element>
template <class Vector1, class Vector2>
Element &_dot<Modular<Element>, typename ZpModule<Element>::Tag>::dot_impl
	(const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	lela_check (x.size () == y.size ());

	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;

	size_t block_size = (M.block_size == 0) ? x.size () : M.block_size;

	typename Vector1::const_iterator block_1_end_x = x.begin () + (x.size () % block_size);
	typename Vector2::const_iterator block_1_end_y = y.begin () + (y.size () % block_size);

	typename ModularTraits<Element>::DoubleFatElement t = 0, s;

	TypeWrapperRing<Element> Rp;

	Subvector<typename Vector1::const_iterator> x_sub_1 (x.begin (), block_1_end_x);
	Subvector<typename Vector2::const_iterator> y_sub_1 (y.begin (), block_1_end_y);

	_dot<TypeWrapperRing<Element>, typename ZpModule<Element>::Tag::TWParent>::op (Rp, M.TWM, t, x_sub_1, y_sub_1);

	ModularTraits<Element>::reduce (t, t, F._modulus);

	for (i = block_1_end_x, j = block_1_end_y; i != x.end (); i += block_size, j += block_size) {
		Subvector<typename Vector1::const_iterator> x_sub (i, i + block_size);
		Subvector<typename Vector2::const_iterator> y_sub (j, j + block_size);

		_dot<TypeWrapperRing<Element>, typename ZpModule<Element>::Tag::TWParent>::op (Rp, M.TWM, s, x_sub, y_sub);

		t += s;
		ModularTraits<Element>::reduce (t, t, F._modulus);
	}

	return res = t;
}

template <class Element>
template <class Vector1, class Vector2>
Element &_dot<Modular<Element>, typename ZpModule<Element>::Tag>::dot_impl
	(const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<Element> > (x, y.size ()));

	size_t block_size = (M.block_size == 0) ? x.size () : M.block_size;

	typename Vector1::const_iterator i;

	typename ModularTraits<Element>::DoubleFatElement t = 0;

	if (x.size () <= block_size) {
		for (i = x.begin (); i != x.end (); ++i)
			t += (typename ModularTraits<Element>::DoubleFatElement) i->second * (typename ModularTraits<Element>::DoubleFatElement) y[i->first];

		return ModularTraits<Element>::reduce (res, t, F._modulus);
	} else {
		typename Vector1::const_iterator iterend = x.begin () + x.size () % block_size;

		for (i = x.begin (); i != iterend; ++i)
			t += (typename ModularTraits<Element>::DoubleFatElement) i->second * (typename ModularTraits<Element>::DoubleFatElement) y[i->first];

		ModularTraits<Element>::reduce (t, t, F._modulus);

		while (iterend != x.end ()) {
			iterend += block_size;

			for (; i != iterend; ++i)
				t += (typename ModularTraits<Element>::DoubleFatElement) i->second * (typename ModularTraits<Element>::DoubleFatElement) y[i->first];

			ModularTraits<Element>::reduce (t, t, F._modulus);
		}

		return res = t;
	}
}

template <class Element>
template <class Vector1, class Vector2>
Element &_dot<Modular<Element>, typename ZpModule<Element>::Tag>::dot_impl
	(const Modular<Element> &F, ZpModule<Element> &M, Element &res, const Vector1 &x, const Vector2 &y,
	 VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	size_t block_size = (M.block_size == 0) ? x.size () : M.block_size;

	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	typename ModularTraits<Element>::DoubleFatElement t = 0;
	size_t count;

	while (i != x.end () && j != y.end ()) {
		for (count = 0; count < block_size && i != x.end () && j != y.end (); ++i) {
			while (j != y.end () && j->first < i->first) ++j;

			if (j != y.end () && i->first == j->first) {
				t += (typename ModularTraits<Element>::DoubleFatElement) i->second * (typename ModularTraits<Element>::DoubleFatElement) j->second;
				++count;
			}
		}

		ModularTraits<Element>::reduce (t, t, F._modulus);
	}

	return res = t;
}

template <class Vector1, class Vector2>
uint32 &_dot<Modular<uint32>, ZpModule<uint32>::Tag>::dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	lela_check (x.size () == y.size ());

	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();
  
	uint64 s = 0;
	uint64 t;

	for (; i != x.end (); ++i, ++j) {
		t = (uint64) *i * (uint64) *j;
		s += t;

		if (s < t)
			s += M.two_64;
	}
  
	s %= (uint64) F._modulus;

	return res = s;
}

template <class Vector1, class Vector2>
uint32 &_dot<Modular<uint32>, ZpModule<uint32>::Tag>::dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	lela_check (VectorUtils::hasDim<Modular<uint32> > (x, y.size ()));

	typename Vector1::const_iterator i;

	uint64 s = 0, t;

	for (i = x.begin (); i != x.end (); ++i) {
		t = (uint64) i->second * (uint64) y[i->first];
		s += t;

		if (s < t)
			s += M.two_64;
	}
  
	s %= (uint64) F._modulus;

	return res = s;
}

template <class Vector1, class Vector2>
uint32 &_dot<Modular<uint32>, ZpModule<uint32>::Tag>::dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	uint64 s = 0, t;

	for (; i != x.end () && j != y.end (); ++i) {
		while (j != y.end () && j->first < i->first) ++j;

		if (j != y.end () && i->first == j->first) {
			t = (uint64) i->second * (uint64) j->second;
			s += t;

			if (s < t)
				s += M.two_64;
		}
	}

	return res = s % (uint64) F._modulus;
}

} // namespace BLAS1

} // namespace LELA

#endif // __BLAS_LEVEL1_MODULAR_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
