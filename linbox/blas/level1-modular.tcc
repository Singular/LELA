/* linbox/blas/level1-modular.tcc
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Implementations of level 1 BLAS interface for Z/p
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLAS_LEVEL1_MODULAR_TCC
#define __BLAS_LEVEL1_MODULAR_TCC

#include "linbox/blas/level1-modular.h"

namespace LinBox
{

namespace BLAS1 
{

template <class Vector1, class Vector2>
uint8 &_dot<Modular<uint8>, ZpModule<uint8>::Tag>::dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
							     VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	linbox_check (x.size () == y.size ());

	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	typename Vector1::const_iterator iterend = x.begin () + (x.size () % M.block_size);

	uint64 t = 0;

	for (; i != iterend; ++i, ++j)
		t += (uint64) *i * (uint64) *j;

	t %= (uint64) F._modulus;

	for (; iterend != x.end (); j += M.block_size) {
		typename Vector1::const_iterator iter_i = iterend;
		typename Vector2::const_iterator iter_j;

		iterend += M.block_size;

		for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
			t += (uint64) *iter_i * (uint64) *j;

		t %= (uint64) F._modulus;
	}

	return res = t;
}

template <class Vector1, class Vector2>
uint8 &_dot<Modular<uint8>, ZpModule<uint8>::Tag>::dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
							     VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	linbox_check (VectorUtils::hasDim<Modular<uint8> > (x, y.size ()));

	typename Vector1::const_iterator i = x.begin ();

	uint64 t = 0;

	if (x.end () - i < (long) M.block_size) {
		for (; i != x.end (); ++i)
			t += (uint64) i->second * (uint64) y[i->first];

		return res = t % (uint64) F._modulus;
	} else {
		// i still points to the beginning
		typename Vector1::const_iterator iterend = i + (x.end () - i) % M.block_size;

		for (; i != iterend; ++i)
			t += (uint64) i->second * (uint64) y[i->first];

		t %= (uint64) F._modulus;

		while (iterend != x.end ()) {
			typename Vector1::const_iterator iter_i = iterend;

			iterend += M.block_size;
			i += M.block_size;

			for (; iter_i != iterend; ++iter_i)
				t += (uint64) iter_i->second * (uint64) y[iter_i->first];

			t %= (uint64) F._modulus;
		}

		return res = t;
	}
}

template <class Vector1, class Vector2>
uint8 &_dot<Modular<uint8>, ZpModule<uint8>::Tag>::dot_impl (const Modular<uint8> &F, ZpModule<uint8> &M, uint8 &res, const Vector1 &x, const Vector2 &y,
							     VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	uint64 t = 0, count;

	while (i != x.end () && j != y.end ()) {
		for (count = 0; count < M.block_size && i != x.end () && j != y.end (); ++i) {
			while (j != y.end () && j->first < i->first) ++j;

			if (j != y.end () && i->first == j->first) {
				t += (uint64) i->second * (uint64) j->second;
				++count;
			}
		}

		t %= (uint64) F._modulus;
	}

	return res = t;
}

template <class Vector1, class Vector2>
uint16 &_dot<Modular<uint16>, ZpModule<uint16>::Tag>::dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	linbox_check (x.size () == y.size ());

	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	typename Vector1::const_iterator iterend = x.begin () + x.size () % M.block_size;

	uint64 t = 0;

	for (; i != iterend; ++i, ++j)
		t += (uint64) *i * (uint64) *j;

	t %= (uint64) F._modulus;

	for (; iterend != x.end (); j += M.block_size) {
		typename Vector1::const_iterator iter_i = iterend;
		typename Vector2::const_iterator iter_j;

		iterend += M.block_size;

		for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
			t += (uint64) *iter_i * (uint64) *j;

		t %= (uint64) F._modulus;
	}

	return res = t;
}

template <class Vector1, class Vector2>
uint16 &_dot<Modular<uint16>, ZpModule<uint16>::Tag>::dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	linbox_check (VectorUtils::hasDim<Modular<uint16> > (x, y.size ()));

	typename Vector1::const_iterator i = x.begin ();

	uint64 t = 0;

	if (x.end () - i < (long) M.block_size) {
		for (; i != x.end (); ++i)
			t += (uint64) i->second * (uint64) y[i->first];

		return res = t % (uint64) F._modulus;
	} else {
		// i still points to the beginning
		typename Vector1::const_iterator iterend = i + (x.end () - i) % M.block_size;

		for (; i != iterend; ++i)
			t += (uint64) i->second * (uint64) y[i->first];

		t %= (uint64) F._modulus;

		while (iterend != x.end ()) {
			typename Vector1::const_iterator iter_i = iterend;

			iterend += M.block_size;
			i += M.block_size;

			for (; iter_i != iterend; ++iter_i)
				t += (uint64) iter_i->second * (uint64) y[iter_i->first];

			t %= (uint64) F._modulus;
		}

		return res = t;
	}
}

template <class Vector1, class Vector2>
uint16 &_dot<Modular<uint16>, ZpModule<uint16>::Tag>::dot_impl (const Modular<uint16> &F, ZpModule<uint16> &M, uint16 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Sparse)
{
	typename Vector1::const_iterator i = x.begin ();
	typename Vector2::const_iterator j = y.begin ();

	uint64 t = 0, count;

	while (i != x.end () && j != y.end ()) {
		for (count = 0; count < M.block_size && i != x.end () && j != y.end (); ++i) {
			while (j != y.end () && j->first < i->first) ++j;

			if (j != y.end () && i->first == j->first) {
				t += (uint64) i->second * (uint64) j->second;
				++count;
			}
		}

		t %= (uint64) F._modulus;
	}

	return res = t;
}

template <class Vector1, class Vector2>
uint32 &_dot<Modular<uint32>, ZpModule<uint32>::Tag>::dot_impl (const Modular<uint32> &F, ZpModule<uint32> &M, uint32 &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	linbox_check (x.size () == y.size ());

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
	linbox_check (VectorUtils::hasDim<Modular<uint32> > (x, y.size ()));

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

template <class Vector1, class Vector2>
float &_dot<Modular<float>, ZpModule<float>::Tag>::dot_impl (const Modular<float> &F, ZpModule<float> &M, float &res, const Vector1 &x, const Vector2 &y,
							     VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	float s = 0.;
	float t = 0.;

	if (x.size () < M.block_size) {
		for (size_t i = 0; i < x.size (); ++i)
			s += x[i] * y[i];

		s = fmod (s, F._modulus);
	} else {
		size_t i = 0;

		for (; i < x.size () - M.block_size ;i = i + M.block_size) {
			for (size_t j = i; j < i + M.block_size; ++j)
				s += x[j] * y[j];

			t += fmod (s, F._modulus);
			s = 0.;
		}

		for (; i < x.size (); ++i)
			t += x[i] * y[i];

		t += fmod (s, F._modulus);
		s = fmod (t, F._modulus);
	}

	return res = s;
}

template <class Vector1, class Vector2>
float &_dot<Modular<float>, ZpModule<float>::Tag>::dot_impl (const Modular<float> &F, ZpModule<float> &M, float &res, const Vector1 &x, const Vector2 &y,
							     VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	float s = 0.;
	float t = 0.;

	if (x.size () < M.block_size) {
		for (size_t i = 0; i < x.size (); ++i)
			s += x[i].second * y[x[i].first];

		s = fmod (s, F._modulus);
	} else {
		size_t i = 0;

		for (; i < x.size() - M.block_size; i = i + M.block_size) {
			for (size_t j = i; j < i + M.block_size; ++j)
				s += x[j].second * y[x[j].first];

			t += fmod (s, F._modulus);
			s = 0.;
		}
		for (; i < x.size (); ++i)
			s += x[i].second * y[x[i].first];

		t += fmod (s, F._modulus);
		s = fmod (t, F._modulus);
	}

	return res = s;
}

template <class Vector1, class Vector2>
double &_dot<Modular<double>, ZpModule<double>::Tag>::dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Dense, VectorRepresentationTypes::Dense)
{
	double s = 0.;
	double t = 0.;

	if (x.size () < M.block_size) {
		for (size_t i = 0; i < x.size(); ++i)
			s += x[i] * y[i];

		s = fmod (s, F._modulus);
	} else {			
		size_t i = 0;

		for (; i < x.size () - M.block_size; i = i + M.block_size) {
			for (size_t j = i; j < i + M.block_size; ++j)
				s += x[j] * y[j];

			t += fmod (s, F._modulus);
			s = 0.;
		}

		for (; i < x.size (); ++i)
			s += x[i] * y[i];

		t += fmod (s, F._modulus);
		s = fmod (t, F._modulus);
	}

	return res = s;
}

template <class Vector1, class Vector2>
double &_dot<Modular<double>, ZpModule<double>::Tag>::dot_impl (const Modular<double> &F, ZpModule<double> &M, double &res, const Vector1 &x, const Vector2 &y,
								VectorRepresentationTypes::Sparse, VectorRepresentationTypes::Dense)
{
	double s = 0.;
	double t = 0.;

	if (x.size () < M.block_size) {
		for (size_t i = 0; i < x.size (); ++i)
			s += x[i].second * y[x[i].first];

		s = fmod (s, F._modulus);
	}
	else {
		size_t i = 0;

		for (; i < x.size () - M.block_size; i = i + M.block_size) {
			for (size_t j = i; j < i + M.block_size; ++j)
				s += x[j].second * y[x[j].first];

			t += fmod (s, F._modulus);
			s = 0.;
		}

		for (; i < x.size (); ++i)
			s += x[i].second * y[x[i].first];

		t += fmod (s, F._modulus);
		s = fmod (t, F._modulus);
	}

	return res = s;
}

} // namespace BLAS1

} // namespace LinBox

#endif // __BLAS_LEVEL1_MODULAR_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
