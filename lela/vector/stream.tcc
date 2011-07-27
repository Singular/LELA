/* lela/vector/stream.tcc
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Stream of vectors mimicing a C++-istream
 *
 * ------------------------------------
 *
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_VECTOR_STREAM_TCC
#define __LELA_VECTOR_STREAM_TCC

#include "lela/blas/level1.h"

namespace LELA
{

template <class Ring, class Vector, class RandIter>
Vector &RandomDenseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Dense>::get (Vector &v) 
{
	typename Vector::iterator i;

	if ( (_m > 0) && (_j++ >= _m) )
		return v;

	for (i = v.begin (); i != v.end (); i++)
		_r.random (*i);

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomSparseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Dense>::get (Vector &v)
{
	double val;

	if (_m > 0 && _j++ >= _m)
		return v;

	for (typename Vector::iterator i = v.begin (); i != v.end (); ++i) {
		val = _MT.randomDouble ();

		if (val < _p)
			_r.random (*i);
		else
			_F.copy (*i, _F.zero ());
	}

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomSparseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Sparse>::get (Vector &v) 
{
	size_t i = (size_t) -1;
	double val;
	int skip;

	if (_m > 0 && _j++ >= _m)
		return v;

	v.clear ();

	while (1) {
		val = _MT.randomDouble ();
		skip = (int) (ceil (log (val) * _1_log_1mp));

		if (skip <= 0)
			++i;
		else
			i += skip;

		if (i >= _n) break;

		v.push_back (std::pair<size_t, typename Ring::Element> (i, typename Ring::Element ()));
		_r.random (v.back ().second);
	}

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomDenseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Dense01>::get (Vector &v) 
{
	typename Vector::word_iterator i;

	if ( (_m > 0) && (_j++ >= _m) )
		return v;

	for (i = v.word_begin (); i != v.word_end (); i++)
		*i = _MT.randomLongLong ();

	typename Vector::word_type t = _MT.randomLongLong ();

	if ((v.size () & WordTraits<typename Vector::word_type>::pos_mask) != 0)
		t &= Vector::Endianness::mask_left (v.size () & WordTraits<typename Vector::word_type>::pos_mask);

	v.back_word () = t;

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomSparseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Sparse01>::get (Vector &v) 
{
	size_t i = (size_t) -1;
	double val;
	int skip;

	if (_m > 0 && _j++ >= _m)
		return v;

	v.clear ();

	while (1) {
		val = _MT.randomDouble ();
		skip = (int) (ceil (log (val) * _1_log_1mp));

		if (skip <= 0)
			i++;
		else
			i += skip;

		if (i >= _n) break;

		v.push_back (i);
	}

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomSparseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Hybrid01>::get (Vector &v) 
{
	size_t i = (size_t) -1;
	double val;
	int skip;

	if (_m > 0 && _j++ >= _m)
		return v;

	v.clear ();

	while (1) {
		val = _MT.randomDouble ();
		skip = (int) (ceil (log (val) * _1_log_1mp));

		if (skip <= 0)
			i++;
		else
			i += skip;

		if (i >= _n) break;

		typename Vector::index_type idx = i >> WordTraits<typename Vector::word_type>::logof_size;
		typename Vector::word_type mask = Vector::Endianness::e_j (i & WordTraits<typename Vector::word_type>::pos_mask);

		if (!v.empty () && idx == v.back ().first)
			v.back ().second |= mask;
		else
			v.push_back (typename Vector::value_type (idx, mask));
	}

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomHybridStream<Ring, Vector, RandIter>::get (Vector &v) 
{
	size_t i = (size_t) -1;
	double val;
	int skip;

	typename Vector::word_type t;

	if (_m > 0 && _j++ >= _m)
		return v;

	v.clear ();

	while (1) {
		val = _MT.randomDouble ();
		skip = (int) (ceil (log (val) * _1_log_1mp));

		if (skip <= 0)
			i++;
		else
			i += skip;

		if (i << WordTraits<typename Vector::word_type>::logof_size >= _n) break;

		t = _MT.randomLongLong ();

		if (i == _n >> WordTraits<typename Vector::word_type>::logof_size)
			t &= Vector::Endianness::mask_left (_n & WordTraits<typename Vector::word_type>::pos_mask);

		v.push_back (typename Vector::value_type (i, t));
	}

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Dense>::get (Vector &v) 
{
	BLAS1::scal (_ctx, _ctx.F.zero (), v);
	_ctx.F.copy (v[_j++], _ctx.F.one ());

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Sparse>::get (Vector &v) 
{
	v.clear ();

	if (_j < _n) {
		v.push_back (typename Vector::value_type (_j++, typename Ring::Element ()));
		_F.copy (v.back ().second, _F.one ());
	}

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Dense01>::get (Vector &v) 
{
	std::fill (v.word_begin (), v.word_end (), 0);
	v.back_word () = 0;

	v[_j] = true;

	_j++;

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Sparse01>::get (Vector &v) 
{
	v.clear ();

	if (_j < _n)
		v.push_back (_j++);

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Hybrid01>::get (Vector &v) 
{
	typedef WordTraits<typename Vector::word_type> WT;

	v.clear ();
	v.push_back (typename Vector::value_type (_j >> WT::logof_size, Vector::Endianness::e_j (_j & WT::pos_mask)));
	++_j;
	return v;
}

} // namespace LELA

#endif // __LELA_VECTOR_STREAM_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
