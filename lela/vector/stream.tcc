/* lela/vector/stream.tcc
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *  
 * See COPYING for license information.
 */

#ifndef __LELA_vector_stream_TCC
#define __LELA_vector_stream_TCC

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
			_F.assign (*i, _zero);
	}

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomSparseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Sparse>::get (Vector &v) 
{
	typename Ring::Element x;
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

		_r.random (x);
		v.push_back (std::pair<size_t, typename Ring::Element> (i, x));
	}

	return v;
}

template <class Ring, class Vector, class RandIter>
Vector &RandomDenseStream<Ring, Vector, RandIter, VectorRepresentationTypes::Dense01>::get (Vector &v) 
{
	typename Vector::iterator i;

	if ( (_m > 0) && (_j++ >= _m) )
		return v;

	for (i = v.begin (); i != v.end (); i++)
		_r.random (*i);

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

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Dense>::get (Vector &v) 
{
	static typename Ring::Element zero;
	typename Vector::iterator i;
	size_t idx;

	for (i = v.begin (), idx = 0; i != v.end (); i++, idx++) {
		if (idx == _j)
			_F.init (*i, 1);
		else
			_F.assign (*i, zero);
	}

	_j++;

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Sparse>::get (Vector &v) 
{
	v.clear ();

	if (_j < _n)
		v.push_back (std::pair <size_t, typename Ring::Element> (_j++, _one));

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Dense01 >::get (Vector &v) 
{
	std::fill (v.word_begin (), v.word_end (), 0);
	v.back_word () = 0;

	v[_j] = true;

	_j++;

	return v;
}

template <class Ring, class Vector>
Vector &StandardBasisStream<Ring, Vector, VectorRepresentationTypes::Hybrid01>::get (Vector &v) 
{
	typedef WordTraits<typename Vector::word_type> WT;

	v.clear ();
	v.push_back (_j >> WT::logof_size, Vector::Endianness::e_j (_j & WT::pos_mask));
	return v;
}

} // namespace LELA

#endif // __LELA_vector_stream_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
