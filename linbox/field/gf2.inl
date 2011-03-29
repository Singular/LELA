/* linbox/field/gf2.inl
 * Copyright 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_gf2_INL
#define __LINBOX_field_gf2_INL

#include <iostream>
#include <time.h>

#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/stream.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/matrix/matrix-domain.h"

#include <cctype> //isdigit

namespace LinBox 
{ 

// Specialization of canonical vector types

template <>
class RawVector<bool>
{
    public:
	typedef BitVector<> Dense;
	typedef std::vector<size_t> Sparse;
	typedef std::vector<size_t> SparseSeq;
	typedef std::vector<size_t> SparseMap;
	typedef std::vector<size_t> SparsePar;
	typedef std::pair<std::vector<size_t>, BitVector<> > Hybrid;
};


// Specialization of RandomDenseStream
template<size_t bitsize> struct MTrandomInt {
    template<typename M32Twister>
    unsigned __LINBOX_INT32 operator() (M32Twister& MT) const {
        return MT.randomInt();
    }
};    

template<> struct MTrandomInt<64> {
    template<typename M32Twister>
    unsigned __LINBOX_INT64 operator() (M32Twister& MT) const {
        unsigned __LINBOX_INT64 tmp = MT.randomInt();
        tmp <<=32;
        return tmp += MT.randomInt();
    }
};

template <class Endianness>
class RandomDenseStreamGF2 : public VectorStream<BitVector<Endianness> >
{
    public:
	typedef BitVector<Endianness> Vector;

	RandomDenseStreamGF2 (const GF2 &, uint32 seed, size_t n, size_t m = 0)
		: MT (seed), _n (n), _m (m), _j (0)
	{}

	Vector &get (Vector &v) 
	{
		typename Vector::word_iterator i;

		if (_m > 0 && _j++ >= _m)
			return v;

		for (i = v.wordBegin (); i != v.wordEnd (); i++)
			*i = MTrandomInt<WordTraits<typename BitVector<Endianness>::word_type>::bits>()(MT);
                
                const size_t zeroing = WordTraits<typename BitVector<Endianness>::word_type>::bits - (v.size() % WordTraits<typename BitVector<Endianness>::word_type>::bits);
                *(v.wordRbegin()) <<= zeroing;
                *(v.wordRbegin()) >>= zeroing;
		return v;
	}

	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

    private:
	MersenneTwister MT;
	size_t          _n;
	size_t          _m;
	size_t          _j;
};

// Specialization of RandomSparseStream

template <class _Vector = Vector<GF2>::Sparse>
class RandomSparseStreamGF2 : public VectorStream<_Vector>
{
    public:
	typedef GF2 Field;
	typedef _Vector Vector;

	RandomSparseStreamGF2 (const GF2 &, uint32 seed, double p, size_t n, size_t m = 0)
		: MT (seed), _n (n), _m (m), _j (0)
	{ setP (p); }

    	RandomSparseStreamGF2 (const GF2 &F, const GF2RandIter& r, double p, size_t n, size_t m = 0)
		: MT (r.getMT()), _n (n), _m (m), _j (0)
	{ setP (p); }

	Vector &get (Vector &v);

	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

	void setP (double p)
	{
		linbox_check ((p >= 0.0) && (p <= 1.0)); 
		_p = p;
		_1_log_1mp   = 1 / log (1 - _p);
	}

    private:
	MersenneTwister MT;
	size_t _n;
	double _p;
	double _1_log_1mp;
	size_t _m;
	size_t _j;
};

template <class _Vector>
_Vector &RandomSparseStreamGF2<_Vector>::get (_Vector &v)
{
	size_t i = (size_t) -1;
	double val;
	int skip;

	if (_m > 0 && _j++ >= _m)
		return v;

	v.clear ();

	while (1) {
		val = (double) MT.randomDouble ();
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

} // namespace LinBox

#include "linbox/vector/vector-domain-gf2.h"
#include "linbox/matrix/matrix-domain-gf2.h"

#ifdef __LINBOX_HAVE_M4RI
#  include "linbox/matrix/m4ri-matrix.h"
#endif

#endif // __LINBOX_field_gf2_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
