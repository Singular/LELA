/* lela/ring/gf2.h
 * Copyright 2003-2007 Bradford Hovinen, Jean-Guillaume Dumas, Clement Pernet
 *
 * Implementation of GF(2)
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RING_GF2_H
#define __LELA_RING_GF2_H

#include <iostream>
#include <climits>
#include <cmath>

#include "lela/lela-config.h"
#include "lela/util/debug.h"
#include "lela/blas/context.h"
#include "lela/integer.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/hybrid.h"

#ifdef __LELA_HAVE_M4RI
#  include "lela/matrix/m4ri-matrix.h"
#else
#  include "lela/algorithms/strassen-winograd.h"
#endif

// Namespace in which all LELA code resides
namespace LELA 
{

class GF2RandIter;

/** 
 * \brief Integers modulo 2
 *
 * This is a tuned implementation of the ring of integers modulo
 * 2. In particular, when one constructs a VectorDomain object over
 * this ring, highly optimized bit operations will be used to make
 * vector arithmetic very fast.
 \ingroup ring
 */

class GF2
{
    public:

	typedef bool Element;

	typedef GF2RandIter RandIter;

	GF2 () {}
	GF2 (int p, int exp = 1)
	{
		if (p != 2)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be 2");

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");
	}

	GF2 (const GF2 &) {}

	const GF2 &operator = (const GF2 &) 
		{ return *this; }

	Element &init (Element &x, const int &y = 0) const { return x = y & 1; }
	Element &init (Element &x, const unsigned int &y = 0) const { return x = y & 1; }
	Element &init (Element &x, const long &y = 0) const { return x = y & 1; }
	Element &init (Element &x, const unsigned long &y = 0) const { return x = y & 1; }
	Element &init (Element &x, const float &y) const { return x = static_cast<unsigned char>(y) & 1; }
	Element &init (Element &x, const double &y) const { return x = static_cast<unsigned char>(y) & 1; }
	Element &init (Element &x, const integer &y) const { return x = y.get_si () & 1; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> init (BitVectorReference<Iterator, Endianness> x, const integer &y = 0) const
		{ return x = y.get_si () & 1; }

	std::_Bit_reference init (std::_Bit_reference x, const integer &y = 0) const
		{ return x = y.get_si () & 1; }

	integer &convert (integer &x, Element y) const { return x = y; }
         std::_Bit_reference convert (std::_Bit_reference x, Element y) const { return x = y; }

	Element &copy (Element &x, Element y) const { return x = y; }

	template <class Iterator, class Endianness>
        BitVectorReference<Iterator, Endianness> copy (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ return x = y; }

        std::_Bit_reference copy (std::_Bit_reference x, Element y) const
		{ return x = y; }
    
	integer &cardinality (integer &c) const { return c = 2; }
	integer &characteristic (integer &c) const { return c = 2; }

	bool areEqual (Element x, Element y) const { return x == y; }
	bool isZero (Element x) const { return !x; }
	bool isOne (Element x) const { return x; }

	std::ostream &write (std::ostream &os) const { return os << "GF(2)"; }
	std::istream &read (std::istream &is) { return is; }

	std::ostream &write (std::ostream &os, Element x) const { return os << x; }
 
	std::istream &read (std::istream &is, Element &x) const
		{ int v; is >> v; if (v == 0) x = false; else x = true; return is; }

	template <class Iterator, class Endianness>
	std::istream &read (std::istream &is, BitVectorReference<Iterator, Endianness> x) const
		{ int v; is >> v; if (v == 0) x = false; else x = true; return is; }

	std::istream &read (std::istream &is, std::_Bit_reference x) const
		{ bool a; is >> a; x=a; return is; }

	size_t elementWidth () const
		{ return 1; }

	Element &add (Element &x, Element y, Element z) const
		{ return x = y ^ z; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> add (BitVectorReference<Iterator, Endianness> x, Element y, Element z) const
		{ return x = y ^ z; }
 
	std::_Bit_reference add (std::_Bit_reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
	Element &sub (Element &x, Element y, Element z) const
		{ return x = y ^ z; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> sub (BitVectorReference<Iterator, Endianness> x, Element y, Element z) const
		{ return x = y ^ z; }
 
	std::_Bit_reference sub (std::_Bit_reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
	Element &mul (Element &x, Element y, Element z) const
		{ return x = y & z; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> mul (BitVectorReference<Iterator, Endianness> x, Element y, Element z) const
		{ return x = y & z; }
 
	std::_Bit_reference mul (std::_Bit_reference x, Element y, Element z) const
		{ return x = y & z; }
 
	bool div (Element &x, Element y, Element z) const
		{ if (z) { x = y; return true; } else return false; }

	template <class Iterator, class Endianness>
	bool div (BitVectorReference<Iterator, Endianness> x, Element y, Element z) const
		{ if (z) { x = y; return true; } else return false; }
 
	bool div (std::_Bit_reference x, Element y, Element z) const
		{ if (z) { x = y; return true; } else return false; }
 
	Element &neg (Element &x, Element y) const
		{ return x = y; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> neg (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ return x = y; }
 
	std::_Bit_reference neg (std::_Bit_reference x, Element y) const
		{ return x = y; }
 
	bool inv (Element &x, Element y) const
		{ if (y) { return x = true; } else return false; }

	template <class Iterator, class Endianness>
	bool inv (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ if (y) { return x = true; } else return false; }

	bool inv (std::_Bit_reference x, Element y) const
		{ if (y) { return x = true; } else return false; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> axpy (BitVectorReference<Iterator, Endianness> r, 
						       Element a, 
						       Element x, 
						       Element y) const
		{ return r = (a & x) ^ y; }

	std::_Bit_reference axpy (std::_Bit_reference r, 
				  Element a, 
				  Element x, 
				  Element y) const
		{ return r = (a & x) ^ y; }

	Element &axpy (Element &r, Element a, Element x, Element y) const
		{ return r = (a & x) ^ y; }

	Element &addin (Element &x, Element y) const
		{ return x ^= y; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> addin (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ return x ^= y; }
 
	std::_Bit_reference addin (std::_Bit_reference x, Element y) const
        	{ return x = x ^ y; }
 
	Element &subin (Element &x, Element y) const
		{ return x ^= y; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> subin (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ return x ^= y; }
 
	std::_Bit_reference subin (std::_Bit_reference x, Element y) const
		{ return x = x ^ y; }
 
	Element &mulin (Element &x, Element y) const
		{ return x &= y; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> mulin (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ return x &= y; }
    
	Element &mulin (std::_Bit_reference &x, Element y) const
		{ return mulin ((bool &) x, y); }
 
	bool divin (Element &x, Element y) const
		{ return y; }

	template <class Iterator, class Endianness>
	bool divin (BitVectorReference<Iterator, Endianness> x, Element y) const
		{ return y; }
 
	bool divin (std::_Bit_reference x, Element y) const
		{ return y; }
 
	Element &negin (Element &x) const
		{ return x; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> negin (BitVectorReference<Iterator, Endianness> x) const
		{ return x; }
 
	std::_Bit_reference negin (std::_Bit_reference x) const
		{ return x; }
 
	bool invin (Element &x) const
		{ return x; }

	template <class Iterator, class Endianness>
	bool invin (BitVectorReference<Iterator, Endianness> x) const
		{ return x; }

	bool invin (std::_Bit_reference x) const
		{ return x; }

	Element &axpyin (Element &r, Element a, Element x) const
		{ return r ^= a & x; }

	template <class Iterator, class Endianness>
	BitVectorReference<Iterator, Endianness> axpyin (BitVectorReference<Iterator, Endianness> r, Element a, Element x) const
		{ return r ^= a & x; }

	std::_Bit_reference axpyin (std::_Bit_reference r, Element a, Element x) const
		{ return r = r ^ (a & x); }

	Element &axpyin (Element &r, const std::_Bit_reference a, Element x) const
		{ return r ^= a & x; }

	std::_Bit_reference axpyin (std::_Bit_reference r, const std::_Bit_reference a, Element x) const
		{ return r = r ^ (a & x); }

	Element &axpyin (Element &r, Element a, const std::_Bit_reference x) const
		{ return r ^= a & static_cast<bool>(x); }

	std::_Bit_reference axpyin (std::_Bit_reference r, Element a, const std::_Bit_reference x) const
		{ return r = r ^ (a & static_cast<bool>(x)); }

	Element &axpyin (Element &r, const std::_Bit_reference a, const std::_Bit_reference x) const
		{ return r ^= a & static_cast<bool>(x); }

	std::_Bit_reference axpyin (std::_Bit_reference r, const std::_Bit_reference a, const std::_Bit_reference x) const
		{ return r = r ^ (a & static_cast<bool>(x)); }

	static inline int getMaxModulus() { return 2; }

	Element zero () const { return false; }
	Element one () const { return true; }
	Element minusOne () const { return true; }

}; // class GF2

// Specialization of canonical vector types

template <>
struct RawVector<bool>
{
    public:
	typedef BitVector<DefaultEndianness<uint64> > Dense;
	typedef std::vector<uint32> Sparse;
	typedef HybridVector<DefaultEndianness<uint64>, uint16, uint64> Hybrid;
};

template <>
struct Vector<GF2>
{
    public:
	typedef bool Element;
	typedef BitVector<DefaultEndianness<uint64> > Dense;
	typedef std::vector<uint32> Sparse;
	typedef HybridVector<DefaultEndianness<uint64>, uint16, uint64> Hybrid;
};

// Calculation-modules

#ifdef __LELA_HAVE_M4RI

struct M4RIModule : public GenericModule<GF2>
{
	struct Tag { typedef GenericModule<GF2>::Tag Parent; };

	M4RIMatrix _tmp;

	M4RIModule (const GF2 &R) : GenericModule<GF2> (R) {}
};

template <>
struct AllModules<GF2> : public M4RIModule
{
	struct Tag { typedef M4RIModule::Tag Parent; };

	AllModules (const GF2 &R) : M4RIModule (R) {}
};

#else // !__LELA_HAVE_M4RI

template <>
struct AllModules<GF2> : public StrassenModule<GF2, GenericModule<GF2> >
{
	struct Tag { typedef StrassenModule<GF2, GenericModule<GF2> >::Tag Parent; };

	AllModules (const GF2 &R) : StrassenModule<GF2, GenericModule<GF2> > (R) {}
};

#endif // __LELA_HAVE_M4RI

} // namespace LELA

#include "lela/blas/level1-gf2.h"
#include "lela/blas/level2-gf2.h"

#ifdef __LELA_HAVE_M4RI
#  include "lela/blas/level3-m4ri.h"
#endif // __LELA_HAVE_M4RI

#include "lela/blas/level1-generic.h"
#include "lela/blas/level2-generic.h"
#include "lela/blas/level3-generic.h"

#include "lela/blas/level1-gf2.tcc"
#include "lela/blas/level2-gf2.tcc"

#ifdef __LELA_HAVE_M4RI
#  include "lela/blas/level3-m4ri.tcc"
#endif // __LELA_HAVE_M4RI

#include "lela/randiter/gf2.h"

#endif // __LELA_RING_GF2_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
