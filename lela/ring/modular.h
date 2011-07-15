/* lela/ring/modular.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LELA_ring_modular_H
#define __LELA_ring_modular_H

#include <iostream>
#include <climits>
#include <cmath>
#include <vector>

#include "lela/lela-config.h"
#include "lela/integer.h"
#include "lela/util/debug.h"
#include "lela/util/property.h"
#include "lela/blas/context.h"
#include "lela/randiter/nonzero.h"
#include "lela/algorithms/strassen-winograd.h"

namespace LELA 
{

/** Traits for the modular ring
 *
 * This contains parameters for computations with the modular ring
 * which depend on the element-type
 */
template <class Element>
struct ModularTraits
{
	/// Type to be used as intermediate when multiplying elements
	typedef Element FatElement;

	/// Type to be used as intermediate when computing dot-products
	typedef Element DoubleFatElement;

	/// Type to be used when running extended Euclidean algorithm
	typedef Element EEAElement;

	/// Return true if the given proposed modulus is valid for this element-type, false otherwise
	static bool valid_modulus (const integer &modulus) { return true; }

	/// Function to reduce by modulus
	template <class FE>
	static Element &reduce (Element &r, FE a, Element m) 
		{ r = a % m; if (r < 0) r += m; return r; }

	/// Version of above which writes into a property
	template <class Iterator, class Accessor, class FE>
	static Property<Iterator, Accessor> &reduce (Property<Iterator, Accessor> &r, FE a, Element m)
		{ r = a % m; return r; }

	/// Initialise an element from an integer; for initialisation of the modulus
	static Element &init_element (Element &elt, integer x)
		{ elt = x; return elt; }
};

// Specialisation for uint8
template <>
struct ModularTraits<uint8>
{
	typedef uint8 Element;
	typedef uint32 FatElement;
	typedef uint32 DoubleFatElement;
	typedef int EEAElement;
	static bool valid_modulus (const integer &modulus) { return modulus < integer (1U << 8); }
	template <class FE>
	static Element &reduce (Element &r, const FE &a, Element m) 
		{ integer t = (integer) a % (integer) m; if (t < 0) t += m; return r = t.get_ui (); }
	static Element &reduce (Element &r, int a, Element m) 
		{ int t = a % (int) m; if (t < 0) t += m; return r = t; }
	template <class Iterator, class Accessor, class FE>
	static Property<Iterator, Accessor> &reduce (Property<Iterator, Accessor> &r, FE a, Element m)
		{ r = a % m; return r; }
	static Element &init_element (Element &elt, integer x)
		{ elt = x.get_ui (); return elt; }
};

// Specialisation for uint16
template <>
struct ModularTraits<uint16>
{
	typedef uint16 Element;
	typedef uint32 FatElement;
	typedef uint64 DoubleFatElement;
	typedef int EEAElement;
	static bool valid_modulus (const integer &modulus) { return modulus < integer (1U << 16); }
	template <class FE>
	static Element &reduce (Element &r, const FE &a, Element m) 
		{ integer t = (integer) a % (integer) m; if (t < 0) t += m; return r = t.get_ui (); }
	static Element &reduce (Element &r, int a, Element m) 
		{ int t = a % (int) m; if (t < 0) t += m; return r = t; }
	template <class Iterator, class Accessor, class FE>
	static Property<Iterator, Accessor> &reduce (Property<Iterator, Accessor> &r, FE a, Element m)
		{ r = a % m; return r; }
	static Element &init_element (Element &elt, integer x)
		{ elt = x.get_ui (); return elt; }
};

// Specialisation for uint32
template <>
struct ModularTraits<uint32>
{
	typedef uint32 Element;
	typedef uint64 FatElement;
	typedef uint64 DoubleFatElement;
	typedef int EEAElement;
	static bool valid_modulus (const integer &modulus) { return modulus < integer (1U << 31) * 2; }
	template <class FE>
	static Element &reduce (Element &r, const FE &a, Element m) 
		{ integer t = (integer) a % (integer) m; if (t < 0) t += m; return r = t.get_ui (); }
	static Element &reduce (Element &r, int a, Element m) 
		{ long long t = (long long) a % (long long) m; if (t < 0) t += m; return r = t; }
	template <class Iterator, class Accessor, class FE>
	static Property<Iterator, Accessor> &reduce (Property<Iterator, Accessor> &r, FE a, Element m)
		{ r = a % m; return r; }
	static Element &init_element (Element &elt, integer x)
		{ elt = x.get_ui (); return elt; }
};

// Specialisation for float
template <>
struct ModularTraits<float>
{
	typedef float Element;
	typedef float FatElement;
	typedef double DoubleFatElement;
	typedef int EEAElement;
	static bool valid_modulus (const integer &modulus) { return modulus.get_d () < 4096.0; }
	template <class FE>
	static Element &reduce (Element &r, const FE &a, Element m) 
		{ integer t = (integer) a % (integer) m; if (t < 0) t += m; return r = t.get_d (); }
	template <class Iterator, class Accessor, class FE>
	static Property<Iterator, Accessor> &reduce (Property<Iterator, Accessor> &r, FE a, Element m)
		{ r = fmod (a, m); return r; }
	static Element &init_element (Element &elt, integer x)
		{ elt = x.get_d (); return elt; }
};

// Specialisation for double
template <>
struct ModularTraits<double>
{
	typedef double Element;
	typedef double FatElement;
	typedef double DoubleFatElement;
	typedef int EEAElement;
	static bool valid_modulus (const integer &modulus) { return modulus.get_d () < 67108864.0; }
	template <class FE>
	static Element &reduce (Element &r, const FE &a, Element m) 
		{ integer t = (integer) a % (integer) m; if (t < 0) t += m; return r = t.get_d (); }
	template <class Iterator, class Accessor, class FE>
	static Property<Iterator, Accessor> &reduce (Property<Iterator, Accessor> &r, FE a, Element m)
		{ r = fmod (a, m); return r; }
	static Element &init_element (Element &elt, integer x)
		{ elt = x.get_d (); return elt; }
};

/** @name ModularBase 
 * \brief Base for prime rings where the elements are represented by various primitive types 
 * (and their operations).
 * Normally use it's children.  This class is of interest for the developer of a new ring representation.
 *
 * 
 * This parameterized ring can be used to construct any prime
 * ring. Typical use would be Modular<integer> for integers modulo a
 * large prime, Modular<long, long long> for integers modulo a wordsize
 * prime, etc. for integers modulo a half-wordsize prime.
 \ingroup ring
*/
template <class _Element>
class Modular
{
public:

	typedef _Element Element;

	class RandIter;

	Element _modulus;

	Modular () {}

	Modular (unsigned long modulus) : _modulus (modulus)
		{ lela_check (ModularTraits<Element>::valid_modulus (integer (modulus))); }

	Modular (const integer &modulus)
	{
		lela_check (ModularTraits<Element>::valid_modulus (modulus));
		ModularTraits<Element>::init_element (_modulus, modulus);
	}

	Modular (const Modular<Element> &F) : _modulus (F._modulus) {}

	integer &convert (integer &x, const Element &y) const { return x = y; }
	double &convert (double &x, const Element &y) const {return  x = (double) y;}
	float &convert (float &x, const Element &y) const {return  x = (float) y;}

	template <class T>
	Element &init (Element &x, const T &y) const
	{
		ModularTraits<Element>::reduce (x, y, _modulus);
		return x;
	}

	template <class Iterator, class Accessor, class T>
	Property<Iterator, Accessor> init (Property<Iterator, Accessor> x, const T &y) const
	{
		ModularTraits<Element>::reduce (x, y, _modulus);
		return x;
	}

	Element &init (Element &x, double y) const
	{ 
		double z = fmod (y, (double) _modulus);
		if (z < 0) z += (double) _modulus;
		return x = (Element) (z+.5);
	}

	template <class Iterator, class Accessor>
	Property<Iterator, Accessor> init (Property<Iterator, Accessor> x, double y) const
	{ 
		double z = fmod (y, (double) _modulus);
		if (z < 0) z += (double) _modulus;
		return x = (Element) (z+.5);
	}

	Element &init (Element &x, float y) const
	{ 
		float z = fmod (y, (float) _modulus);
		if (z < 0) z += (float) _modulus;
		return x = (Element) (z+.5);
	}

	template <class Iterator, class Accessor>
	Property<Iterator, Accessor> init (Property<Iterator, Accessor> x, float y) const
	{ 
		float z = fmod (y, (float) _modulus);
		if (z < 0) z += (float) _modulus;
		return x = (Element) (z+.5);
	}

	Element &assign (Element &x, Element y) const
		{ return x = y; }

	template <class Iterator, class Accessor>
	Property<Iterator, Accessor> assign (Property<Iterator, Accessor> x, Element y) const
		{ return x = y; }

	integer &cardinality (integer &c) const { return c = _modulus; }
	integer &characteristic (integer &c) const { return c = _modulus; }

	bool areEqual (const Element &x, const Element &y) const { return x == y; }
	bool isZero (const Element &x) const { return x == 0; }
	bool isOne (const Element &x) const { return x == 1; }

	std::ostream &write (std::ostream &os) const { return os << "ZZ/" << integer (_modulus).get_str (); }
	std::istream &read (std::istream &is) { return is >> _modulus; }

	std::ostream &write (std::ostream &os, const Element &x) const { return os << integer (x).get_str (); }

	std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;

		ModularTraits<Element>::reduce (x, abs (tmp.get_si ()), _modulus);
		if (tmp < 0) x = _modulus - x;

		return is; 
	}

	Element &add (Element &x, const Element &y, const Element &z) const
	{
		typename ModularTraits<Element>::FatElement ty = y;
		ty += z;
		if (ty >= _modulus) ty -= _modulus;
		return x = ty;
	}
 
	Element &sub (Element &x, const Element &y, const Element &z) const
	{ 
		typename ModularTraits<Element>::FatElement ty = y;
		ty += _modulus - z;
		if (ty >= _modulus) ty -= _modulus;
		return x = ty;
	}
 
	Element &mul (Element &x, const Element &y, const Element &z) const
	{
		typename ModularTraits<Element>::FatElement ty = y;
		ty *= z;
		return ModularTraits<Element>::reduce (x, ty, _modulus);
	}
 
	bool div (Element &x, const Element &y, const Element &z) const
	{ 
		typename ModularTraits<Element>::EEAElement a, b;
		typename ModularTraits<Element>::EEAElement gcd;
		Element yremgcd;

		eea (gcd, a, z, b, _modulus);

		if (ModularTraits<Element>::reduce (yremgcd, y, Element (gcd)) == 0) {
			if (a < 0) a += _modulus;
			mul (x, Element (a), y / Element (gcd));
			return true;
		} else
			return false;
	}
 
	Element &neg (Element &x, const Element &y) const
		{ if (y == 0) return x = y; else return x = _modulus - y; }
 
	bool inv (Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::EEAElement ty, tm;
		typename ModularTraits<Element>::EEAElement gcd;

		eea (gcd, ty, y, tm, _modulus);

		if (gcd == 1) {
			if (ty < 0) ty += _modulus;
			x = ty;
			return true;
		} else
			return false;
	}

	Element &axpy (Element &r, const Element &a, const Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::FatElement t = a;
		t *= x; t += y;
		return ModularTraits<Element>::reduce (r, t, _modulus);
	}

	Element &addin (Element &x, const Element &y) const
	{ 
		typename ModularTraits<Element>::FatElement tx = x;
		tx += y;
		if (tx >= _modulus) tx -= _modulus;
		return x = tx;
	}

	Element &subin (Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::FatElement tx = x;
		tx += _modulus - y;
		if (tx >= _modulus) tx -= _modulus;
		return x = tx;
	}
 
	Element &mulin (Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::FatElement tx = x;
		tx *= y;
		return ModularTraits<Element>::reduce (x, tx, _modulus);
	}

	template <class Iterator, class Accessor>
	Property<Iterator, Accessor> &mulin (Property<Iterator, Accessor> &x, const Element &y) const
	{
		typename ModularTraits<Element>::FatElement tx = x;
		tx *= y;
		return ModularTraits<Element>::reduce (x, tx, _modulus);
	}
 
	bool divin (Element &x, const Element &y) const
	{
		Element temp;

		if (!div (temp, x, y))
			return false;
		else {
			x = temp;
			return true;
		}
	}
 
	Element &negin (Element &x) const
		{ if (x == 0) return x; else return x = _modulus - x; }
 
	bool invin (Element &x) const
		{ return inv (x, x); }

	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		typename ModularTraits<Element>::FatElement t = a;
		t *= x; t += r;
		return ModularTraits<Element>::reduce (r, t, _modulus);
	}

	Element zero () const { return 0; }
	Element one () const { return 1; }
	Element minusOne () const { return _modulus - 1; }

private:
	// The extended Euclidean algoritm
	typename ModularTraits<Element>::EEAElement &eea (typename ModularTraits<Element>::EEAElement &gcd,
							  typename ModularTraits<Element>::EEAElement &a,
							  typename ModularTraits<Element>::EEAElement x,
							  typename ModularTraits<Element>::EEAElement &b,
							  typename ModularTraits<Element>::EEAElement y) const
	{
		if (y == 0) {
			gcd = x;
			a = 1;
			b = 0;
		} else {
			typename ModularTraits<Element>::EEAElement q = x / y, r = x % y;
			eea (gcd, b, y, a, r);
			b -= q * a;
		}

		return gcd;
	}

}; // class Modular

template <class Element>
struct ZpModule : public GenericModule<Modular<Element> >
{
	struct Tag { typedef typename GenericModule<Modular<Element> >::Tag Parent; };

	/// Number of times a product of two elements can be added before it is necessary to reduce by the modulus; 0 for unlimited
	size_t block_size;

	mutable std::vector<typename ModularTraits<Element>::DoubleFatElement> _tmp;

	ZpModule (const Modular<Element> &R) : block_size (((typename ModularTraits<Element>::DoubleFatElement) -1LL) / ((R._modulus - 1) * (R._modulus - 1))) {}
};

template <>
struct ZpModule<integer> : public GenericModule<Modular<integer> >
{
	struct Tag { typedef GenericModule<Modular<integer> >::Tag Parent; };

	/// Number of times a product of two elements can be added before it is necessary to reduce by the modulus; 0 for unlimited
	size_t block_size;

	mutable std::vector<integer> _tmp;

	ZpModule (const Modular<integer> &R) : block_size (0) {}
};

template <>
struct ZpModule<uint32> : public GenericModule<Modular<uint32> >
{
	struct Tag { typedef GenericModule<Modular<uint32> >::Tag Parent; };

	/// 2^64 modulo modulus
	uint32 two_64;

	mutable std::vector<ModularTraits<uint32>::DoubleFatElement> _tmp;

	ZpModule (const Modular<uint32> &R) : two_64 (init_two_64 (R._modulus)) {}

private:
	uint32 init_two_64 (uint32 modulus)
	{
		uint64 two_64 = 2;

		for (int i = 0; i < 6; ++i)
			two_64 = (two_64 * two_64) % modulus;

		return two_64;
	}
};

template <class Element>
struct AllModules<Modular<Element> > : public StrassenModule<Modular<Element>, ZpModule<Element> >
{
	struct Tag { typedef typename StrassenModule<Modular<Element>, ZpModule<Element> >::Tag Parent; };

	AllModules (const Modular<Element> &R) : StrassenModule<Modular<Element>, ZpModule<Element> > (R) {}
};

} // namespace LELA

#include "lela/blas/level1-modular.h"
#include "lela/blas/level2-modular.h"

#include "lela/blas/level1-modular.tcc"
#include "lela/blas/level2-modular.tcc"

#include "lela/blas/level3-sw.h"

#include "lela/randiter/modular.h"

#endif // __LELA_ring_modular_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
