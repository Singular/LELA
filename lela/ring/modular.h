/* lela/ring/modular.h
 * Copyright 1999-2001 William J Turner,
 *           2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Integers modulo n
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RING_MODULAR_H
#define __LELA_RING_MODULAR_H

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
#include "lela/ring/type-wrapper.h"

#define FLOAT_MANTISSA 24
#define DOUBLE_MANTISSA 53

namespace LELA 
{

/** Traits for the modular ring
 *
 * This contains parameters for computations with the modular ring
 * which depend on the element-type
 *
 * \ingroup ring
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
		{ r = a % m; shift_up (r, m); return r; }

	/// Initialise an element from an integer; for initialisation of the modulus
	static Element &init_modulus (Element &elt, integer x)
		{ elt = x; return elt; }

	/// Write an element to the given stream
	static std::ostream &write (std::ostream &os, const Element &x)
		{ return os << x; }

	/// Get the width of a typical element of the ring
	static size_t element_width (const integer &modulus)
		{ return (size_t) ceil (log (modulus.get_d ()) / M_LN10); }

	template <class T>
	static size_t element_width (const T &modulus)
		{ return (size_t) ceil (log (double (modulus)) / M_LN10); }

	/// If v is too low for the window of validity, shift it up by the value of the modulus
	template <class T>
	static T &shift_up (T &v, const Element &modulus)
		{ if (v < 0) v += modulus; return v; }

	/// If v is too high for the window of validity, shift it down by the value of the modulus
	template <class T>
	static T &shift_down (T &v, const Element &modulus)
		{ if (v >= modulus) v -= modulus; return v; }

	/// Move v to within the valid range for the modulus
	template <class T>
	static T &valid_rep (T &v, const Element &modulus)
		{ return shift_down (v, modulus); }

	/// Negate the given element
	template <class T>
	static T &neg (T &v, const Element &modulus)
		{ v = modulus - v; return valid_rep (v, modulus); }

	/// Subtract the one element from the other
	template <class T>
	static T &sub (T &v, const Element &y, const Element &modulus)
		{ v += modulus - y; return valid_rep (v, modulus); }

	/// Get what is in any event a positive representation of the given element
	template <class T>
	static T &positive_rep (T &v, const Element &modulus)
		{ return v; }
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
		{ integer t = (integer) a % (integer) m; shift_up (t, m); return r = t.get_ui (); }
	static Element &reduce (Element &r, FatElement a, Element m) 
		{ return r = a % m; }
	static DoubleFatElement &reduce (DoubleFatElement &r, DoubleFatElement a, Element m) 
		{ return r = a % m; }
	static Element &reduce (Element &r, int a, Element m) 
		{ int t = a % (int) m; if (t < 0) t += m; return r = t; }
	static Element &init_modulus (Element &elt, integer x)
		{ elt = x.get_ui (); return elt; }
	static std::ostream &write (std::ostream &os, const Element &x)
		{ return os << (int) x; }
	static size_t element_width (Element modulus)
		{ return (size_t) ceil (log (double (modulus)) / M_LN10); }
	template <class T>
	static T &shift_up (T &v, uint8 modulus)
		{ if (v < 0) v += modulus; return v; }
	template <class T>
	static T &shift_down (T &v, uint8 modulus)
		{ if (v >= modulus) v -= modulus; return v; }
	template <class T>
	static T &valid_rep (T &v, const Element &modulus)
		{ return shift_down (v, modulus); }
	template <class T>
	static T &neg (T &v, const Element &modulus)
		{ v = modulus - v; return valid_rep (v, modulus); }
	template <class T>
	static T &sub (T &v, const Element &y, const Element &modulus)
		{ v += modulus - y; return valid_rep (v, modulus); }
	template <class T>
	static T &positive_rep (T &v, const Element &modulus)
		{ return v; }
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
		{ integer t = (integer) a % (integer) m; shift_up (t, m); return r = t.get_ui (); }
	static Element &reduce (Element &r, FatElement a, Element m) 
		{ return r = a % m; }
	static Element &reduce (Element &r, DoubleFatElement a, Element m) 
		{ return r = a % m; }
	static DoubleFatElement &reduce (DoubleFatElement &r, DoubleFatElement a, Element m) 
		{ return r = a % m; }
	static Element &reduce (Element &r, int a, Element m) 
		{ int t = a % (int) m; if (t < 0) t += m; return r = t; }
	static Element &init_modulus (Element &elt, integer x)
		{ elt = x.get_ui (); return elt; }
	static std::ostream &write (std::ostream &os, const Element &x)
		{ return os << x; }
	static size_t element_width (Element modulus)
		{ return (size_t) ceil (log (double (modulus)) / M_LN10); }
	template <class T>
	static T &shift_up (T &v, uint16 modulus)
		{ if (v < 0) v += modulus; return v; }
	template <class T>
	static T &shift_down (T &v, uint16 modulus)
		{ if (v >= modulus) v -= modulus; return v; }
	template <class T>
	static T &neg (T &v, const Element &modulus)
		{ v = modulus - v; return valid_rep (v, modulus); }
	template <class T>
	static T &sub (T &v, const Element &y, const Element &modulus)
		{ v += modulus - y; return valid_rep (v, modulus); }
	template <class T>
	static T &valid_rep (T &v, const Element &modulus)
		{ return shift_down (v, modulus); }
	template <class T>
	static T &positive_rep (T &v, const Element &modulus)
		{ return v; }
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
	static Element &reduce (Element &r, FatElement a, Element m) 
		{ return r = a % m; }
	static DoubleFatElement &reduce (DoubleFatElement &r, DoubleFatElement a, Element m) 
		{ return r = a % m; }
	static Element &reduce (Element &r, int a, Element m) 
		{ long long t = (long long) a % (long long) m; shift_up (t, m); return r = t; }
	static Element &init_modulus (Element &elt, integer x)
		{ elt = x.get_ui (); return elt; }
	static std::ostream &write (std::ostream &os, const Element &x)
		{ return os << x; }
	static size_t element_width (Element modulus)
		{ return (size_t) ceil (log (double (modulus)) / M_LN10); }
	template <class T>
	static T &shift_up (T &v, uint32 modulus)
		{ if (v < 0) v += modulus; return v; }
	template <class T>
	static T &shift_down (T &v, uint32 modulus)
		{ if ((FatElement) v >= (FatElement) modulus) v -= modulus; return v; }
	template <class T>
	static T &neg (T &v, const Element &modulus)
		{ v = modulus - v; return valid_rep (v, modulus); }
	template <class T>
	static T &sub (T &v, const Element &y, const Element &modulus)
		{ v += modulus - y; return valid_rep (v, modulus); }
	template <class T>
	static T &valid_rep (T &v, const Element &modulus)
		{ return shift_down (v, modulus); }
	template <class T>
	static T &positive_rep (T &v, const Element &modulus)
		{ return v; }
};

// Specialisation for float
template <>
struct ModularTraits<float>
{
	typedef float Element;
	typedef float FatElement;
	typedef float DoubleFatElement;
	typedef int EEAElement;
	static bool valid_modulus (const integer &modulus) { return modulus.get_d () < 4096.0; }
	template <class FE>
	static Element &reduce (Element &r, const FE &a, Element m) 
		{ integer t = (integer) a % (integer) m; shift_down (t, m); return r = t.get_d (); }
	static Element &reduce (Element &r, FatElement a, Element m) 
		{ r = fmod (a, m); return shift_down (shift_up (r, m), m); }
	static double &reduce (double &r, FatElement a, Element m) 
		{ r = fmod (a, m); return shift_down (shift_up (r, m), m); }
	static Element &init_modulus (Element &elt, integer x)
		{ elt = x.get_d (); return elt; }
	static std::ostream &write (std::ostream &os, const Element &x)
		{ return os << (int) x; }
	static size_t element_width (Element modulus)
		{ return (size_t) ceil (log (modulus / 2) / M_LN10) + 1; }
	template <class T>
	static T &shift_up (T &v, float modulus)
		{ if (v < -modulus / 2) v += modulus; return v; }
	template <class T>
	static T &shift_down (T &v, float modulus)
		{ if (v > modulus / 2) v -= modulus; return v; }
	template <class T>
	static T &neg (T &v, const Element &modulus)
		{ return v = -v; }
	template <class T>
	static T &sub (T &v, const Element &y, const Element &modulus)
		{ v -= y; return valid_rep (v, modulus); }
	template <class T>
	static T &valid_rep (T &v, const Element &modulus)
		{ return shift_down (shift_up (v, modulus), modulus); }
	template <class T>
	static T positive_rep (T v, const Element &modulus)
		{ if (v < 0) return v + modulus; else return v; }
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
		{ integer t = (integer) a % (integer) m; shift_up (shift_down (t, m), m); return r = t.get_d (); }
	static Element &reduce (Element &r, FatElement a, Element m) 
		{ r = fmod (a, m); return shift_up (shift_down (r, m), m); }
	static Element &init_modulus (Element &elt, integer x)
		{ elt = x.get_d (); return elt; }
	static std::ostream &write (std::ostream &os, const Element &x)
		{ return os << (long long) x; }
	static size_t element_width (Element modulus)
		{ return (size_t) ceil (log (modulus / 2) / M_LN10) + 1; }
	template <class T>
	static T &shift_up (T &v, double modulus)
		{ if (v < -modulus / 2) v += modulus; return v; }
	template <class T>
	static T &shift_down (T &v, double modulus)
		{ if (v > modulus / 2) v -= modulus; return v; }
	template <class T>
	static T &neg (T &v, const Element &modulus)
		{ return v = -v; }
	template <class T>
	static T &sub (T &v, const Element &y, const Element &modulus)
		{ v -= y; return valid_rep (v, modulus); }
	template <class T>
	static T &valid_rep (T &v, const Element &modulus)
		{ return shift_up (shift_down (v, modulus), modulus); }
	template <class T>
	static T positive_rep (T v, const Element &modulus)
		{ if (v < 0) return v + modulus; else return v; }
};

/** Integers modulo n
 * 
 * \ingroup ring
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
		ModularTraits<Element>::init_modulus (_modulus, modulus);
	}

	Modular (const Modular<Element> &F) : _modulus (F._modulus) {}

	integer &convert (integer &x, const Element &y) const { return x = y; }
	double &convert (double &x, const Element &y) const {return  x = (double) y;}
	float &convert (float &x, const Element &y) const {return  x = (float) y;}

	template <class Iterator, class Accessor, class T>
	Element &init (Property<Iterator, Accessor> x, const T &y) const
		{ return init (x.ref (), y); }

	template <class T>
	Element &init (Element &x, const T &y) const
	{
		ModularTraits<Element>::reduce (x, y, _modulus);
		return x;
	}

	Element &init (Element &x, double y) const
	{ 
		double z = fmod (y, (double) _modulus);
		ModularTraits<Element>::shift_up (z, _modulus);
		return x = (Element) (z+.5);
	}

	Element &init (Element &x, float y) const
	{ 
		float z = fmod (y, (float) _modulus);
		ModularTraits<Element>::shift_up (z, _modulus);
		return x = (Element) (z+.5);
	}

	Element &copy (Element &x, Element y) const
		{ return x = y; }

	template <class Iterator, class Accessor>
	Element &copy (Property<Iterator, Accessor> x, Element y) const
		{ return copy (x.ref (), y); }

	integer &cardinality (integer &c) const { return c = _modulus; }
	integer &characteristic (integer &c) const { return c = _modulus; }

	bool areEqual (const Element &x, const Element &y) const { return x == y; }
	bool isZero (const Element &x) const { return x == 0; }
	bool isOne (const Element &x) const { return x == 1; }

	std::ostream &write (std::ostream &os) const { os << "ZZ/"; return ModularTraits<Element>::write (os, _modulus); }
	std::istream &read (std::istream &is) { return is >> _modulus; }

	std::ostream &write (std::ostream &os, const Element &x) const { return ModularTraits<Element>::write (os, x); }

	std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;

		ModularTraits<Element>::reduce (x, tmp.get_si (), _modulus);

		return is; 
	}

	size_t elementWidth () const
		{ return ModularTraits<Element>::element_width (_modulus); }

	Element &add (Element &x, const Element &y, const Element &z) const
	{
		typename ModularTraits<Element>::FatElement ty = y;
		ty += z;
		ModularTraits<Element>::valid_rep (ty, _modulus);
		return x = ty;
	}
 
	Element &sub (Element &x, const Element &y, const Element &z) const
	{ 
		typename ModularTraits<Element>::FatElement ty = y;
		ModularTraits<Element>::sub (ty, z, _modulus);
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

		eea (gcd, a, ModularTraits<Element>::positive_rep (z, _modulus), b, _modulus);

		if (ModularTraits<Element>::reduce (yremgcd, y, Element (gcd)) == 0) {
			ModularTraits<Element>::shift_down (ModularTraits<Element>::shift_up (a, _modulus), _modulus);
			mul (x, Element (a), y / Element (gcd));
			return true;
		} else
			return false;
	}
 
	Element &neg (Element &x, const Element &y) const
		{ x = y; return ModularTraits<Element>::neg (x, _modulus); }
 
	bool inv (Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::EEAElement ty, tm;
		typename ModularTraits<Element>::EEAElement gcd;

		eea (gcd, ty, ModularTraits<Element>::positive_rep (y, _modulus), tm, _modulus);

		if (gcd == 1) {
			ModularTraits<Element>::shift_up (ty, _modulus);
			x = ty;
			ModularTraits<Element>::shift_down (x, _modulus);
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
		ModularTraits<Element>::valid_rep (tx, _modulus);
		return x = tx;
	}

	Element &subin (Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::FatElement tx = x;
		ModularTraits<Element>::sub (tx, y, _modulus);
		return x = tx;
	}
 
	Element &mulin (Element &x, const Element &y) const
	{
		typename ModularTraits<Element>::FatElement tx = x;
		tx *= y;
		return ModularTraits<Element>::reduce (x, tx, _modulus);
	}

	template <class Iterator, class Accessor>
	Element &mulin (Property<Iterator, Accessor> &x, const Element &y) const
		{ return mulin (x.ref (), y); }
 
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
		{ return ModularTraits<Element>::neg (x, _modulus); }
 
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
	Element minusOne () const { Element t = _modulus - 1; return ModularTraits<Element>::shift_down (t, _modulus); }

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
	struct Tag
	{
		typedef typename GenericModule<Modular<Element> >::Tag Parent;
		typedef typename AllModules<TypeWrapperRing<Element> >::Tag TWParent;
	};

	/// Number of times a product of two elements can be added before it is necessary to reduce by the modulus; 0 for unlimited
	size_t block_size;

	/// Modules for the switch over to TypeWrapperRing
	AllModules<TypeWrapperRing<Element> > TWM;

	mutable std::vector<typename ModularTraits<Element>::DoubleFatElement> _tmp;

	ZpModule (const Modular<Element> &R)
		: block_size (((typename ModularTraits<Element>::DoubleFatElement) -1LL) / ((R._modulus - 1) * (R._modulus - 1))),
		  TWM (TypeWrapperRing<Element> ())
		{}
};

template <>
struct ZpModule<integer> : public GenericModule<Modular<integer> >
{
	struct Tag
	{
		typedef GenericModule<Modular<integer> >::Tag Parent;
		typedef AllModules<TypeWrapperRing<integer> >::Tag TWParent;
	};

	/// Number of times a product of two elements can be added before it is necessary to reduce by the modulus; 0 for unlimited
	size_t block_size;

	/// Modules for the switch over to TypeWrapperRing
	AllModules<TypeWrapperRing<integer> > TWM;

	mutable std::vector<integer> _tmp;

	ZpModule (const Modular<integer> &R) : block_size (0), TWM (TypeWrapperRing<integer> ()) {}
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

template <>
struct ZpModule<float> : public GenericModule<Modular<float> >
{
	struct Tag
	{
		typedef GenericModule<Modular<float> >::Tag Parent;
		typedef AllModules<TypeWrapperRing<float> >::Tag TWParent;
	};

	/// Number of times a product of two elements can be added before it is necessary to reduce by the modulus; 0 for unlimited
	size_t block_size;

	/// Modules for the switch over to TypeWrapperRing
	AllModules<TypeWrapperRing<float> > TWM;

	mutable std::vector<ModularTraits<float>::DoubleFatElement> _tmp;

	ZpModule (const Modular<float> &R)
		: block_size (floor (float (1 << FLOAT_MANTISSA) / ((R._modulus - 1) * (R._modulus - 1)))),
		  TWM (TypeWrapperRing<float> ()) {}
};

template <>
struct ZpModule<double> : public GenericModule<Modular<double> >
{
	struct Tag
	{
		typedef GenericModule<Modular<double> >::Tag Parent;
		typedef AllModules<TypeWrapperRing<double> >::Tag TWParent;
	};

	/// Number of times a product of two elements can be added before it is necessary to reduce by the modulus; 0 for unlimited
	size_t block_size;

	/// Modules for the switch over to TypeWrapperRing
	AllModules<TypeWrapperRing<double> > TWM;

	mutable std::vector<ModularTraits<double>::DoubleFatElement> _tmp;

	ZpModule (const Modular<double> &R)
		: block_size (floor (double (1ULL << DOUBLE_MANTISSA) / ((R._modulus - 1) * (R._modulus - 1)))),
		  TWM (TypeWrapperRing<double> ()) {}
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
#include "lela/blas/level3-modular.h"

#include "lela/blas/level1-generic.h"
#include "lela/blas/level2-generic.h"
#include "lela/blas/level3-generic.h"

#include "lela/blas/level1-modular.tcc"
#include "lela/blas/level2-modular.tcc"
#include "lela/blas/level3-modular.tcc"

#include "lela/blas/level3-sw.h"

#include "lela/randiter/modular.h"

#endif // __LELA_RING_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
