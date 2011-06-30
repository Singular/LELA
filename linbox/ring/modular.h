/* linbox/ring/modular.h
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

#ifndef __LINBOX_ring_modular_H
#define __LINBOX_ring_modular_H

#include <iostream>
#include <climits>
#include <cmath>
#include <vector>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/util/debug.h"
#include "linbox/util/property.h"
#include "linbox/blas/context.h"
#include "linbox/ring/traits.h"
#include "linbox/randiter/nonzero.h"

namespace LinBox 
{

template <class Element>
class Modular;

template <class Ring>
struct ClassifyRing; 

template <class Element>
struct ClassifyRing<Modular<Element> >
{
	typedef RingCategories::ModularTag categoryTag;
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
class ModularBase 
{
public:

	/*- Element type
	 */
	typedef _Element Element;

	/*- Random iterator generator type.
	 * It must meet the common object interface of random element generators
	 * as given in the the archetype RandIterArchetype.
	 */
	class RandIter;

	/*- @name Object Management
	 */
	//@{
 
	/*- Default constructor.
	 */
	ModularBase (void) {}

	/*- Constructor from an element type.
	 * Sets the modulus of the ring throug the static member of the 
	 * element type.
	 * @param modulus constant reference to integer prime modulus
	 */
	ModularBase (unsigned long modulus) : _modulus (modulus) {}

	/*- Constructor from an integer.
	 * Sets the modulus of the ring throug the static member of the 
	 * element type.
	 * @param modulus constant reference to integer prime modulus
	 */
	ModularBase (const integer &modulus) : _modulus (modulus) {}

	/*- Copy constructor.
	 * Constructs Modular object by copying the ring.
	 * This is required to allow ring objects to be passed by value
	 * into functions.
	 * @param  F Modular object.
	 */
	ModularBase (const ModularBase<Element> &F) : _modulus (F._modulus) {}
 
	/*- Conversion of ring base element to a template class T.
	 * This function assumes the output ring base element x has already been
	 * constructed, but that it is not already initialized.
	 * @return reference to template class T.
	 * @param x template class T to contain output (reference returned).
	 * @param y constant ring base element.
	 */
	integer &convert (integer &x, const Element &y) const
		{ return x = y; }

	double &convert (double &x, const Element &y) const
		{return  x = (double) y;}
 
	float &convert (float &x, const Element &y) const
		{return  x = (float) y;}
	
	/*- Assignment of one ring base element to another.
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &assign (Element &x, Element y) const { return x = y; }

	/*- Cardinality.
	 * Return integer representing cardinality of the domain.
	 * Returns a non-negative integer for all domains with finite
	 * cardinality, and returns -1 to signify a domain of infinite
	 * cardinality.
	 * @return integer representing cardinality of the domain
	 */
	integer &cardinality (integer &c) const
		{ return c = _modulus; }

	/*- Characteristic.
	 * Return integer representing characteristic of the domain.
	 * Returns a positive integer to all domains with finite characteristic,
	 * and returns 0 to signify a domain of infinite characteristic.
	 * @return integer representing characteristic of the domain.
	 */
	integer &characteristic (integer &c) const
		{ return c = _modulus; }

	//@} Object Management

	/*- @name Arithmetic Operations
	 * x <- y op z; x <- op y
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized ring base elements will
	 * give undefined results.
	 */
	//@{

	/*- Equality of two elements.
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return boolean true if equal, false if not.
	 * @param  x ring base element
	 * @param  y ring base element
	 */
	bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	/*- Zero equality.
	 * Test if ring base element is equal to zero.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals zero, false if not.
	 * @param  x ring base element.
	 */
	bool isZero (const Element &x) const
		{ return x == 0; }
 
	/*- One equality.
	 * Test if ring base element is equal to one.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals one, false if not.
	 * @param  x ring base element.
	 */
	bool isOne (const Element &x) const
		{ return x == 1; }


	//@} Arithmetic Operations

	/*- @name Input/Output Operations */
	//@{

	/*- Print ring.
	 * @return output stream to which ring is written.
	 * @param  os  output stream to which ring is written.
	 */
	std::ostream &write (std::ostream &os) const 
		{ return os << "Modular ring, mod " << _modulus; }

	/*- Read ring.
	 * @return input stream from which ring is read.
	 * @param  is  input stream from which ring is read.
	 */
	std::istream &read (std::istream &is) { return is >> _modulus; }


	/*- Print ring base element.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return output stream to which ring base element is written.
	 * @param  os  output stream to which ring base element is written.
	 * @param  x   ring base element.
	 */
	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << integer (x); }

	/*- Read ring base element.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return input stream from which ring base element is read.
	 * @param  is  input stream from which ring base element is read.
	 * @param  x   ring base element.
	 */
	std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;

		x = abs (tmp.get_si ()) % _modulus;
		if (tmp < 0) x = _modulus - x;

		return is; 
	}

	//@}

	/// Return the zero-element of the ring
	Element zero () const { return 0; }

	/// Return the one-element of the ring
	Element one () const { return 1; }

	/// Return the negative of the one-element of the ring
	Element minusOne () const { return _modulus - 1; }

	/// Modulus
	Element _modulus;

}; // class ModularBase

/* .. such comments as here should be on specialization...
 * @param element Element type, e.g. long or integer
 * @param Intermediate Type to use for intermediate computations. This
 *                     should be a data type that can support integers
 *                     twice the length of the maximal modulus used.
 *
 * The primality of the modulus will not be checked, so it is the
 * programmer's responsibility to supply a prime modulus. This class
 * implements a ring of unparameterized integers modulo a prime
 * integer. Ring has (non-static) member to contain modulus of ring.
 */

/** @brief Prime rings of positive characteristic implemented directly in LinBox.
 * 
 * This parameterized ring can be used to construct prime
 * rings. Typical use would be Modular<integer> for integers modulo a
 * large prime, Modular<uint32>, modular<int>, or modular<double> for
 * integers modulo a wordsize prime.  Each of those has specialized
 * performance features suitable to certain applications.
 */
template <class _Element>
class Modular : public ModularBase<_Element>
{
    public:
	typedef _Element Element;
	typedef typename ModularBase<_Element>::RandIter RandIter;

	/*- @name Object Management
	 * @brief see \ref{RingArchetype} for member specs.
	 */
	//@{
 
	//private:
	/*- Default constructor.
	 */
	Modular () {}

	/*- Constructor from an element type
	 * Sets the modulus of the ring throug the static member of the 
	 * element type.
	 * @param modulus constant reference to integer prime modulus
	 */
	Modular (unsigned long modulus) : ModularBase<_Element> (modulus) {}

	/*- Constructor from an integer
	 * Sets the modulus of the ring throug the static member of the 
	 * element type.
	 * @param modulus constant reference to integer prime modulus
	 */
	Modular (const integer &modulus) : ModularBase<_Element> (modulus) {}

	/* Assignment operator
	 * Required by the archetype
	 *
	 * @param F constant reference to Modular object
	 * @return reference to Modular object for self
	 */
	const Modular &operator=(const Modular &F) 
	{
		ModularBase<Element>::_modulus = F._modulus;
		return *this;
	}

	static inline Element getMaxModulus()
		{ return Element ((1ULL << (sizeof(Element) * 8 - 1)) - 1); } 


	/*- Initialization of ring base element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output ring base element x has already been
	 * constructed, but that it is not already initialized.
	 * This is not a specialization of the template function because
	 * such a specialization is not allowed inside the class declaration.
	 * @return reference to ring base element.
	 * @param x ring base element to contain output (reference returned).
	 * @param y integer.
	 */
	Element &init (Element &x, const integer &y = 0) const
	{ 
		x = y % ModularBase<Element>::_modulus;
		if (x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}

	Element &init (Element &x, const int y) const
	{ 
		x = y % ModularBase<Element>::_modulus;
		if (x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}

	Element &init (Element &x, const unsigned int y) const
	{ 
		x = y % ModularBase<Element>::_modulus;
		if (x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}

	/*- Initialization of ring base element from a double.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output ring base element x has already been
	 * constructed, but that it is not already initialized.
	 * This is not a specialization of the template function because
	 * such a specialization is not allowed inside the class declaration.
	 * @return reference to ring base element.
	 * @param x ring base element to contain output (reference returned).
	 * @param y integer.
	 */
	Element &init (Element &x, const double &y) const
	{ 
		double z = fmod(y, (double)ModularBase<Element>::_modulus);
		if (z < 0) z += (double) ModularBase<Element>::_modulus;
		return x = (Element) (z+.5);
	}

	Element &init (Element &x, const float &y) const
	{ 
		float z = fmod(y, (float)ModularBase<Element>::_modulus);
		if (z < 0) z += (float) ModularBase<Element>::_modulus;
		return x = (Element) (z+.5);
	}

	//@}  
	/*- @name Arithmetic Operations
	 * @brief see \ref{RingArchetype} for member specs.
	 * x <- y op z; x <- op y
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized ring base elements will
	 * give undefined results.
	 */
	//@{

	/*- Addition.
	 * x = y + z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	Element &add (Element &x, const Element &y, const Element &z) const
	{
		x = y + z;
		if (x >= ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
		return x;
	}
 
	/* Subtraction.
	 * x = y - z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	Element &sub (Element &x, const Element &y, const Element &z) const
	{ 
		x = y - z;
		if (x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}
 
	/* Multiplication.
	 * x = y * z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = (y * z) % ModularBase<Element>::_modulus; }
 
	/* Division.
	 * x = y / z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	Element &div (Element &x, const Element &y, const Element &z) const
	{ 
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	/* Additive Inverse (Negation).
	 * x = - y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &neg (Element &x, const Element &y) const
		{ if (y == 0) return x = y; else return x = ModularBase<Element>::_modulus - y; }
 
	/* Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &inv (Element &x, const Element &y) const
	{
		// The extended Euclidean algoritm
		Element x_int, y_int, q, tx, ty, temp;
		x_int = ModularBase<Element>::_modulus; 
		y_int = y;
		tx = 0; 
		ty = 1;

		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int;  y_int  = x_int  - q*y_int;  x_int  = temp;
			temp = ty; ty = tx - q*ty; tx = temp;
		}

		// now x_int = gcd (modulus,residue)
		x = tx;
		if (x < 0) x += ModularBase<Element>::_modulus;

		return x;
	}

	/* Natural AXPY.
	 * r  = a * x + y
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 * @param  y ring element.
	 */
	Element &axpy (Element &r, 
		       const Element &a, 
		       const Element &x, 
		       const Element &y) const
	{ 
		r = (a * x + y) % ModularBase<Element>::_modulus;
		if (r < 0) r += ModularBase<Element>::_modulus;
		return r;
	}

	//@} Arithmetic Operations
 
	/*- @name Inplace Arithmetic Operations
	 * @brief see \ref{RingArchetype} for member specs.
	 * x <- x op y; x <- op x
	 */
	//@{

	/*- Inplace Addition.
	 * x += y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &addin (Element &x, const Element &y) const
	{ 
		x += y;
		if (x >= ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
		return x;
	}
 
	/* Inplace Subtraction.
	 * x -= y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &subin (Element &x, const Element &y) const
	{
		x -= y;
		if (x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}
 
	/* Inplace Multiplication.
	 * x *= y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &mulin (Element &x, const Element &y) const
	{
		x *= y;
		x %= ModularBase<Element>::_modulus;
		return x;
	}

	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
	{
		x *= y;
		x %= ModularBase<Element>::_modulus;
		return x;
	}
 
	/* Inplace Division.
	 * x /= y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	Element &divin (Element &x, const Element &y) const
	{
		Element temp;
		inv (temp, y);
		return mulin (x, temp);
	}
 
	/* Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 */
	Element &negin (Element &x) const
		{ if (x == 0) return x; else return x = ModularBase<Element>::_modulus - x; }
 
	/* Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the ring base elementhas already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 */
	Element &invin (Element &x) const
		{ return inv (x, x); }

	/* Inplace AXPY.
	 * r  += a * x
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * Purely virtual
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 */
	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{ 
		r = (r + a * x) % ModularBase<Element>::_modulus;
		if (r < 0) r += ModularBase<Element>::_modulus;
		return r;
	}

	//@} Inplace Arithmetic Operations

}; // class Modular

/** @brief Allows compact storage when the modulus is less than 2^8. 
	
    Requires 1 < the modulus < 2^8, normally prime.
    See RingArchetype for member specifications.
*/
template <>
class Modular<uint8> : public ModularBase<uint8>
{
    public:

	typedef uint8 Element;

	Modular () : _k (0) {}
	Modular (uint32 modulus)
		: ModularBase<uint8> (modulus),
		  _k (((uint64) -1LL) / ((modulus - 1) * (modulus - 1))),
		  _pinv (1.0 / (double) ((uint8) modulus))
		{}
	Modular (const integer &modulus)
		: ModularBase<uint8> (modulus.get_ui ()),
		  _k (((uint64) -1LL) / (((uint8) modulus.get_ui () - 1) * ((uint8) modulus.get_ui () - 1))),
		  _pinv (1.0 / (double) (modulus.get_ui ()))
		{}

	const Modular &operator=(const Modular &F) 
	{
		ModularBase<uint8>::_modulus = F._modulus;
		_k = F._k;
		_pinv = F._pinv;
		return *this;
	}

	Element &init (Element &x, const integer &y = 0) const
	{
		x = abs (y.get_ui ()) % ModularBase<Element>::_modulus;
		if (y < 0) x = ModularBase<Element>::_modulus - x;
		return x;
	}

	Element &init (Element &x, const double &y) const
	{ 
		double z = fmod(y, (double)_modulus);
		if (z < 0) z += (double) _modulus;
		return x = (Element) (z);
	}

	Element &add (Element &x, const Element &y, const Element &z) const
	{
		uint32 t = (uint32) y + (uint32) z;
		if (t >= (uint32) ModularBase<Element>::_modulus) t -= ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &sub (Element &x, const Element &y, const Element &z) const
	{ 
		int32 t = (int32) y - (int32) z;
		if (t < 0) t += ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = ((uint32) y * (uint32) z) % (uint32) ModularBase<Element>::_modulus; }
 
	Element &div (Element &x, const Element &y, const Element &z) const
	{ 
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	Element &neg (Element &x, const Element &y) const
		{ if (y == 0) return x = y; else return x = ModularBase<Element>::_modulus - y; }
 
	Element &inv (Element &x, const Element &y) const
	{
		// The extended Euclidean algoritm 
		int32 x_int, y_int, q, tx, ty, temp;
		x_int = ModularBase<Element>::_modulus;
		y_int = y;
		tx = 0; 
		ty = 1;

		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int; y_int = x_int - q * y_int;
			x_int = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}

		if (tx < 0) tx += ModularBase<Element>::_modulus;

		// now x_int = gcd (modulus,residue)
		return x = tx;
	}

	Element &axpy (Element &r, 
		       const Element &a, 
		       const Element &x, 
		       const Element &y) const
	{
		r = ((uint32) a * (uint32) x + (uint32) y) % (uint32) ModularBase<Element>::_modulus;
		return r;
	}

	Element &addin (Element &x, const Element &y) const
	{ 
		uint32 t = (long) x + (long) y;
		if (t >= (uint32) ModularBase<Element>::_modulus) t -= ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &subin (Element &x, const Element &y) const
	{
		long t = x - y;
		if (t < 0) t += ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &mulin (Element &x, const Element &y) const
	{
		x = ((uint32) x * (uint32) y) % (uint32) ModularBase<Element>::_modulus;
		return x;
	}
 
	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
	{
		x = ((uint32) x * (uint32) y) % (uint32) ModularBase<Element>::_modulus;
		return x;
	}

	Element &divin (Element &x, const Element &y) const
	{
		Element temp;
		inv (temp, y);
		return mulin (x, temp);
	}
 
	Element &negin (Element &x) const
		{ if (x == 0) return x; else return x = ModularBase<Element>::_modulus - x; }
 
	Element &invin (Element &x) const
		{ return inv (x, x); }

	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{ 
		r = ((uint32) r + (uint32) a * (uint32) x) % (uint32) ModularBase<Element>::_modulus;
		return r;
	}

	static inline Element getMaxModulus()
		{ return Element ((1ULL << (sizeof(Element) * 8 - 1)) - 1); } 

	// Number of times one can perform an axpy into a long long
	// before modding out is mandatory.
	size_t _k;

    private:

	// Inverse of modulus in floating point
	double _pinv;

}; // class Modular<uint8>

/** @brief Specialization of class Modular for uint16 element type */
template <>
class Modular<uint16> : public ModularBase<uint16>
{
    public:

	typedef uint16 Element;

	Modular () : _k (0) {}
	Modular (uint32 modulus)
		: ModularBase<uint16> (modulus),
		  _k (((uint64) -1LL) / ((ModularBase<Element>::_modulus - 1) * (ModularBase<Element>::_modulus - 1))),
		  _pinv (1.0 / (double) ((uint16) ModularBase<Element>::_modulus))
		{}
	Modular (const integer &modulus)
		: ModularBase<uint16> (modulus.get_ui ()),
		  _k (((uint64) -1LL) / ((ModularBase<Element>::_modulus - 1) * (ModularBase<Element>::_modulus - 1))),
		  _pinv (1.0 / (double) ((uint16) ModularBase<Element>::_modulus))
		{}

	const Modular &operator = (const Modular &F) 
	{
		ModularBase<Element>::_modulus = F._modulus;
		_k = F._k;
		_pinv = F._pinv;
		return *this;
	}

	Element &init (Element &x, const integer &y = 0) const
	{
		x = abs (y.get_ui ()) % ModularBase<Element>::_modulus;
		if (y < 0) x = ModularBase<Element>::_modulus - x;
		return x;
	}

	Element &init (Element &x, const double &y) const
	{ 
		double z = fmod(y, (double)_modulus);
		if (z < 0) z += (double) _modulus;
		return x = (Element) (z);
	}

	Element &add (Element &x, const Element &y, const Element &z) const
	{
		uint32 t = (uint32) y + (uint32) z;
		if (t >= (uint32) ModularBase<Element>::_modulus) t -= ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &sub (Element &x, const Element &y, const Element &z) const
	{ 
		int32 t = (int32) y - (int32) z;
		if (t < 0) t += ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = ((uint32) y * (uint32) z) % (uint32) ModularBase<Element>::_modulus; }
 
	Element &div (Element &x, const Element &y, const Element &z) const
	{ 
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	Element &neg (Element &x, const Element &y) const
		{ if (y == 0) return x = y; else return x = ModularBase<Element>::_modulus - y; }
 
	Element &inv (Element &x, const Element &y) const
	{
		// The extended Euclidean algoritm 
		int32 x_int, y_int, q, tx, ty, temp;
		x_int = ModularBase<Element>::_modulus;
		y_int = y;
		tx = 0; 
		ty = 1;

		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int; y_int = x_int - q * y_int;
			x_int = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}

		if (tx < 0) tx += ModularBase<Element>::_modulus;

		// now x_int = gcd (modulus,residue)
		return x = tx;
	}

	Element &axpy (Element &r, 
		       const Element &a, 
		       const Element &x, 
		       const Element &y) const
	{
		r = ((uint32) a * (uint32) x + (uint32) y) % (uint32) ModularBase<Element>::_modulus;
		return r;
	}

	Element &addin (Element &x, const Element &y) const
	{ 
		uint32 t = (long) x + (long) y;
		if (t >= (uint32) ModularBase<Element>::_modulus) t -= ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &subin (Element &x, const Element &y) const
	{
		long t = x - y;
		if (t < 0) t += ModularBase<Element>::_modulus;
		return x = t;
	}
 
	Element &mulin (Element &x, const Element &y) const
	{
		x = ((uint32) x * (uint32) y) % (uint32) ModularBase<Element>::_modulus;
		return x;
	}

	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
	{
		x = ((uint32) x * (uint32) y) % (uint32) ModularBase<Element>::_modulus;
		return x;
	}
 
	Element &divin (Element &x, const Element &y) const
	{
		Element temp;
		inv (temp, y);
		return mulin (x, temp);
	}
 
	Element &negin (Element &x) const
		{ if (x == 0) return x; else return x = ModularBase<Element>::_modulus - x; }
 
	Element &invin (Element &x) const
		{ return inv (x, x); }

	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{ 
		r = ((uint32) r + (uint32) a * (uint32) x) % (uint32) ModularBase<Element>::_modulus;
		return r;
	}

	static inline Element getMaxModulus()
		{ return Element ((1ULL << (sizeof(Element) * 8 - 1)) - 1); } 

	// Number of times one can perform an axpy into a long long
	// before modding out is mandatory.
	size_t _k;

    private:

	// Inverse of modulus in floating point
	double _pinv;

}; // class Modular<uint16>

/** @brief Specialization of class Modular for uint32 element type */
template <>
class Modular<uint32> : public ModularBase<uint32>
{
    public:

	typedef uint32 Element;

	Modular () {}
	Modular (uint32 modulus)  : ModularBase<uint32> (modulus) { init_two_64 (); }
	Modular (const integer &modulus)
		: ModularBase<uint32> (modulus.get_ui ())
		{ init_two_64 (); }

	const Modular &operator=(const Modular &F) 
	{
		ModularBase<Element>::_modulus = F._modulus;
		_two_64 = F._two_64;
		return *this;
	}

	Element &init (Element &x, const integer &y = 0) const
	{
		x = abs (y.get_ui ()) % ModularBase<Element>::_modulus;
		if (y < 0) x = ModularBase<Element>::_modulus - x;
		return x;
	}

	Element &init (Element &x, const double &y) const
	{ 
		double z = fmod(y, (double)_modulus);
		if (z < 0) z += (double) _modulus;
		return x = (Element) (z);
	}

	Element &add (Element &x, const Element &y, const Element &z) const
	{
		x = y + z;
		if ((uint32) x >= (uint32) ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
		return x;
	}
 
	Element &sub (Element &x, const Element &y, const Element &z) const
	{
		x = y - z;
		if ((int32) x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}
 
	Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = ((uint64) y * (uint64) z) % (uint64) ModularBase<Element>::_modulus; }
 
	Element &div (Element &x, const Element &y, const Element &z) const
	{ 
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	Element &neg (Element &x, const Element &y) const
		{ if (y == 0) return x = y; else return x = ModularBase<Element>::_modulus - y; }
 
	Element &inv (Element &x, const Element &y) const
	{
		// The extended Euclidean algoritm
		int64 x_int, y_int, q, tx, ty, temp;
		x_int = ModularBase<Element>::_modulus;
		y_int = y;
		tx = 0; 
		ty = 1;

		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int;  y_int  = x_int  - q * y_int;
			x_int  = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}

		if (tx < 0) tx += ModularBase<Element>::_modulus;

		// now x_int = gcd (modulus,residue)
		return x = tx;
	}

	Element &axpy (Element &r, 
		       const Element &a, 
		       const Element &x, 
		       const Element &y) const
	{
		r = ((uint64) a * (uint64) x + (uint64) y) % (uint64) ModularBase<Element>::_modulus;
		if ((int32) r < 0) r += ModularBase<Element>::_modulus;
		return r;
	}

	Element &addin (Element &x, const Element &y) const
	{ 
		x += y;
		if ((uint32) x >= (uint32) ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
		return x;
	}
 
	Element &subin (Element &x, const Element &y) const
	{
		x -= y;
		if ((int32) x < 0) x += ModularBase<Element>::_modulus;
		return x;
	}
 
	Element &mulin (Element &x, const Element &y) const
	{
		x = ((uint64) x * (uint64) y) % (uint64) ModularBase<Element>::_modulus;
		return x;
	}
 
	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
	{
		x = ((uint64) x * (uint64) y) % (uint64) ModularBase<Element>::_modulus;
		return x;
	}

	Element &divin (Element &x, const Element &y) const
	{
		Element temp;
		inv (temp, y);
		return mulin (x, temp);
	}
 
	Element &negin (Element &x) const
		{ if (x == 0) return x; else return x = ModularBase<Element>::_modulus - x; }
 
	Element &invin (Element &x) const
		{ return inv (x, x); }

	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{ 
		r = ((uint64) r + (uint64) a * (uint64) x) % (uint64) ModularBase<Element>::_modulus;
		if ((int32) r < 0) r += ModularBase<Element>::_modulus;
		return r;
	}

	static inline Element getMaxModulus()
		{ return Element ((1ULL << (sizeof(Element) * 8 - 1)) - 1); } 

	Element _two_64;

    private:

	void init_two_64 () 
	{
		uint64 two_64 = 2;

		for (int i = 0; i < 6; ++i)
			two_64 = (two_64 * two_64) % ModularBase<Element>::_modulus;

		_two_64 = two_64;
	}

}; // class Modular<uint32>

template <class Element>
class ModularRandIter;

template <>
class Modular<float>
{

public:

	float modulus;
	unsigned long lmodulus;

	typedef float Element;
	typedef ModularRandIter<float> RandIter;
	typedef NonzeroRandIter<Modular<float>, ModularRandIter<float> > NonZeroRandIter;

	static ClassifyRing<Modular<float> >::categoryTag getCategory ()
		{ return ClassifyRing<Modular<float> >::categoryTag (); }

	Modular () {}

	Modular (int32 p, int exp = 1)
		: modulus ((float) p), lmodulus (p)
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");

		integer max;

		if (modulus > RingTraits<Modular<float> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (float p)
		: modulus (p), lmodulus ((unsigned long) p)
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (modulus > RingTraits<Modular<float> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (long int p)
		: modulus ((float) p), lmodulus (p)
	{
		if ((float) modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if ((float) modulus > RingTraits<Modular<float> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (const integer &p)
		: modulus (p.get_d ()),
		  lmodulus (p.get_ui ())
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (modulus > RingTraits<Modular<float> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (const Modular<float> &mf)
		: modulus (mf.modulus),
		  lmodulus (mf.lmodulus)
		{}

	const Modular &operator=(const Modular<float> &F)
	{
		modulus = F.modulus;
		lmodulus = F.lmodulus;
		return *this;
	}

	integer &cardinality (integer &c) const
		{ return c = integer (modulus); }

	integer &characteristic (integer &c) const
		{ return c = integer (modulus); }

	integer &convert (integer &x, const Element &y) const
		{ return x = integer (y); }

	float &convert (float &x, const Element &y) const
		{ return x = y; }
		
	std::ostream &write (std::ostream &os) const
		{ return os << "float mod " << (int) modulus; }

	std::istream &read (std::istream &is)
	{
		is >> modulus;

		if (modulus <= 1) 
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		if (modulus > 94906265) 
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");

		return is;
	}
		
	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }

	std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;
		init (x, tmp);

		return is;
	}

	Element &init (Element &x, const integer &y) const
		{ return x = (Element) (y.get_ui () % lmodulus); }

	inline Element &init (Element &x, float y = 0.0) const
	{
		x = fmod (y, modulus);
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &init (Element &x, double y) const
	{
		x = fmod (y, double (modulus));
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &init (Element &x, unsigned long y) const
	{
		x = fmod (float (y), modulus);
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &init (Element &x, int y) const
	{
		x = fmod (float (y), modulus);
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &init (Element &x, unsigned int y) const
	{
		x = fmod (float (y), modulus);
		if (x < 0) x += modulus;
		return x;
	}
		
	inline Element &assign (Element &x, const Element &y) const
		{ return x = y; }

	inline bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	inline  bool isZero (const Element &x) const
		{ return x == 0.; }

	inline bool isOne (const Element &x) const
		{ return x == 1.; }

	inline Element &add (Element &x, const Element &y, const Element &z) const
	{
		x = y + z;
		if (x >= modulus) x -= modulus;
		return x;
	}
 
	inline Element &sub (Element &x, const Element &y, const Element &z) const
	{
		x = y - z;
		if (x < 0) x += modulus;
		return x;
	}
		
	inline Element &mul (Element &x, const Element &y, const Element &z) const
	{
		float tmp = y * z;
		x = fmod (tmp, modulus);
		return x;
	}
 
	template <class Iterator>
	inline Property<Iterator> &mul (Property<Iterator> &x, const Element &y, const Element &z) const
	{
		float tmp = y * z;
		x = fmod (tmp, modulus);
		return x;
	}

	inline Element &div (Element &x, const Element &y, const Element &z) const
	{
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	inline Element &neg (Element &x, const Element &y) const
	{
		if (y == 0)
			return x = 0;
		else 
			return x = modulus - y;
	}
 
	inline Element &inv (Element &x, const Element &y) const
	{
		// The extended Euclidean algoritm 
		int x_int, y_int, q, tx, ty, temp;

		x_int = int (modulus);
		y_int = int (y);
		tx = 0; 
		ty = 1;
		  
		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int; y_int = x_int - q * y_int;
			x_int = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}

		if (tx < 0) tx += (int) modulus;

		// now x_int = gcd (modulus,residue)
		return x = (float) tx;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const
	{
		float tmp = a * x + y;
		return r = fmod (tmp, modulus); 
	}

	inline Element &addin (Element &x, const Element &y) const
	{
		x += y;
		if (x >= modulus) x -= modulus;
		return x;
	}
 
	inline Element &subin (Element &x, const Element &y) const
	{
		x -= y;
		if (x < 0.) x += modulus;
		return x;
	}
 
	inline Element &mulin (Element &x, const Element &y) const
		{ return mul (x, x, y); }

	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
		{ return mul (x, x, y); }
 
	inline Element &divin (Element &x, const Element &y) const
		{ return div (x, x, y); }
 
	inline Element &negin (Element &x) const
	{
		if (x == 0.)
			return x; 
		else
			return x = modulus - x; 
	}
		
	inline Element &invin (Element &x) const
		{ return inv (x, x); }
		
	inline Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		float tmp = r + a * x;
		return r = fmod (tmp, modulus); 
	}

	static inline float getMaxModulus()
		{ return 4096.0; } // floor( 2^12 )

	float zero () const { return 0.0; }
	float one () const { return 1.0; }
	float minusOne () const { return modulus - 1.0; }
}; // class Modular<float>

template <>
class Modular<double>
{

public:	       

	double modulus;
	unsigned long lmodulus;

	typedef double Element;

	typedef ModularRandIter<double> RandIter;
	typedef NonzeroRandIter<Modular<double>, ModularRandIter<double> > NonZeroRandIter;

	static ClassifyRing<Modular<double> >::categoryTag getCategory ()
		{ return ClassifyRing<Modular<double> >::categoryTag (); }

	Modular () {}

	Modular (int32 p, int exp = 1)
		: modulus ((double) p),
		  lmodulus (p)
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");

		integer max;

		if (modulus > RingTraits<Modular<double> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (double p)
		: modulus (p),
		  lmodulus ((unsigned long) p)
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (modulus > RingTraits<Modular<double> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (long int p)
		: modulus ((double) p),
		  lmodulus (p)
	{
		if ((double) modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if ((double) modulus > RingTraits<Modular<double> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (const integer &p)
		: modulus (p.get_d ()),
		  lmodulus (p.get_d ())
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");
		if (modulus > 94906265)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (const Modular<double> &mf)
		: modulus (mf.modulus),
		  lmodulus (mf.lmodulus)
		{}

	const Modular &operator = (const Modular<double> &F)
	{
		modulus = F.modulus;
		lmodulus = F.lmodulus;
		return *this;
	}

	integer &cardinality (integer &c) const
		{ return c = integer (modulus); }

	integer &characteristic (integer &c) const
		{ return c = integer (modulus); }

	integer &convert (integer &x, const Element &y) const
		{ return x = integer (y); }

	double &convert (double &x, const Element &y) const
		{ return x=y; }
		
	std::ostream &write (std::ostream &os) const
		{ return os << "double mod " << (int) modulus; }

	std::istream &read (std::istream &is)
	{
		is >> modulus;

		if (modulus <= 1) 
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		if (modulus > 94906265) 
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");

		return is;
	}

	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << (int) x; }

	std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;
		init (x, tmp);

		return is;
	}

	Element &init (Element &x, const integer &y) const
		{ return x = (Element) (y.get_ui () % lmodulus); }

	inline Element &init (Element &x, double y = 0) const
	{
		x = fmod (y, modulus);
		if (x < 0) x += modulus;
		return x;
	}
		
	inline Element &assign(Element &x, const Element &y) const
		{ return x = y; }

	inline bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	inline bool isZero (const Element &x) const
		{ return x == 0.; }

	inline bool isOne (const Element &x) const
		{ return x == 1.; }

	inline bool isMinusOne (const Element &x) const
		{ return (x == modulus - 1.); }

	inline Element &add (Element &x, const Element &y, const Element &z) const
	{
		x = y + z;
		if (x >= modulus) x -= modulus;
		return x;
	}
 
	inline Element &sub (Element &x, const Element &y, const Element &z) const
	{
		x = y - z;
		if (x < 0) x += modulus;
		return x;
	}
		
	inline Element &mul (Element &x, const Element &y, const Element &z) const
	{
		double tmp = y * z;
		x = fmod (tmp, modulus);
		return x;
	}

	template <class Iterator>
	inline Property<Iterator> &mul (Property<Iterator> &x, const Element &y, const Element &z) const
	{
		double tmp = y * z;
		x = fmod (tmp, modulus);
		return x;
	}
 
	inline Element &div (Element &x, const Element &y, const Element &z) const
	{
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	inline Element &neg (Element &x, const Element &y) const
	{
		if (y == 0)
			return x = 0;
		else
			return x = modulus - y;
	}
 
	inline Element &inv (Element &x, const Element &y) const
	{
		// The extended Euclidean algoritm 
		int x_int, y_int, q, tx, ty, temp;
		x_int = int (modulus);
		y_int = int (y);
		tx = 0; 
		ty = 1;
		  
		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int; y_int = x_int - q * y_int;
			x_int = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}
		  
		if (tx < 0)
			tx += (int) modulus;
		  
		// now x_int = gcd (modulus,residue)
		return x = (double) tx;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const
	{
		double tmp = a * x + y;
		return r = fmod (tmp, modulus); 
	}

	inline Element &addin (Element &x, const Element &y) const
	{
		x += y;
		if (x >= modulus) x -= modulus;
		return x;
	}
 
	inline Element &subin (Element &x, const Element &y) const
	{
		x -= y;
		if (x < 0.) x += modulus;
		return x;
	}
 
	inline Element &mulin (Element &x, const Element &y) const
		{ return mul (x, x, y); }

	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
		{ return mul (x, x, y); }
 
	inline Element &divin (Element &x, const Element &y) const
		{ return div (x, x, y); }
 
	inline Element &negin (Element &x) const
	{
		if (x == 0.)
			return x; 
		else 
			return x = modulus - x; 
	}
		
	inline Element &invin (Element &x) const
	{
		return inv (x, x);
	}
		
	inline Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		double tmp = r + a * x;
		return r = fmod (tmp, modulus); 
	}

	static inline double getMaxModulus()
		{ return 67108864.0; } // 2^26 

	double zero () const { return 0.0; }
	double one () const { return 1.0; }
	double minusOne () const { return modulus - 1.0; }
}; // class Modular<double>

template <>
inline std::ostream& ModularBase<integer>::write (std::ostream &os) const 
	{ return os << "GMP integers mod " << _modulus; }

template <>
inline integer& Modular<integer>::init (integer& x, const double& y) const 
{
	integer tmp = (integer) y % _modulus;
	if (tmp<0) tmp += _modulus;
	return x = tmp;
}

template <class Element>
struct ZpModule : public GenericModule
{
	struct Tag { typedef GenericModule::Tag Parent; };

	mutable std::vector<uint64> _tmp;
};

template <>
struct ZpModule<float> : public GenericModule
{
	struct Tag { typedef GenericModule::Tag Parent; };

	size_t _nmax;
};

template <>
struct ZpModule<double> : public GenericModule
{
	struct Tag { typedef GenericModule::Tag Parent; };

	size_t _nmax;
};

template <class Element>
struct AllModules<Modular<Element> > : public ZpModule<Element>
{
	struct Tag { typedef typename ZpModule<Element>::Tag Parent; };
};

} // namespace LinBox

#include "linbox/blas/level1-modular.h"
#include "linbox/blas/level2-modular.h"

#include "linbox/blas/level1-generic.h"
#include "linbox/blas/level2-generic.h"
#include "linbox/blas/level3-generic.h"

#include "linbox/blas/level1-modular.tcc"
#include "linbox/blas/level2-modular.tcc"

#include "linbox/randiter/modular.h"

#endif // __LINBOX_ring_modular_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
