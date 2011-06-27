/* linbox/ring/unparametric.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
 
#ifndef __LINBOX_ring_unparametric_H
#define __LINBOX_ring_unparametric_H
#include <typeinfo>

#include <string>
#include <algorithm>
#include <cmath>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/randiter/unparametric.h"
#include "linbox/ring/ring-interface.h"
#include "linbox/ring/traits.h"
#include "linbox/blas/context.h"

namespace LinBox 
{

/** \brief Unparameterized ring adapter.
    \ingroup ring

    A ring having an interface similar to that of floats is adapted to LinBox.

    Used to generate efficient ring classes for unparameterized rings (or hidden parameter rings).
	 
    Some rings are implemented by definition of the C++ arithmetic operators, such as z = x*y, 
    for z, y, z instances of a type K.   The LinBox ring 
    Unparametric<K> is the adaptation to LinBox.

    For a typical unparametric ring, some of the methods must be defined in a specialization. 

*/

template <class Ring>
struct ClassifyRing;

template <class K>
class UnparametricRing;

template <class K>
struct ClassifyRing<UnparametricRing<K> >
{
	typedef RingCategories::GenericTag categoryTag;
};

template <class K>
class UnparametricRing : public RingInterface
{
protected:
	integer _p;
	integer _card;

public:
    
	/** @name Common Object Interface for a LinBox Ring.
	 * These methods and member types are required of all LinBox rings.
	 * See \ref{RingArchetype} for detailed specifications.
	 */
	//@{
    
	/** The ring's element type.
	 * Type K must provide a default constructor, 
	 * a copy constructor, a destructor, and an assignment operator.
	 */

	typedef K Element;    

	/// Type of random ring element generators.
	typedef UnparametricRandIter<K> RandIter;

	/** @name Ring Object Basics.
	 */
	//@{
    
	/** Builds this ring to have characteristic q and cardinality q<sup>e</sup>.
	 *  This constructor must be defined in a specialization.
	 */
	UnparametricRing (integer q = 0, size_t e = 1)
		: _p (q),
		  _card ((q == 0) ? integer (-1) : pow (q.get_d (), e))
		{}  // assuming q is a prime or zero.

	/// construct this ring as copy of F.
	UnparametricRing (const UnparametricRing &F)
		: _p (F._p),
		  _card (F._card)
		{}
    
	/// 
	~UnparametricRing () {}
    
	/* Assignment operator.
	 * Assigns UnparametricRing object F to ring.
	 * @param  F UnparametricRing object.
	 */
	const UnparametricRing &operator=(const UnparametricRing &F) const { return *this; }

	//@} Ring Object Basics.
    
	/** @name Data Object Management.
	 * first argument is set and the value is also returned.
	 */

	//@{

	/// x := y.  Caution: it is via cast to long.  Good candidate for specialization.
	Element &init (Element &x, const integer &y = 0) const 
		{ return x = (const Element &) (static_cast<const Element> (y.get_ui ())); }

	Element &init (Element &x, const double &y) const 
		{ return x = (const Element &) (y); }
   
	/// x :=  y.  Caution: it is via cast to long.  Good candidate for specialization.
	integer &convert (integer &x, const Element &y) const 
	{
		Element temp (y);
		return x = static_cast<integer> (temp); 
	}
    
	/// x :=  y.  Caution: it is via cast to long.  Good candidate for specialization. --dpritcha
	double &convert (double &x, const Element &y) const 
	{ 
		Element temp (y);
		return x = static_cast<double> (temp); 
	}
 
	float &convert (float &x, const Element &y) const 
	{ 
		Element temp (y);
		return x = static_cast<float> (temp); 
	}

	///
	Element &assign (Element &x, const Element &y) const
		{ return x = y; }

	/// c := cardinality of this ring (-1 if infinite).
	integer &cardinality (integer &c) const
		{ return c = _card; }
    
	/// c := characteristic of this ring (zero or prime).
	integer &characteristic (integer &c) const
		{ return c = _p; }

	//@} Data Object Management
    
	/// @name Comparison Predicates
	//@{
	///  x == y
	bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	///  x == 0
	bool isZero (const Element &x) const
		{ return x == Element (0); }

	///  x == 1
	bool isOne (const Element &x) const
		{ return x == Element (1); }
	//@} Comparison Predicates
		 
	/** @name Arithmetic Operations 
	 * The first argument is set and is also the return value.
	 */
	//@{
    
	/// x := y + z
	Element &add (Element &x, const Element &y, const Element &z) const
		{ return x = y + z; }
    
	/// x := y - z
	Element &sub (Element &x, const Element &y, const Element &z) const
		{ return x = y - z; }
    
	/// x := y*z
	Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = y * z; }
    
	/// x := y/z
	Element &div (Element &x, const Element &y, const Element &z) const
		{ return x = y / z; }
    
	/// x := -y
	Element &neg (Element &x, const Element &y) const { return x = - y; }
    
	/// x := 1/y
	Element &inv (Element &x, const Element &y) const 
		{ return x = Element (1) / y; }
    
	/// z := a*x + y 
	// more optimal implementation, if available, can be defined in a template specialization.
	Element &axpy (Element &z, 
		       const Element &a, 
		       const Element &x, 
		       const Element &y) const
		{ return z = a * x + y; }
 
	//@} Arithmetic Operations
    
	/** @name Inplace Arithmetic Operations 
	 * The first argument is modified and the result is the return value.
	 */
	//@{
    
	/// x := x + y
	Element &addin (Element &x, const Element &y) const
		{ return x += y; }
    
	/// x := x - y
	Element &subin (Element &x, const Element &y) const
		{ return x -= y; }
    
	/// x := x*y
	Element &mulin (Element &x, const Element &y) const
		{ return x *= y; }
    
	/// x := x/y
	Element &divin (Element &x, const Element &y) const
		{ return x /= y; }
    
	/// x := -x
	Element &negin (Element &x) const
		{ return x = - x; }
    
	/// x := 1/x
	Element &invin (Element &x) const
		{ return x = Element (1) / x; }
    
	/// y := a*x + y
	Element &axpyin (Element &y, const Element &a, const Element &x) const
		{ return y += a * x; }
 
	//@} Inplace Arithmetic Operations

	/** @name Input/Output Operations */
	//@{
    
	/** Print ring.
	 * @return output stream to which ring is written.
	 * @param  os  output stream to which ring is written.
	 */
	std::ostream &write (std::ostream &os) const
		{ return os << "unparameterized ring(" << sizeof (Element) <<',' << typeid (Element).name() << ')'; }
    
	/** Read ring.
	 * @return input stream from which ring is read.
	 * @param  is  input stream from which ring is read.
	 */
	std::istream &read (std::istream &is) const
		{ return is; }
    
	/** Print ring element.
	 * @return output stream to which ring element is written.
	 * @param  os  output stream to which ring element is written.
	 * @param  x   ring element.
	 */
	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }
    
	/** Read ring element.
	 * @return input stream from which ring element is read.
	 * @param  is  input stream from which ring element is read.
	 * @param  x   ring element.
	 */
	std::istream &read (std::istream &is, Element &x) const
		{ return is >> x; }
    
	//@}
    
	//@} Common Object Interface
    
	/** @name Implementation-Specific Methods.
	 * These methods are not required of all LinBox rings
	 * and are included only for the implementation of this ring
	 * template.
	 */
	//@{
    
	/// Default constructor
	//UnparametricRing (void) {}
    
	/** Constructor from ring object.
	 * @param  A unparameterized ring object
	 */
	UnparametricRing (const K &A) {} 
    
	/** Constant access operator.
	 * @return constant reference to ring object
	 */
	const K &operator () (void) const { return Element (); }
    
	/** Access operator.
	 * @return reference to ring object
	 */
	K &operator () (void) { return Element (); }
    
	//@} Implementation-Specific Methods

}; // template <class K> class UnparametricRing

template<class Ring>
class RingAXPY;

template<>
class RingAXPY<UnparametricRing<integer> >
{
public:
	typedef UnparametricRing<integer> Ring;
	typedef integer Element;

	/** Constructor.
	 * A faxpy object if constructed from a Ring and a ring element.
	 * Copies of this objects are stored in the faxpy object.
	 * @param F ring F in which arithmetic is done
	 */
	RingAXPY (const Ring &F)
		: _F (F)
		{ _y = 0; }
 
	/** Copy constructor.
	 * @param faxpy
	 */
	RingAXPY (const RingAXPY<Ring> &faxpy)
		: _F (faxpy._F), _y (faxpy._y)
		{}
 
	/** Assignment operator
	 * @param faxpy
	 */
	RingAXPY<Ring> &operator = (const RingAXPY &faxpy)
		{ _y = faxpy._y; return *this; }
 
	/** Add a*x to y
	 * y += a*x.
	 * @param a constant reference to element a
	 * @param x constant reference to element x
	 * allow optimal multiplication, such as integer * int
	 */
	template<class Element1>
	inline Element&  mulacc (const Element &a, const Element1 &x)
		{ return _y += (a * x); }
 
	template<class Element1>
	inline Element& accumulate (const Element1 &t)
		{ return _y += t; }
 
	/** Retrieve y
	 *
	 * Performs the delayed modding out if necessary
	 */
	inline Element &get (Element &y)
		{ y = _y; return y; }
 
	/** Assign method.
	 * Stores new ring element for arithmetic.
	 * @return reference to self
	 * @param y_init constant reference to element a
	 */
	inline RingAXPY &assign (const Element& y)
	{
		_y = y;
		return *this;
	}
		
	inline void reset ()
		{ _y = 0; }
			
private:
 
	/// Ring in which arithmetic is done
	Ring _F;
 
	/// Ring element for arithmetic
	Element _y;

};
	
struct BLASModule : public GenericModule
{
	struct Tag { typedef GenericModule::Tag Parent; };
};

template <>
struct AllModules<UnparametricRing<float> > : public BLASModule
{
	struct Tag { typedef BLASModule::Tag Parent; };
};

template <>
struct AllModules<UnparametricRing<double> > : public BLASModule
{
	struct Tag { typedef BLASModule::Tag Parent; };
};

} // namespace LinBox

#ifdef __LINBOX_BLAS_AVAILABLE
#  include "linbox/blas/level1-cblas.h"
#  include "linbox/blas/level2-cblas.h"
#  include "linbox/blas/level3-cblas.h"
#endif

#include "linbox/blas/level1-generic.tcc"
#include "linbox/blas/level2-generic.tcc"
#include "linbox/blas/level3-generic.tcc"

#endif // __LINBOX_ring_unparametric_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
