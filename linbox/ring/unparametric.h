/* linbox/ring/unparametric.h
 * Copyright 1999-2005 William J Turner,
 *           2001 Bradford Hovinen
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
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

template <class K>
class UnparametricRing
{
protected:
	integer _p;
	integer _card;

public:
    
	typedef K Element;    

	typedef UnparametricRandIter<K> RandIter;

	UnparametricRing (integer q = 0, size_t e = 1)
		: _p (q),
		  _card ((q == 0) ? integer (-1) : pow (q.get_d (), e))
		{}  // assuming q is a prime or zero.

	UnparametricRing (const UnparametricRing &F)
		: _p (F._p),
		  _card (F._card)
		{}
    
	~UnparametricRing () {}
    
	const UnparametricRing &operator=(const UnparametricRing &F) const { return *this; }

	Element &init (Element &x, const integer &y = 0) const 
		{ return x = (const Element &) (static_cast<const Element> (y.get_ui ())); }

	Element &init (Element &x, const double &y) const 
		{ return x = (const Element &) (y); }
   
	integer &convert (integer &x, const Element &y) const 
	{
		Element temp (y);
		return x = static_cast<integer> (temp); 
	}
    
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

	Element &assign (Element &x, const Element &y) const
		{ return x = y; }

	integer &cardinality (integer &c) const
		{ return c = _card; }
    
	integer &characteristic (integer &c) const
		{ return c = _p; }

	bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	bool isZero (const Element &x) const
		{ return x == Element (0); }

	bool isOne (const Element &x) const
		{ return x == Element (1); }

	Element &add (Element &x, const Element &y, const Element &z) const
		{ return x = y + z; }
    
	Element &sub (Element &x, const Element &y, const Element &z) const
		{ return x = y - z; }
    
	Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = y * z; }
    
	bool div (Element &x, const Element &y, const Element &z) const
		{ if (!isZero (y)) { x = y / z; return true; } else return false; }
    
	Element &neg (Element &x, const Element &y) const { return x = - y; }
    
	bool inv (Element &x, const Element &y) const 
		{ if (!isZero (y)) { x = Element (1) / y; return true; } else return false; }
    
	Element &axpy (Element &z, const Element &a, const Element &x, const Element &y) const
		{ return z = a * x + y; }

	Element &addin (Element &x, const Element &y) const
		{ return x += y; }
    
	Element &subin (Element &x, const Element &y) const
		{ return x -= y; }
    
	Element &mulin (Element &x, const Element &y) const
		{ return x *= y; }
    
	bool divin (Element &x, const Element &y) const
		{ if (!isZero (y)) { x /= y; return true; } else return false; }
    
	Element &negin (Element &x) const
		{ return x = - x; }
    
	bool invin (Element &x) const
		{ if (!isZero (x)) { x = Element (1) / x; return true; } else return false; }
    
	Element &axpyin (Element &y, const Element &a, const Element &x) const
		{ return y += a * x; }

	std::ostream &write (std::ostream &os) const
		{ return os << "unparameterized ring(" << sizeof (Element) <<',' << typeid (Element).name() << ')'; }
    
	std::istream &read (std::istream &is) const
		{ return is; }
    
	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }
    
	std::istream &read (std::istream &is, Element &x) const
		{ return is >> x; }

	Element zero () const { return Element (0); }
	Element one () const { return Element (1); }
	Element minusOne () const { return Element (-1); }
    
	/** @name Implementation-Specific Methods.
	 * These methods are not required of all LinBox rings
	 * and are included only for the implementation of this ring
	 * template.
	 */
	//@{

	/** Constructor from ring object.
	 * @param  A unparameterized ring object
	 */
	UnparametricRing (const K &A) {} 
    
	/** Constant access operator.
	 * @return constant reference to ring object
	 */
	const K &operator () () const { return Element (); }
    
	/** Access operator.
	 * @return reference to ring object
	 */
	K &operator () () { return Element (); }
    
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

template <class Element>
struct BLASModule : public GenericModule<UnparametricRing<Element> >
{
	struct Tag { typedef typename GenericModule<UnparametricRing<Element> >::Tag Parent; };

	BLASModule (const UnparametricRing<Element> &R) : GenericModule<UnparametricRing<Element> > (R) {}
};

template <>
struct AllModules<UnparametricRing<float> > : public BLASModule<float>
{
	struct Tag { typedef BLASModule<float>::Tag Parent; };

	AllModules (const UnparametricRing<float> &R) : BLASModule<float> (R) {}
};

template <>
struct AllModules<UnparametricRing<double> > : public BLASModule<double>
{
	struct Tag { typedef BLASModule<double>::Tag Parent; };

	AllModules (const UnparametricRing<double> &R) : BLASModule<double> (R) {}
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
