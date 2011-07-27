/* lela/ring/type-wrapper.h
 * Copyright 1999-2005 William J Turner,
 *           2001 Bradford Hovinen
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */
 
#ifndef __LELA_RING_TYPE_WRAPPER_H
#define __LELA_RING_TYPE_WRAPPER_H

#include <typeinfo>
#include <string>
#include <algorithm>
#include <cmath>

#include "lela/lela-config.h"
#include "lela/integer.h"
#include "lela/randiter/type-wrapper.h"
#include "lela/blas/context.h"
#include "lela/util/property.h"

namespace LELA 
{

/** \brief Ring-type which wraps a C++ type
 * \ingroup ring
 *
 * This class takes any C++ type which supports all basic
 * arithmetic-operators (+, -, *, /, %) and wraps it in a ring so that
 * it can be used in LELA as such.
 */

template <class K>
class TypeWrapperRing
{
protected:
	integer _p;
	integer _card;

public:
    
	typedef K Element;    

	typedef TypeWrapperRandIter<K> RandIter;

	TypeWrapperRing (integer q = 0, size_t e = 1)
		: _p (q),
		  _card ((q == 0) ? integer (-1) : pow (q.get_d (), (double) e))
		{}  // assuming q is a prime or zero.

	TypeWrapperRing (const TypeWrapperRing &F)
		: _p (F._p),
		  _card (F._card)
		{}
    
	~TypeWrapperRing () {}
    
	const TypeWrapperRing &operator=(const TypeWrapperRing &F) const { return *this; }

	Element &init (Element &x, const integer &y = 0) const 
		{ return x = (const Element &) (static_cast<const Element> (y.get_ui ())); }

	Element &init (Element &x, double y) const 
		{ return x = Element (y); }

	Element &init (Element &x, int y) const 
		{ return x = Element (y); }

	template <class Iterator, class Accessor>
	Element &init (Property<Iterator, Accessor> &x, int y) const 
		{ return init (x.ref (), y); }

	template <class T>
	T &copy (T &x, const Element &y) const
		{ return x = y; }

	template <class Iterator, class Accessor>
	Element &copy (Property<Iterator, Accessor> &x, const Element &y) const
		{ return copy (x.ref (), y); }

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

	template <class T>
	T &addin (T &x, const Element &y) const
		{ return x += y; }
    
	template <class T>
	T &subin (T &x, const Element &y) const
		{ return x -= y; }
    
	Element &mulin (Element &x, const Element &y) const
		{ return x *= y; }

	template <class Iterator, class Accessor>
	Element &mulin (Property<Iterator, Accessor> &x, const Element &y) const
		{ return mulin (x.ref (), y); }

	template <class T>
	T &mulin (T &x, const Element &y) const
		{ return x *= y; }

	bool divin (Element &x, const Element &y) const
		{ if (!isZero (y)) { x /= y; return true; } else return false; }

	template <class T>
	T &negin (T &x) const
		{ return x = - x; }
    
	bool invin (Element &x) const
		{ if (!isZero (x)) { x = Element (1) / x; return true; } else return false; }

	template <class T>
	T &axpyin (T &y, const Element &a, const Element &x) const
		{ return y += T (a) * T (x); }

	std::ostream &write (std::ostream &os) const
		{ return os << "type-wrapper-ring"; }
    
	std::istream &read (std::istream &is) const
		{ return is; }
    
	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }
    
	std::istream &read (std::istream &is, Element &x) const
		{ return is >> x; }

	size_t elementWidth () const
		{ return 10; }

	Element zero () const { return Element (0); }
	Element one () const { return Element (1); }
	Element minusOne () const { return Element (-1); }
    
	/** @name Implementation-Specific Methods.
	 * These methods are not required of all LELA rings
	 * and are included only for the implementation of this ring
	 * template.
	 */
	//@{

	/** Constructor from ring object.
	 * @param  A unparameterized ring object
	 */
	TypeWrapperRing (const K &A) {} 
    
	/** Constant access operator.
	 * @return constant reference to ring object
	 */
	const K &operator () () const { return Element (); }
    
	/** Access operator.
	 * @return reference to ring object
	 */
	K &operator () () { return Element (); }
    
	//@} Implementation-Specific Methods

}; // template <class K> class TypeWrapperRing

#ifdef __LELA_BLAS_AVAILABLE

template <class Element>
struct BLASModule : public GenericModule<TypeWrapperRing<Element> >
{
	struct Tag { typedef typename GenericModule<TypeWrapperRing<Element> >::Tag Parent; };

	BLASModule (const TypeWrapperRing<Element> &R) : GenericModule<TypeWrapperRing<Element> > (R) {}
};

template <>
struct AllModules<TypeWrapperRing<float> > : public BLASModule<float>
{
	struct Tag { typedef BLASModule<float>::Tag Parent; };

	AllModules (const TypeWrapperRing<float> &R) : BLASModule<float> (R) {}
};

template <>
struct AllModules<TypeWrapperRing<double> > : public BLASModule<double>
{
	struct Tag { typedef BLASModule<double>::Tag Parent; };

	AllModules (const TypeWrapperRing<double> &R) : BLASModule<double> (R) {}
};

#endif // __LELA_BLAS_AVAILABLE

} // namespace LELA

#ifdef __LELA_BLAS_AVAILABLE
#  include "lela/blas/level1-cblas.h"
#  include "lela/blas/level2-cblas.h"
#  include "lela/blas/level3-cblas.h"
#endif

#include "lela/blas/level1-generic.tcc"
#include "lela/blas/level2-generic.tcc"
#include "lela/blas/level3-generic.tcc"

#endif // __LELA_RING_TYPE_WRAPPER_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
