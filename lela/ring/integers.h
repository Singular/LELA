/* lela/ring/integers.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Ring of integers, using LELA::integer as an element-type
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */
 
#ifndef __LELA_RING_INTEGERS_H
#define __LELA_RING_INTEGERS_H

#include <typeinfo>
#include <string>
#include <algorithm>
#include <cmath>

#include "lela/lela-config.h"
#include "lela/ring/interface.h"
#include "lela/integer.h"
#include "lela/blas/context.h"
#include "lela/util/property.h"

namespace LELA 
{

// Forward-declaration
class IntegerRandIter;

/** Ring of integers
 *
 * \ingroup ring
 */

class Integers : RingInterface<integer>
{
public:
    
	typedef integer Element;    

	typedef IntegerRandIter RandIter;

	Integers () : _zero (0), _one (1), _minus_one (-1) {}

	Element &init (Element &x, const integer &y = 0) const 
		{ return x = y; }

	Element &init (Element &x, double y) const 
		{ return x = Element (y); }

	Element &init (Element &x, int y) const 
		{ return x = Element (y); }

	template <class Iterator, class Accessor>
	Element &init (Property<Iterator, Accessor> &x, Element &y) const 
		{ return init (x.ref (), y); }

	Element &copy (Element &x, const Element &y) const
		{ return x = y; }

	template <class Iterator, class Accessor>
	Element &copy (Property<Iterator, Accessor> &x, const Element &y) const
		{ return copy (x.ref (), y); }

	integer &cardinality (integer &c) const
		{ return c = 0; }
    
	integer &characteristic (integer &c) const
		{ return c = 0; }

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
		{ if (y % z == 0) { x = y / z; return true; } else return false; }
    
	Element &neg (Element &x, const Element &y) const { return x = - y; }
    
	bool inv (Element &x, const Element &y) const 
		{ if (y == 1 || y == -1) { x = y; return true; } else return false; }
    
	Element &axpy (Element &z, const Element &a, const Element &x, const Element &y) const
		{ return z = a * x + y; }

	Element &addin (Element &x, const Element &y) const
		{ return x += y; }
    
	Element &subin (Element &x, const Element &y) const
		{ return x -= y; }
    
	Element &mulin (Element &x, const Element &y) const
		{ return x *= y; }

	template <class Iterator, class Accessor>
	Element &mulin (Property<Iterator, Accessor> &x, const Element &y) const
		{ return mulin (x.ref (), y); }

	bool divin (Element &x, const Element &y) const
		{ if (x % y == 0) { x /= y; return true; } else return false; }
    
	Element &negin (Element &x) const
		{ return x = - x; }
    
	bool invin (Element &x) const
		{ return x == 1 || x == -1; }

	Element &axpyin (Element &y, const Element &a, const Element &x) const
		{ return y += a * x; }

	std::ostream &write (std::ostream &os) const
		{ return os << "ZZ"; }
    
	std::istream &read (std::istream &is)
		{ return is; }
    
	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }
    
	std::istream &read (std::istream &is, Element &x) const
		{ return is >> x; }

	size_t elementWidth () const
		{ return 10; }

	const Element &zero () const { return _zero; }
	const Element &one () const { return _one; }
	const Element &minusOne () const { return _minus_one; }

private:
	Element _zero, _one, _minus_one;

}; // class Integers

} // namespace LELA

#include "lela/randiter/integers.h"

#include "lela/blas/level1-generic.tcc"
#include "lela/blas/level2-generic.tcc"
#include "lela/blas/level3-generic.tcc"

#endif // __LELA_RING_INTEGERS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
