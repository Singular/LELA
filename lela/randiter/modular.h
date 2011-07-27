/* lela/randiter/modular.h
 * Copyright 1999-2005 William J Turner,
 *           2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RANDITER_MODULAR_H
#define __LELA_RANDITER_MODULAR_H

#include <iostream>
#include <vector>

#include "lela/integer.h"
#include "lela/ring/modular.h"
#include "lela/util/commentator.h"
#include "lela/randiter/mersenne-twister.h"
#include "lela/util/property.h"

namespace LELA 
{ 

/** Random field base element generator.
 * This encapsulated class is a generator of random field base elements for 
 * the encapsulating field.
 * It is required to contain constructors from a field object and
 * two integers.  The first integer being a cardinality of a set to 
 * draw the random elements from, and the second being a seed for the 
 * random number generator.
 * It is also required to contain a copy constructor, a destructor, and
 * an operator () which acts on a reference to a field base element.  In this 
 * operator (), the random element is placed into the input field base element 
 * and also returned as a reference.
 */
template <class Element>
class ModularRandIter
{
public:

	ModularRandIter (const Modular<Element> &F, 
			 const integer &size = 0, 
			 const integer &seed = 0)
		: _MT (seed.get_ui ()), _F (F), _size (size), _seed (seed.get_ui ())
	{
		integer cardinality;

		F.cardinality (cardinality);

		if ((_size == 0) || (_size > cardinality))
			_size = cardinality;

		commentator.report (10, INTERNAL_DESCRIPTION)
			<< "Created random generator with size " << _size 
			<< " and seed " << _seed << std::endl;
	}

	ModularRandIter (const ModularRandIter<Element> &R) 
		: _F (R._F), _size (R._size), _seed (R._seed) {}

	~ModularRandIter () {}
    
	ModularRandIter<Element> &operator = (const ModularRandIter<Element> &R)
	{
		if (this != &R) { // guard against self-assignment
			_size = R._size;
			_seed = R._seed;
			_MT.setSeed (_seed);
		}

		return *this;
	}

	Element &random (Element &a) const
		{ return _F.init (a, _MT.randomIntRange (0, _size.get_ui ())); }

	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }

private:
	MersenneTwister _MT;

	/// Field in which arithmetic is done
	Modular<Element> _F;

	/// Sampling size
	integer _size;
    
	/// Seed
	long _seed;

}; // class ModularRandIter

template <class Element>
class Modular<Element>::RandIter
{
	ModularRandIter<Element> _r;

public:
	RandIter (const Modular<Element> &F, const integer &size = 0, const integer &seed = 0)
		: _r (F, size, seed) {}
	RandIter (const Modular<Element>::RandIter &r)
		: _r (r._r) {}

	~RandIter () {}
	RandIter &operator= (const RandIter &r)
		{ _r = r._r; return *this; }

	Element &random (Element &a) const
		{ return _r.random (a); }

	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }
};

template <>
class Modular<uint16>::RandIter
{
	MersenneTwister _r;
	uint16 _size;
	uint16 _seed;

public:
	typedef uint16 Element;

	RandIter (const Modular<Element> &F, const integer &size = 0, const integer &seed = 0)
		: _r (seed.get_ui ()), _size (size.get_ui ()), _seed (seed.get_ui ())
	{
		integer c;

		F.cardinality (c);

		lela_check (c != -1);

		if ((_size == 0) || (_size > c.get_d ()))
			_size = c.get_ui ();
	}

	RandIter (const Modular<Element>::RandIter &r)
		: _r (r._r), _size (r._size), _seed (r._seed) {}

	~RandIter () {}
	RandIter &operator= (const RandIter &r)
		{ _r = r._r; return *this; }

	Element &random (Element &a) const
		{ return a = _r.randomIntRange (0, _size); }

	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }
};

template <>
class Modular<uint32>::RandIter
{
	MersenneTwister _r;
	uint32 _size;
	uint32 _seed;

public:
	typedef uint32 Element;

	RandIter (const Modular<Element> &F, const integer &size = 0, const integer &seed = 0)
		: _r (seed.get_ui ()), _size (size.get_ui ()), _seed (seed.get_ui ())
	{
		integer c;

		F.cardinality (c);

		lela_check (c != -1);

		if ((_size == 0) || (_size > c.get_d ()))
			_size = c.get_ui ();
	}

	RandIter (const Modular<Element>::RandIter &r)
		: _r (r._r), _size (r._size), _seed (r._seed) {}

	~RandIter () {}

	RandIter &operator= (const RandIter &r)
		{ _r = r._r; return *this; }

	Element &random (Element &a) const
		{ return a = _r.randomIntRange (0, _size); }

	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }
};

} // namespace LELA 

#endif // __LELA_RANDITER_MODULAR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
