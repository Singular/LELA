/* lela/randiter/integers.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RANDITER_INTEGERS_H
#define __LELA_RANDITER_INTEGERS_H

#include <sys/time.h>
#include <stdlib.h>

#include "lela/lela-config.h"
#include "lela/ring/integers.h"
#include "lela/integer.h"
#include "lela/randiter/interface.h"
#include "lela/util/property.h"

namespace LELA
{

class IntegerRandIter : public RandIterInterface<integer>
{
    public:
    
	typedef integer Element;
    
	IntegerRandIter (const Integers &__R,
			 const integer &size = 0,
			 const integer &seed = 0)
		: R (__R), _size (size), _seed (seed)
	{
		if (seed == 0)
			_seed = time (NULL);
	}

	IntegerRandIter (const IntegerRandIter &ri)
		: R (ri.R), _size (ri._size), _seed (ri._seed) {}

	~IntegerRandIter() 
	{}
    
	IntegerRandIter &operator = (const IntegerRandIter &ri)
	{
		if (this != &ri) { // guard against self-assignment
			R = ri.R;
			_seed = ri._seed;
			_size = ri._size;
		}
		return *this;
	}
 
	Element &random (Element &a) const
	{
		unsigned int s;
		int value = 0;

		if (_size == 0) {
			s = _seed.get_ui ();
			value = rand_r (&s);
			a = value;
			_seed = s;
		} else {
			unsigned int s;
			int num;

			s = _seed.get_ui ();
			num = rand_r (&s);

			if (_size > 0) {
				unsigned long tmp = _size.get_ui ();
				num %= tmp;
			}

			_seed = s;

			a = num;
		}

		return a;
	}
 
	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }

    private:

	Integers R;

	integer _size;
	mutable integer _seed;
     
}; // class IntegerRandIter
 
} // namespace LELA

#endif // __LELA_RANDITER_INTEGERS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
