/* lela/randiter/rationals.h
 * Copyright 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RANDITER_RATIONALS_H
#define __LELA_RANDITER_RATIONALS_H

#include <sys/time.h>
#include <stdlib.h>

#include "lela/lela-config.h"
#include "lela/ring/rationals.h"
#include "lela/element/rational.h"
#include "lela/randiter/interface.h"
#include "lela/util/property.h"

namespace LELA
{

class RationalRandIter : public RandIterInterface<RationalElement>
{
    public:
    
	typedef RationalElement Element;
    
	RationalRandIter (const Rationals &F,
			  const integer &size = 0,
			  const integer &seed = 0)
		: _F (F), _size (size), _seed (seed)
	{
		if (seed == 0)
			_seed = time (NULL);
	}

	RationalRandIter (const RationalRandIter &R)
		: _F (R._F), _size (R._size), _seed (R._seed) {}

	~RationalRandIter() 
	{}
    
	RationalRandIter &operator = (const RationalRandIter &R)
	{
		if (this != &R) { // guard against self-assignment
			_F = R._F;
			_seed = R._seed;
			_size = R._size;
		}
		return *this;
	}
 
	Element &random (Element &a)  const
	{
		unsigned int s;
		int value = 0;

		if (_size == 0) {
			s = _seed.get_ui ();

			value = rand_r (&s);

			mpz_set_si (mpq_numref (a.rep), value);

			do {
				value = rand_r (&s);
			} while (value == 0);

			const_cast<integer&>(_seed) = s;
			mpz_set_si (mpq_denref (a.rep), value);
		}
		else {
			unsigned int s;
			int num, den;

			s = _seed.get_ui ();
			num = rand_r (&s);

			if (_size > 0) {
				unsigned long tmp = _size.get_ui ();
				num %= tmp;
				den = 1L;
			} else {
				den = rand_r (&s);
			}

			const_cast<integer&>(_seed) = s;

			mpz_set_si (mpq_numref (a.rep), num);
			mpz_set_si (mpq_denref (a.rep), den);
		}

		mpq_canonicalize (a.rep);

		return a;
	}
 
	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }

    private:

	Rationals _F;

	integer _size;
	integer _seed;
     
}; // class RationalRandIter
 
} // namespace LELA

#endif // __LELA_RANDITER_RATIONALS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
