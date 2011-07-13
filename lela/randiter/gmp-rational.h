/* lela/randiter/gmp-rational.h
 * Copyright 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LELA_RANDITER_GMP_RATIONAL_H
#define __LELA_RANDITER_GMP_RATIONAL_H

#include "lela/lela-config.h"
#include "lela/ring/gmp-rational.h"
#include "lela/element/gmp-rational.h"

#include <sys/time.h>
#include <stdlib.h>

namespace LELA
{

class GMPRationalRandIter
{
    public:
    
	typedef GMPRationalElement Element;
    
	GMPRationalRandIter (const GMPRationalField &F,
			     const integer &size = 0,
			     const integer &seed = 0)
		: _F (F), _size (size), _seed (seed)
	{
		if (seed == 0)
			_seed = time (NULL);
	}

	GMPRationalRandIter (const GMPRationalRandIter &R)
		: _F (R._F), _size (R._size), _seed (R._seed) {}

	~GMPRationalRandIter() 
	{}
    
	GMPRationalRandIter &operator = (const GMPRationalRandIter &R)
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
 
    private:

	GMPRationalField _F;

	integer _size;
	integer _seed;
     
}; // class GMPRationalRandIter
 
} // namespace LELA

#endif // __LELA_RANDITER_GMP_RANDOM_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
