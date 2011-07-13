/* lela/randiter/nonzero.h
 * Copyright 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LELA_randiter_nonzero_H
#define __LELA_randiter_nonzero_H

#include "lela/lela-config.h"
#include "lela/ring/interface.h"
#include "lela/randiter/interface.h"

namespace LELA
{

/** Random iterator for nonzero random numbers
 *
 * Wraps around an existing random iterator and ensures that the output
 * is entirely nonzero elements.
 **/
template <class Ring, class RandIter = typename Ring::RandIter>
class NonzeroRandIter
{
public:
    
	typedef typename Ring::Element Element;

	NonzeroRandIter (const Ring &F, const RandIter &r)
		: _F (F), _r (r) {}

	NonzeroRandIter (const NonzeroRandIter &R)
		: _F (R._F), _r (R._r) {}

	~NonzeroRandIter() 
		{}

	NonzeroRandIter &operator = (const NonzeroRandIter &R)
	{
		if (this != &R) { // guard against self-assignment
			_F = R._F;
			_r = R._r;
		}

		return *this;
	}

	typename Ring::Element &random (typename Ring::Element &a)  const
	{
		do _r.random (a); while (_F.isZero (a));
		return a;
	}

private:

	Ring     _F;
	RandIter _r;
     
}; // class NonzeroRandIter
 
} // namespace LELA

#endif // __LELA_randiter_nonzero_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

