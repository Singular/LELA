/* lela/randiter/mersenne-twister.h
 * Copyright 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * Header file for implementation of the Mersenne twister pseudo-random number
 * generator, crated by Makoto Matsumoto and Takuji Nishimura. Further
 * information on the Mersenne Twister can be found at
 * http://en.wikipedia.org/wiki/Mersenne_Twister
 *
 * This forms the basic underlying algorithm for most psuedo-random number
 * generation in LELA.
 *
 * N.B. This module is tested in test-modular under the random number
 * generator tests.
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_mersenne_twister_H
#define __LELA_mersenne_twister_H

#include <vector>

#include "lela/integer.h"

namespace LELA
{

class MersenneTwister 
{
    public:
	MersenneTwister (uint32 seed = 0);

	uint32 reload ();
	uint32 randomInt ();
	uint32 randomInt () const
		{ return const_cast<MersenneTwister&>(*this).randomInt();}

	uint64 randomLongLong ();
	uint64 randomLongLong () const
		{ return const_cast<MersenneTwister&>(*this).randomLongLong();}

	uint32 randomIntRange (uint32 start, uint32 end);
	uint32 randomIntRange (uint32 start, uint32 end) const
		{ return const_cast<MersenneTwister&>(*this).randomIntRange(start,end); }

	double randomDouble ();
	double randomDouble ()  const
		{ return const_cast<MersenneTwister&>(*this).randomDouble(); }

	double randomDoubleRange (double start, double end)
		{ return randomDouble () * (end - start) + start; }
	double randomDoubleRange (double start, double end) const
		{ return const_cast<MersenneTwister&>(*this).randomDoubleRange(start,end); }

	void setSeed (uint32 seed);

	static int getSeed ();

    private:
	std::vector<uint32>           _state;
	std::vector<uint32>::iterator _next;
	int                           _left;
};
 
} // namespace LELA

#endif // __LELA_mersenne_twister_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

