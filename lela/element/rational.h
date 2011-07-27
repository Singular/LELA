/* lela/element/rational.h
 * Copyright 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ----------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ELEMENT_RATIONAL_H
#define __LELA_ELEMENT_RATIONAL_H

#include "lela/integer.h"

#include <gmp.h>

namespace LELA
{

// Forward declarations
class Rationals;
class RationalRandIter;

/** \brief Element of the field Rationals
\ingroup element
 */
class RationalElement
{
    public:

	RationalElement () { mpq_init (rep); }

	RationalElement (const RationalElement &a) 
		{ mpq_init (rep); mpq_set (rep, a.rep); }

	~RationalElement () { mpq_clear (rep); }

	RationalElement &operator = (const RationalElement& a)
	{
		if (this != &a) // guard against self-assignment
			mpq_set (rep, a.rep);

		return *this;
	}

	/** Constructor.
	 * Constructs ring element from an mpq_t
	 * Not part of the interface.
	 * Creates new copy of element object in dynamic memory.
	 * @param  elem_ptr  pointer to \ref{ElementAbstract}
	 */
	RationalElement (mpq_t _rep)
	{
		mpq_init (rep);
		mpq_set (rep, _rep);
	}

	/** Constructor
	 * Initialize from numerator and denominator
	 */
	RationalElement (const integer &num, const integer &den) 
	{
		mpq_init (rep);
		mpz_set (mpq_numref (rep), num.get_mpz_t ());
		mpz_set (mpq_denref (rep), den.get_mpz_t ());
	}
	
	// Added by Rich Seagraves to take care of some headaches
	/** Constructor
	 *  Initalizes from a single integer, (which is assumed to be the
	 *  numerator, with the denominator being 1)
	 */
	RationalElement (const integer &num)
	{
		mpq_init (rep);
		mpq_set_z(rep, num.get_mpz_t ());
	}

	mpq_ptr get_rep () { return rep; }

	//@}p
    
    private:

	friend class Rationals;
	friend class RationalRandIter;

	mutable mpq_t rep;
};

} // namespace LELA

#endif // __LELA_ELEMENT_RATIONAL_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

