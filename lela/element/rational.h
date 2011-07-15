/* lela/element/rational.h
 * Copyright 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ----------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LELA_ELEMENT_RATIONAL_H
#define __LELA_ELEMENT_RATIONAL_H

#include "lela/integer.h"
#include "lela/element/interface.h"

#include <gmp.h>

namespace LELA
{

// Forward declarations
class Rationals;
class RationalRandIter;

/** \brief elements of GMP_Rationals.
\ingroup element
 */
class RationalElement
{
    public:

	/** @name Common Object Interface for LELA Ring elements.
	 * These methods are required of all \ref{LELA} 
	 * {@link Rings ring} elements.
	 */
	//@{

	/** Default constructor.
	 * This constructor is required to allow 
	 * {@link Rings ring} elements to be primitive C++ types.
	 * Because constructor does not know what {@link Rings ring} 
	 * the element belongs to, it cannot actually construct the element.
	 * In this implementation, the constructor it sets _elem_ptr
	 * to the null pointer.  Initialization of the element is done through
	 * the ring function init where the ring is known.
	 */
	RationalElement () { mpq_init (rep); }

	/** Copy constructor.
	 * This constructor is required to allow 
	 * {@link Rings ring} elements to be primitive C++ types, 
	 * and to allow ring elements to be passed by value into 
	 * functions.
	 * Constructs {@link Rings ring} element by copying the 
	 * {@link Rings ring} element.
	 * In this implementation, this means copying the element to
	 * which a._elem_ptr points.
	 * @param  a ring element.
	 */
	RationalElement (const RationalElement &a) 
		{ mpq_init (rep); mpq_set (rep, a.rep); }

	/** Destructor.
	 * In this implementation, this destroys element by deleting ring 
	 * element to which _elem_ptr points.
	 */
	~RationalElement () { mpq_clear (rep); }

	/** Assignment operator.
	 * Assigns element a to element.  
	 * In this implementation, this is done 
	 * by copying ring element to which _elem_ptr points.
	 * @param  a ring element.
	 */
	RationalElement &operator = (const RationalElement& a)
	{
		if (this != &a) // guard against self-assignment
			mpq_set (rep, a.rep);

		return *this;
	}

	//@} Common Object Interface

	/** @name Implementation-Specific Methods.
	 * These methods are not required of all LELA ring elements
	 * and are included only for this implementation of the archetype.
	 */
	//@{

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
		/*
		//potential error, error occurs when |num| is bigger than largest int
		mpz_set_si (mpq_numref(rep), num);
		mpz_set_si (mpq_denref(rep), integer(1));
		*/
		mpq_set_z(rep, num.get_mpz_t ());
	}

	mpq_ptr get_rep() { return rep; }

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

