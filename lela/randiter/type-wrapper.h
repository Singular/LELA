/* lela/ring/type-wrapper.h
 * Copyright 1999-2005 William J Turner,
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RANDITER_TYPE_WRAPPER_H
#define __LELA_RANDITER_TYPE_WRAPPER_H

#include <ctime>
#include <vector>

namespace LELA
{

// forward declarations
template <class K> class TypeWrapperRing;
		
/** Unparameterized random ring element generator template.
 * Implements LELA random ring element generator common object interface 
 * for unparameterized rings.
 * Used to generate efficient ring classes for unparameterized rings.
 * Constructs LELA unparameterized random ring element generators from 
 * ring types K.
 * In particular, constructs LELA random ring element generators for
 * unparameterized rings from ring types that
 * adhere to the operations for double, for
 * example TypeWrapperRandIter< float >.
 * Can be used as a pattern to write a particular
 * ring interface, such as, TypeWrapperRandIter< SaclibQ > as
 * a template specialization.
 * This implementation uses the standard C++ random number generator.  Thus,
 * only one random ring element generator can be used at a time since 
 * creating a new one will re-seed the built-in generator and affect all 
 * current LELA generators.
 * @param  K unparameterized ring class
 */
template <class K>
class TypeWrapperRandIter
{
public:
		
	/** @name Common Object Interface.
	 * These methods are required of all LELA random ring element generators.
	 */
	//@{
	 
	/** Ring element type.
	 * The ring element must contain a default constructor, 
	 * a copy constructor, a destructor, and an assignment operator.
	 */
	typedef K Element;    

	/** Constructor from ring, sampling size, and seed.
	 * The random ring element iterator works in the ring F, is seeded
	 * by seed, and it returns any one element with probability no more
	 * than 1/min(size, F.cardinality(c)).
	 * A sampling size of zero means to sample from the entire ring.
	 * A seed of zero means to use some arbitrary seed for the generator.
	 * This implementation sets the sampling size to be no more than the
	 * cardinality of the ring.
	 * @param F LELA ring archetype object in which to do arithmetic
	 * @param size constant integer reference of sample size from which to 
	 *             sample (default = 0)
	 * @param seed constant integer reference from which to seed random number
	 *             generator (default = 0)
	 */
	TypeWrapperRandIter (const TypeWrapperRing<K> &F,
			      const integer &size = 0,
			      const integer &seed = 0)
		: _size(size), _seed(seed)
	{
		if (_seed == integer(0))
			_seed = integer (time (NULL));
			
		integer cardinality;
		F.cardinality (cardinality);

		if ((cardinality != integer (-1)) && (_size > cardinality) )
			_size = cardinality;

#ifdef TRACE
		cout << "created random generator with size " << _size 
		     << " and seed " << _seed << endl;
#endif // TRACE
			
		// Seed random number generator
		srand (_seed.get_si ());
	}

	/** Copy constructor.
	 * Constructs TypeWrapperRandIter object by copying the random ring
	 * element generator.
	 * This is required to allow generator objects to be passed by value
	 * into functions.
	 * In this implementation, this means copying the random ring element
	 * generator to which R._randIter_ptr points.
	 * @param  R TypeWrapperRandIter object.
	 */
	TypeWrapperRandIter(const TypeWrapperRandIter& R)
		: _size(R._size), _seed(R._seed) {}

	/** Destructor.
	 * This destructs the random ring element generator object.
	 * In this implementation, this destroys the generator by deleting 
	 * the random generator object to which _randIter_ptr points.
	 */
	~TypeWrapperRandIter(void) {}
		
	/** Assignment operator.
	 * Assigns TypeWrapperRandIter object R to generator.
	 * In this implementation, this means copying the generator to
	 * which R._randIter_ptr points.
	 * @param  R TypeWrapperRandIter object.
	 */
	TypeWrapperRandIter& operator=(const TypeWrapperRandIter& R)
	{
		if (this != &R) { // guard against self-assignment
			_size = R._size;
			_seed = R._seed;
		}

		return *this;
	}
 
	/** Random ring element creator.
	 * This returns a random ring element from the information supplied
	 * at the creation of the generator.
	 * @return random ring element
	 */
	Element& random (Element& x) const
	{
		// Create new random elements
		if (_size == 0)
			return x = rand();
		else
			return x = static_cast<Element> ((double (rand ()) / RAND_MAX) * _size.get_d ());
	}

	//@} Common Object Iterface
	 
	/** @name Implementation-Specific Methods.
	 * These methods are not required of all 
	 * \Ref{LELA Random ring element generators}
	 * and are included only for this implementation of the archetype.
	 */
	//@{

	/// Default constructor
	TypeWrapperRandIter(void) : _size(0), _seed(0) { time(NULL); }
		
	//@}

private:

	/// Sampling size
	integer _size;
		
	/// Seed
	integer _seed;

}; // template <class K> class TypeWrapperRandIter

} // namespace LELA

#endif // __LELA_RANDITER_TYPE_WRAPPER_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

