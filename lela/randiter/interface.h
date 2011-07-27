/* lela/randiter/interface.h
 * Copyright 1999-2001 William J Turner,
 *           2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RANDITER_INTERFACE_H
#define __LELA_RANDITER_INTERFACE_H

#include "lela/util/property.h"

namespace LELA
{

/** Random ring element generator interface.
 *
 * This abstract base-class defines the interface for obtaining random
 * elements of a ring.
 */
template <class _Element>
class RandIterInterface
{
public:
    
	/** @name Common Object Interface.
	 * These methods are required of all \Ref{LELA} ring element generators.
	 */
	//@{
    
	/// element type
	typedef _Element Element;
 
	/** Random ring element creator.
	 * This returns a random ring element from the information supplied
	 * at the creation of the generator.
	 * @return reference to random ring element
	 */
	virtual Element &random (Element &a) const = 0;

	/** Version of random taking a property as input */
	template <class Iterator, class Accessor>
	Element &random (Property<Iterator, Accessor> a) const
		{ return random (a.ref ()); }

	//@} Common Object Iterface
    
}; // class RandIterInterface
 
} // namespace LELA

#endif // __LELA_randiter_archetype_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

