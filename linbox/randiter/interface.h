/* linbox/randiter/archetype.h
 * Copyright 1999-2001 William J Turner,
 *           2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_randiter_archetype_H
#define __LINBOX_randiter_archetype_H

#include "linbox/element/interface.h"

namespace LinBox
{

/** Random ring element generator interface.
 *
 * This abstract base-class defines the interface for obtaining random
 * elements of a ring.
 */
class RandIterInterface
{
public:
    
	/** @name Common Object Interface.
	 * These methods are required of all \Ref{LinBox} ring element generators.
	 */
	//@{
    
	/// element type
	typedef ElementInterface Element;
 
	/** Random ring element creator.
	 * This returns a random ring element from the information supplied
	 * at the creation of the generator.
	 * @return reference to random ring element
	 */
	virtual Element &random (Element &a) const = 0;

	//@} Common Object Iterface
    
}; // class RandIterInterface
 
} // namespace LinBox

#endif // __LINBOX_randiter_archetype_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

