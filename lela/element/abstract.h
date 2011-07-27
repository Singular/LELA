/* lela/element/abstract.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_ELEMENT_ABSTRACT_H
#define __LELA_ELEMENT_ABSTRACT_H

#include "lela/ring/interface.h"

namespace LELA 
{ 

/** Abstract element
 *
 * This class can be used as a generic element-type against which
 * RingInterface may be instantiated. At a cost of some overhead, this
 * permits that LELA-functions are instantiated only once when
 * multiple rings with multiple element-types are used.
 *
 * This class contains two public members: ref, which is a pointer to
 * arbitrary data and ring, which is a pointer to the element's ring
 * from which ref and unref are invoked. An implementation of
 * RingInterface can set ref to point to element-data.
 *
 * A ring-implementation which uses this class must set the pointer
 * dispose whenever the element is initialised.
 *
 * \ingroup ring
 */
class AbstractElement
{
public:
	/// Pointer back to the ring
	RingInterface *ring;

	/// Pointer to arbitrary data
	void *ref;

	AbstractElement ()
		: ring (NULL), ref (NULL) {}

	AbstractElement (const AbstractElement &e)
		: ring (e.ring), ref (e.ref)
		{ if (ring != NULL) ring->ref (*this); }

	~AbstractElement () { if (ring != NULL) ring->unref (*this); }

	AbstractElement &operator = (const AbstractElement &e)
	{
		ring = e.ring;
		ref = e.ref;
		if (ring != NULL) ring->ref (*this);
	}
};

} // namespace LELA

#endif // __LELA_ELEMENT_ABSTRACT_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
