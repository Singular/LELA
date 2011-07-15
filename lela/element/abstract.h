/* lela/element/abstract.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LELA_ELEMENT_ABSTRACT_H
#define __LELA_ELEMENT_ABSTRACT_H

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
 * arbitrary data and dispose, which is a pointer to a function which
 * is invoked by the destructor. An implementation of RingInterface
 * can set ref to point to element-data.
 *
 * A ring-implementation which uses this class must set the pointer
 * dispose whenever the element is initialised.
 */
class AbstractElement
{
public:
	/// Pointer to a function which disposes of the ring element
	void (*dispose) (AbstractElement &);

	/// Pointer to arbitrary data
	void *ref;

	AbstractElement ()
		: dispose (NULL), ref (NULL) {}

	~AbstractElement () { if (dispose != NULL) dispose (*this); }
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
