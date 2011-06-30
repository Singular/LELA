/* linbox/element/interface.h
 * Copyright 1999-2001 William J Turner,
 *           2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ELEMENT_INTERFACE_H
#define __LINBOX_ELEMENT_INTERFACE_H

namespace LinBox 
{ 

/** \brief Element base class
 *
 * This is the element-type for \ref RingInterface. Element-types for
 * rings inheriting from RingInterface should inherit from this type.
 */
class ElementInterface {};

} // namespace LinBox

#endif // __LINBOX_ELEMENT_INTERFACE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
