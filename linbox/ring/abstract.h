/* linbox/ring/abstract.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
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

#ifndef __LINBOX_abstract_H
#define __LINBOX_abstract_H

#include <iostream>

#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"
#include "linbox/integer.h"
#include "linbox/linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#endif

namespace LinBox 
{ 

/** \brief ring base class.
 *  \ingroup ring
 *
 * Found in the file \URL{linbox/ring/abstract.h}.
 * Abstract base class used to implement the ring archetype to minimize
 * code bloat.  All public member functions of this class are purely
 * virtual and must be implemented by all derived classes.
 *
 * If a template is instantiated on the ring archetype, we can change the
 * ring it is using by changing the derived class of this class.  This allows
 * us to change the ring used in a template without having to reinstantiate
 * it.  This minimizes code bloat, but it also introduces indirection through
 * the use of pointers and virtual functions which is inefficient.
 */
class RingAbstract
{
public:

	/// element type.
	typedef ElementAbstract Element;

	/// Random iterator generator type.
	typedef RandIterAbstract RandIter;
 
	/** @name Object Management
	 * There are no public constructors for this class.
	 * It should only be used in tandem with \ref{RingArchetype}.
	 */
	//@{

	/** Destructor.
	 * Required because of virtual member functions.
	 * Virtual.
	 */
	virtual ~RingAbstract (void) {}

	/** Virtual copy constructor.
	 * Required because constructors cannot be virtual.
	 * Passes construction on to derived classes.
	 * Purely virtual.
	 * This function is not part of the common object interface.
	 * @return pointer to new object in dynamic memory.
	 */
	virtual RingAbstract* clone () const = 0;

	/** Assignment operator.
	 * Purely virtual.
	 * @return reference to self
	 * @param F constant reference to RingAbstract object
	 */
	virtual RingAbstract &operator= (const RingAbstract &F) = 0;

	/** Initialization of ring element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output ring element x has already been
	 * constructed, but that it is not already initialized.
	 * Purely virtual.
	 * @return reference to ring element.
	 * @param x ring element to contain output (reference returned).
	 * @param y integer.
	 */
	virtual Element &init (Element &x, const integer &y) const = 0;
 
	/** Conversion of ring element to an integer.
	 * This function assumes the output ring element x has already been
	 * constructed, but that it is not already initialized.
	 * Purely virtual.
	 * @return reference to integer.
	 * @param x reference to interger to contain output (reference returned).
	 * @param y constant ring element.
	 */
	virtual integer &convert (integer &x, const Element &y) const = 0;
 
	/** Assignment of one ring element to another.
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &assign (Element &x, const Element &y) const = 0;

	/** Cardinality.
	 * Return integer representing cardinality of the domain.
	 * Returns a non-negative integer for all domains with finite
	 * cardinality, and returns -1 to signify a domain of infinite
	 * cardinality.
	 * Purely virtual.
	 * @return integer representing cardinality of the domain
	 */
	virtual integer &cardinality (integer &c) const = 0;
 
	/** Characteristic.
	 * Return integer representing characteristic of the domain.
	 * Returns a positive integer to all domains with finite characteristic,
	 * and returns 0 to signify a domain of infinite characteristic.
	 * Purely virtual.
	 * @return integer representing characteristic of the domain.
	 */
	virtual integer &characteristic (integer &c) const = 0;

	//@} Object Management

	/** @name Arithmetic Operations
	 * x <- y op z; x <- op y
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized ring elements will
	 * give undefined results.
	 */
	//@{

	/** Equality of two elements.
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return boolean true if equal, false if not.
	 * @param  x ring element
	 * @param  y ring element
	 */
	virtual bool areEqual (const Element &x, const Element &y) const = 0;

	/** Addition.
	 * x = y + z
	 * This function assumes all the ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 * @param  z ring element.
	 */
	virtual Element &add (Element       &x,
			      const Element &y,
			      const Element &z) const = 0;
 
	/** Subtraction.
	 * x = y - z
	 * This function assumes all the ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 * @param  z ring element.
	 */
	virtual Element &sub (Element       &x,
			      const Element &y,
			      const Element &z) const = 0;
 
	/** Multiplication.
	 * x = y * z
	 * This function assumes all the ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 * @param  z ring element.
	 */
	virtual Element &mul (Element       &x,
			      const Element &y,
			      const Element &z) const = 0;
 
	/** Division.
	 * x = y / z
	 * This function assumes all the ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 * @param  z ring element.
	 */
	virtual Element &div (Element &x,
			      const Element &y, const Element &z) const = 0;
 
	/** Additive Inverse (Negation).
	 * x = - y
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &neg (Element &x, const Element &y) const = 0;
 
	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &inv (Element &x, const Element &y) const = 0;

	/** Natural AXPY.
	 * r  = a * x + y
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 * @param  y ring element.
	 */
	virtual Element &axpy (Element       &r, 
			       const Element &a, 
			       const Element &x, 
			       const Element &y) const = 0;

	//@} Arithmetic Operations
 
	/** @name Inplace Arithmetic Operations
	 * x <- x op y; x <- op x
	 */
	//@{

	/** Zero equality.
	 * Test if ring element is equal to zero.
	 * This function assumes the ring element has already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return boolean true if equals zero, false if not.
	 * @param  x ring element.
	 */
	virtual bool isZero (const Element &x) const = 0;
 
	/** One equality.
	 * Test if ring element is equal to one.
	 * This function assumes the ring element has already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return boolean true if equals one, false if not.
	 * @param  x ring element.
	 */
	virtual bool isOne (const Element &x) const = 0;

	/** Inplace Addition.
	 * x += y
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &addin (Element &x, const Element &y) const = 0;
 
	/** Inplace Subtraction.
	 * x -= y
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &subin (Element &x, const Element &y) const = 0;
 
	/** Inplace Multiplication.
	 * x *= y
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &mulin (Element &x, const Element &y) const = 0;

	/** Inplace Division.
	 * x /= y
	 * This function assumes both ring elements have already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &divin (Element &x, const Element &y) const = 0;
 
	/** Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the ring element has already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 */
	virtual Element &negin (Element &x) const = 0;
 
	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the ring elementhas already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 */
	virtual Element &invin (Element &x) const = 0;

	/** Inplace AXPY.
	 * r  += a * x
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * Purely virtual
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 */
	virtual Element &axpyin (Element       &r, 
				 const Element &a, 
				 const Element &x) const = 0;
 
	//@} Inplace Arithmetic Operations

#ifndef __LINBOX_XMLENABLED
	/** @name Input/Output Operations */
	//@{

	/** Print ring.
	 * Purely virtual.
	 * @return output stream to which ring is written.
	 * @param  os  output stream to which ring is written.
	 */
	virtual std::ostream &write (std::ostream &os) const = 0;
 
	/** Read ring.
	 * Purely virtual.
	 * @return input stream from which ring is read.
	 * @param  is  input stream from which ring is read.
	 */
	virtual std::istream &read (std::istream &is) = 0;

	/** Print ring element.
	 * This function assumes the ring element has already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return output stream to which ring element is written.
	 * @param  os  output stream to which ring element is written.
	 * @param  x   ring element.
	 */
	virtual std::ostream &write (std::ostream &os, const Element &x) const = 0;
 
	/** Read ring element.
	 * This function assumes the ring element has already been
	 * constructed and initialized.
	 * Purely virtual.
	 * @return input stream from which ring element is read.
	 * @param  is  input stream from which ring element is read.
	 * @param  x   ring element.
	 */
	virtual std::istream &read (std::istream &is, Element &x) const = 0;

	//@}
#else
	virtual std::ostream &write(std::ostream &os) const = 0;
	virtual bool toTag(Writer &W) const = 0;
	virtual std::ostream &write(std::ostream &os, const Element &e) const = 0;
	virtual bool toTag(Writer &W, const Element &e) const = 0;
	virtual std::istream &read(std::istream &is, Element &e) const = 0;
	virtual bool fromTag(Reader &R, Element &e) const = 0;
#endif

protected:

	/** Default Constructor.
	 * Required by derived classes, but protected because this class should
	 * never be constructed by itself.
	 */
	RingAbstract () {}
 
private:

	/// RingArchetype is friend.
	friend class RingArchetype;

}; // class RingAbstract

} // namespace LinBox

#endif // __LINBOX_abstract_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
