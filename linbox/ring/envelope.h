/* linbox/ring/envelope.h
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
 *
 * ------------------------------------
 * 2002-05-14 William J. Turner <wjturner@acm.org>
 * 
 * changed randIter to RandIter.
 * ------------------------------------
 */

#ifndef __LINBOX_ring_envelope_H
#define __LINBOX_ring_envelope_H

#include <iostream>

#include "linbox/integer.h"
#include "linbox/element/envelope.h"
#include "linbox/ring/abstract.h"
#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"
#include "linbox/randiter/envelope.h"

#include "linbox/linbox-config.h"
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#endif //__LINBOX_XMLENABLED

// Namespace in which all LinBox code resides
namespace LinBox 
{

// Forward declarations
template <class Ring> class RandIterEnvelope;

/** \brief Derived class used to implement the ring archetype
 *  \ingroup ring
 *
 * Helps to minimize
 * code bloat.  This class implements all purely virtual member functions
 * of the abstract base class.  This class is used to wrap a
 * LinBox
 * ring so that it might be used with the Ring archetype.
 */
template <class Ring>
class RingEnvelope : public RingAbstract
{
public:

	/** element type.
	 * It is derived from the class ElementAbstract, and it must contain
	 * a wrapped ring element.
	 */
	typedef ElementEnvelope<Ring> Element;

	/** Random iterator generator type.
	 * It is derived from the class RandIterAbstract, and it must contain
	 * a wrapped ring random iterator generator.
	 */
	typedef RandIterEnvelope<Ring> RandIter;

	/** @name Object Management
	 */
	//@{
 
	/** Default constructor.
	 * In this implementation, this means copying the ring {\tt E.\_ring}.
	 */
	RingEnvelope (void) {}

	/** Constructor from ring to be wrapped.
	 * @param F Ring object to be wrapped.
	 */
	RingEnvelope (const Ring &F) : _ring (F) {}
 
	/** Copy constructor.
	 * Constructs RingEnvelope object by copying the ring.
	 * This is required to allow ring objects to be passed by value
	 * into functions.
	 * In this implementation, this means copying the ring {\tt E.\_ring}.
	 * @param  E RingEnvelope object.
	 */
	RingEnvelope (const RingEnvelope &E) : _ring (E._ring) {}

#ifdef __LINBOX_XMLENABLED
RingEnvelope(Reader &R) : _ring(R) {}
#endif

 
	/** Virtual copy constructor.
	 * Required because constructors cannot be virtual.
	 * Passes construction on to derived classes.
	 * This function is not part of the common object interface.
	 * @return pointer to new object in dynamic memory.
	 */
	RingAbstract* clone () const
		{ return new RingEnvelope (*this); }

	/** Assignment operator.
	 * Required by abstract base class.
	 * @return reference to RingAbstract object for self
	 * @param F constant reference to RingAbstract object
	 */
	RingAbstract &operator= (const RingAbstract &F)
	{
		if (this != &F) // guard against self-assignment
			_ring = static_cast<const RingEnvelope&> (F)._ring;

		return *this;
	}

	/** Initialization of ring base element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output ring base element x has already been
	 * constructed, but that it is not already initialized.
	 * This is not a specialization of the template function because
	 * such a specialization is not allowed inside the class declaration.
	 * @return reference to ring base element.
	 * @param x ring base element to contain output (reference returned).
	 * @param y integer.
	 */
	ElementAbstract &init (ElementAbstract &x, const integer &y = 0) const
	{
		_ring.init (static_cast<ElementEnvelope<Ring>&> (x)._elem, y);
		return x;
	}
 
	/** Conversion of ring base element to a template class T.
	 * This function assumes the output ring base element x has already been
	 * constructed, but that it is not already initialized.
	 * @return reference to template class T.
	 * @param x template class T to contain output (reference returned).
	 * @param y constant ring base element.
	 */
	integer &convert (integer &x, const ElementAbstract &y) const
	{
		_ring.convert (x, static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}
 
	/** Assignment of one ring base element to another.
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &assign (ElementAbstract &x, const ElementAbstract &y) const
	{
		_ring.assign (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			       static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}

	/** Cardinality.
	 * Return integer representing cardinality of the domain.
	 * Returns a non-negative integer for all domains with finite
	 * cardinality, and returns -1 to signify a domain of infinite
	 * cardinality.
	 * @return integer representing cardinality of the domain
	 */
	integer &cardinality (integer &c) const
		{ return _ring.cardinality (c); }
 
	/** Characteristic.
	 * Return integer representing characteristic of the domain.
	 * Returns a positive integer to all domains with finite characteristic,
	 * and returns 0 to signify a domain of infinite characteristic.
	 * @return integer representing characteristic of the domain.
	 */
	integer &characteristic (integer &c) const
		{ return _ring.characteristic (c); }

	//@} Object Management

	/** @name Arithmetic Operations
	 * x <- y op z; x <- op y
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized ring base elements will
	 * give undefined results.
	 */
	//@{

	/** Equality of two elements.
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return boolean true if equal, false if not.
	 * @param  x ring base element
	 * @param  y ring base element
	 */
	bool areEqual (const ElementAbstract &x, const ElementAbstract &y) const
	{
		return _ring.areEqual (static_cast<const ElementEnvelope<Ring>&> (x)._elem,
					static_cast<const ElementEnvelope<Ring>&> (y)._elem);
	}

	/** Addition.
	 * x = y + z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	ElementAbstract &add (ElementAbstract &x,
			      const ElementAbstract &y,
			      const ElementAbstract &z) const
	{
		_ring.add (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (y)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (z)._elem);
		return x;
	}
 
	/** Subtraction.
	 * x = y - z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	ElementAbstract &sub (ElementAbstract &x,
			      const ElementAbstract &y,
			      const ElementAbstract &z) const
	{
		_ring.sub (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (y)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (z)._elem);
		return x;
	}
 
	/** Multiplication.
	 * x = y * z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	ElementAbstract &mul (ElementAbstract &x,
			      const ElementAbstract &y,
			      const ElementAbstract &z) const
	{
		_ring.mul (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (y)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (z)._elem);
		return x;
	}
 
	/** Division.
	 * x = y / z
	 * This function assumes all the ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 * @param  z ring base element.
	 */
	ElementAbstract &div (ElementAbstract &x,
			      const ElementAbstract &y,
			      const ElementAbstract &z) const
	{
		_ring.div (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (y)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (z)._elem);
		return x;
	}
 
	/** Additive Inverse (Negation).
	 * x = - y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &neg (ElementAbstract &x, const ElementAbstract &y) const
	{
		_ring.neg (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}
 
	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &inv (ElementAbstract &x, const ElementAbstract &y) const
	{
		_ring.inv (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			    static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}

	/** Natural AXPY.
	 * r  = a * x + y
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 * @param  y ring element.
	 */
	ElementAbstract &axpy (ElementAbstract &r, 
			       const ElementAbstract &a, 
			       const ElementAbstract &x, 
			       const ElementAbstract &y) const
	{
		_ring.axpy (static_cast<ElementEnvelope<Ring>&> (r)._elem,
			     static_cast<const ElementEnvelope<Ring>&> (a)._elem,
			     static_cast<const ElementEnvelope<Ring>&> (x)._elem,
			     static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return r;
	}
 
	//@} Arithmetic Operations
 
	/** @name Inplace Arithmetic Operations
	 * x <- x op y; x <- op x
	 */
	//@{

	/** Zero equality.
	 * Test if ring base element is equal to zero.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals zero, false if not.
	 * @param  x ring base element.
	 */
	bool isZero (const ElementAbstract &x) const
	{ return _ring.isZero (static_cast<const ElementEnvelope<Ring>&> (x)._elem); }
 
	/** One equality.
	 * Test if ring base element is equal to one.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals one, false if not.
	 * @param  x ring base element.
	 */
	bool isOne (const ElementAbstract &x) const
		{ return _ring.isOne (static_cast<const ElementEnvelope<Ring>&> (x)._elem); }

	/** Inplace Addition.
	 * x += y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &addin (ElementAbstract &x, const ElementAbstract &y) const
	{
		_ring.addin (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			      static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}
 
	/** Inplace Subtraction.
	 * x -= y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &subin (ElementAbstract &x, const ElementAbstract &y) const
	{
		_ring.subin (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			      static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}
 
	/** Inplace Multiplication.
	 * x *= y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &mulin (ElementAbstract &x, const ElementAbstract &y) const
	{
		_ring.mulin (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			      static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}

	/** Inplace Division.
	 * x /= y
	 * This function assumes both ring base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 * @param  y ring base element.
	 */
	ElementAbstract &divin (ElementAbstract &x, 
				const ElementAbstract &y) const
	{
		_ring.divin (static_cast<ElementEnvelope<Ring>&> (x)._elem,
			      static_cast<const ElementEnvelope<Ring>&> (y)._elem);
		return x;
	}
 
	/** Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 */
	ElementAbstract &negin (ElementAbstract &x) const
	{
		_ring.negin (static_cast<ElementEnvelope<Ring>&> (x)._elem);
		return x;
	}
 
	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the ring base elementhas already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x ring base element (reference returned).
	 */
	ElementAbstract &invin (ElementAbstract &x) const
	{
		_ring.invin (static_cast<ElementEnvelope<Ring>&> (x)._elem);
		return x;
	}

	/** Inplace AXPY.
	 * r  += a * x
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 */
	ElementAbstract &axpyin (ElementAbstract &r, 
				 const ElementAbstract &a, 
				 const ElementAbstract &x) const
	{
		_ring.axpyin (static_cast<ElementEnvelope<Ring>&> (r)._elem,
			       static_cast<const ElementEnvelope<Ring>&> (a)._elem,
			       static_cast<const ElementEnvelope<Ring>&> (x)._elem);
		return r;
	}
 
	//@} Inplace Arithmetic Operations

#ifndef __LINBOX_XMLENABLED
	/** @name Input/Output Operations */
	//@{

	/** Print ring.
	 * @return output stream to which ring is written.
	 * @param  os  output stream to which ring is written.
	 */
	std::ostream &write (std::ostream &os) const
		{ return _ring.write (os); }
 
	/** Read ring.
	 * @return input stream from which ring is read.
	 * @param  is  input stream from which ring is read.
	 */
	std::istream &read (std::istream &is)
		{ return _ring.read (is); }

	/** Print ring base element.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return output stream to which ring base element is written.
	 * @param  os  output stream to which ring base element is written.
	 * @param  x   ring base element.
	 */
	std::ostream &write (std::ostream &os, const ElementAbstract &x) const
		{ return _ring.write (os, static_cast<const ElementEnvelope<Ring>&> (x)._elem); }
 
	/** Read ring base element.
	 * This function assumes the ring base element has already been
	 * constructed and initialized.
	 * @return input stream from which ring base element is read.
	 * @param  is  input stream from which ring base element is read.
	 * @param  x   ring base element.
	 */
	std::istream &read (std::istream &is, ElementAbstract &x) const
		{ return _ring.read (is, static_cast<ElementEnvelope<Ring>&> (x)._elem); }

	//@}
#else
	std::ostream &write (ostream &os) const
		{ return _ring.write(os); }

	bool toTag (Writer &W) const
		{ return _ring.toTag (W); }

	std::ostream &write(ostream &os, const ElementAbstract &x) const
		{ return _ring.write (os, static_cast<const ElementEnvelope<Ring>&>(x)._elem); }

	bool toTag(Writer &W, const ElementAbstract &x) const
		{ return _ring.toTag (W, static_cast<const ElementEnvelope<Ring>&>(x)._elem); }

	std::istream &read(istream &is, ElementAbstract &x) const
		{ return _ring.read (is, static_cast<ElementEnvelope<Ring>&>(x)._elem); }

	bool fromTag(Reader &R, ElementAbstract &x) const
		{ return _ring.fromTag (R, static_cast<ElementEnvelope<Ring>&>(x)._elem); }
#endif

protected:

	friend class RandIterEnvelope<Ring>;

	/// Wrapped ring.
	Ring _ring;

}; // class RingEnvelope

} // namespace LinBox

#include "linbox/randiter/envelope.h"

#endif // __LINBOX_ring_envelope_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

