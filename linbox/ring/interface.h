/* linbox/ring/archetype.h
 * Copyright 1999-2005 William J Turner,
 *           2001 Bradford Hovinen
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */


#ifndef __LINBOX_RING_INTERFACE_H
#define __LINBOX_RING_INTERFACE_H

#include <iostream>

#include "linbox/linbox-config.h"
#include "linbox/element/interface.h"
#include "linbox/randiter/interface.h"
#include "linbox/integer.h"
#include "linbox/util/property.h"
#include "linbox/util/error.h"

namespace LinBox
{

/** Abstract ring
 *
 * This class defines the ring-interface. It is an abstract base-class
 * from which rings may be derived, though this is not required
 * provided that a ring satisfy this interface.
 */
class RingInterface
{
public:

	/** @name Common Object Interface for a LinBox Ring.
	 * These methods are required of all \ref{LinBox} rings.
	 */
	//@{
    
	/// the type in which ring elements are represented.
	typedef ElementInterface Element;

	/// An object of this type is a generator of random ring elements.
	typedef RandIterInterface RandIter;
    
	/// @name Object Management
	//@{
    
	/** \brief Initialization of ring element from an integer.
	 *
	 * x becomes the image of n under the natural map from the integers
	 * to the prime subring.  It is the result obtained from adding n 1's 
	 * in the ring.
	 *
	 * This function assumes the output ring element x
	 * has already been constructed, but that it is not
	 * necessarily already initialized. In this
	 * archetype implementation, this means the <tt> _elem_ptr</tt> of
	 * x exists, but that it may be the null pointer.
	 *
	 * @return reference to x.
	 * @param x output ring element.
	 * @param n input integer.
	 */
	virtual Element &init (Element &x, const integer &n = 0) const = 0;
  
	/** \brief Conversion of ring element to an integer.  The
	 * meaning of conversion is specific to each ring
	 * class. However, if x is in the image of the integers in the
	 * ring, the integer n returned is such that an init from n
	 * will reproduce x. Most often 0 &leq; n &lt; characteristic.
	 *
	 * @return reference to n.
	 * @param n output integer.
	 * @param x input ring element.
	 */
	virtual integer &convert (integer &n, const Element &y) const = 0;
    
	/** \brief  Assignment of one ring element to another.
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x
	 * @param  x destination ring element.
	 * @param  y source ring element.
	 */
	virtual Element &assign (Element &x, const Element &y) const = 0;
    
	/** \brief Cardinality.
	 *
	 * Return c, integer representing cardinality of the ring.  c
	 * becomes a positive integer for all rings with finite
	 * cardinality, and 0 to signify a ring of infinite
	 * cardinality.
	 */
	virtual integer &cardinality (integer &c) const = 0;
    
	/** \brief Characteristic.
	 *
	 * Return c, integer representing characteristic of the ring
	 * (the least positive n such that the sum of n copies of x is
	 * 0 for all ring elements x).  c becomes a positive integer
	 * for all rings with finite characteristic, and 0 to signify
	 * a ring of infinite characteristic.
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
    
	/** \brief Equality of two elements.
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return boolean true if equal, false if not.
	 * @param  x ring element
	 * @param  y ring element
	 */
	virtual bool areEqual (const Element &x, const Element &y) const = 0;
    
	/** \brief Addition, x <-- y + z.
	 *
	 * This function assumes all the ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 */
	virtual Element &add (Element &x, const Element &y, const Element &z) const = 0;
    
	/** \brief Subtraction, x <-- y - z.
	 *
	 * This function assumes all the ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 */
	virtual Element &sub (Element &x, const Element &y, const Element &z) const = 0;
    
	/** \brief Multiplication, x <-- y * z.
	 *
	 * This function assumes all the ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 */
	virtual Element &mul (Element &x, const Element &y, const Element &z) const = 0;
    
	/** Division, x <-- y / z.
	 *
	 * This function attempts to compute the quotient of y and z
	 * and, if it exists, stores the result in x. If the quotient
	 * does not exist in the ring, the value of x is left
	 * unchangd.
	 *
	 * This function assumes all the ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return true if the quotient exists in the ring, false if not
	 */
	virtual bool div (Element &x, const Element &y, const Element &z) const = 0;
    
	/** \brief Additive Inverse (Negation), x <-- - y.
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 */
	virtual Element &neg (Element &x, const Element &y) const = 0;
    
	/** \brief Multiplicative Inverse, x <-- 1 / y.
	 *
	 * This computes the multiplicative inverse if possible and if
	 * it exists, stores it in x. If not, then x is left
	 * unchanged.
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return true if the inverse exists in the ring, false if not
	 */
	virtual bool inv (Element &x, const Element &y) const = 0;
    
	/** \brief Ring element AXPY, r  <-- a * x + y.
	 *
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to r.
	 */
	virtual Element &axpy (Element &r, const Element &a, const Element &x, const Element &y) const = 0;

	//@} Arithmetic Operations
    
	/** @name Predicates
	 */
	//@{
	/** Equality with the zero-element
	 *
	 * Test if ring element is equal to zero.
	 * This function assumes the ring element has already been 
	 * constructed and initialized.
	 *
	 * @return boolean true if equals zero, false if not.
	 * @param  x ring element.
	 */
	virtual bool isZero (const Element &x) const = 0;
    
	/** Equality with the one-element
	 *
	 * Test if ring element is equal to one.
	 * This function assumes the ring element has already been 
	 * constructed and initialized.
	 *
	 * @return boolean true if equals one, false if not.
	 * @param  x ring element.
	 */
	virtual bool isOne (const Element &x) const = 0;
	//@}

	/** @name Inplace Arithmetic Operations 
	 * x <- x op y; x <- op x
	 *
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized ring elements will
	 * give undefined results.
	 */
	//@{
    
	/** Inplace Addition; x += y
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &addin (Element &x, const Element &y) const = 0;
    
	/** Inplace Subtraction; x -= y
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &subin (Element &x, const Element &y) const = 0;
 
	/** Inplace Multiplication; x *= y
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual Element &mulin (Element &x, const Element &y) const = 0;
   
	/** Inplace Multiplication using a property; x *= y
	 *
	 * This function does the same thing as above but takes a
	 * reference to a Property rather than to an element. This is
	 * required for that the scal operation on sparse matrices
	 * work.
	 *
	 * @return reference to x.
	 * @param  x Property (reference returned).
	 * @param  y ring element.
	 */
	template <class Iterator>
	Property<Iterator> &mulin (Property<Iterator> &x, const Element &y) const
	{
		Element t = x;
		mulin (t, y);
		x = t;
	}
   
	/** Inplace Division; x /= y
	 *
	 * This function attempts to compute the quotient of x and y
	 * and, if it exists, stores the result in x. If the quotient
	 * does not exist in the ring, the value of x is left
	 * unchanged.
	 *
	 * This function assumes both ring elements have already been 
	 * constructed and initialized.
	 *
	 * @return true if the quotient exists in the ring, false otherwise
	 * @param  x ring element (reference returned).
	 * @param  y ring element.
	 */
	virtual bool divin (Element &x, const Element &y) const = 0;
    
	/** Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the ring element has already been 
	 * constructed and initialized.
	 *
	 * @return reference to x.
	 * @param  x ring element (reference returned).
	 */
	virtual Element &negin (Element &x) const = 0;
    
	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the ring element has already been 
	 * constructed and initialized.
	 *
	 * @return true if the inverse of x exists in the ring, false if not
	 * @param  x ring element.
	 */
	virtual bool invin (Element &x) const = 0;
    
	/** Inplace AXPY.
	 * r  += a * x
	 * This function assumes all ring elements have already been 
	 * constructed and initialized.
	 * @return reference to r.
	 * @param  r ring element (reference returned).
	 * @param  a ring element.
	 * @param  x ring element.
	 */
	virtual Element &axpyin (Element &r, const Element &a, const Element &x) const = 0;

	//@} Inplace Arithmetic Operations
    
	/** @name Input/Output Operations */
	//@{
    
	/** Print ring.
	 * @return output stream to which ring is written.
	 * @param  os  output stream to which ring is written.
	 */
	virtual std::ostream &write (std::ostream &os) const = 0;
    
	/** Read ring.
	 * @return input stream from which ring is read.
	 * @param  is  input stream from which ring is read.
	 */
	virtual std::istream &read (std::istream &is) = 0;
    
	/** Print ring element.
	 * This function assumes the ring element has already been 
	 * constructed and initialized.
	 *
	 * In this implementation, this means for the <tt>
	 * _elem_ptr</tt> for x exists and does not point to
	 * null.
	 *
	 * @return output stream to which ring element is written.
	 * @param  os  output stream to which ring element is written.
	 * @param  x   ring element.
	 */
	virtual std::ostream &write (std::ostream &os, const Element &x) const = 0;
    
	/** Read ring element.
	 * This function assumes the ring element has already been 
	 * constructed and initialized.
	 *
	 * In this implementation, this means for the <tt>
	 * _elem_ptr</tt> for x exists and does not point to
	 * null.
	 *
	 * @return input stream from which ring element is read.
	 * @param  is  input stream from which ring element is read.
	 * @param  x   ring element.
	 */
	virtual std::istream &read (std::istream &is, Element &x) const = 0;
    
	//@} Input/Output Operations
	/** @name Standard elements
	 */
	//@{
	
	/// Return a reference to the zero-element of the ring
	virtual const Element &zero () const = 0;

	/// Return a reference to the one-element of the ring
	virtual const Element &one () const = 0;

	/// Return a reference to the negative of the one-element of the ring
	virtual const Element &minusOne () const = 0;

	//@}
	//@} Common Object Interface
    
}; // class RingInterface
  
} // namespace LinBox

#endif // __LINBOX_RING_INTERFACE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
