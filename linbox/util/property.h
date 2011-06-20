/* linbox/util/property.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * C++-structure which mimics C#-style properties
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_UTIL_PROPERTY_H
#define __LINBOX_UTIL_PROPERTY_H

#include <vector>

namespace LinBox
{

/// This structure mimics a C#-style property
template <class Iterator>
struct Property
{
	Iterator _i;

	typedef typename std::iterator_traits<Iterator>::value_type value_type;

	Property () {}
	Property (Iterator i) : _i (i) {}

	Property &operator = (const value_type &v)
		{ *_i = v; return *this; }

	Property &operator = (const Property &v)
		{ *_i = *v._i; return *this; }

	Property &operator += (const value_type &v)
		{ *_i += v; return *this; }

	Property &operator += (const Property &v)
		{ *_i += *v._i; return *this; }

	Property &operator -= (const value_type &v)
		{ *_i -= v; return *this; }

	Property &operator -= (const Property &v)
		{ *_i -= *v._i; return *this; }

	Property &operator *= (const value_type &v)
		{ *_i *= v; return *this; }

	Property &operator *= (const Property &v)
		{ *_i *= *v._i; return *this; }

	Property &operator %= (const value_type &v)
		{ *_i %= v; return *this; }

	Property &operator %= (const Property &v)
		{ *_i %= *v._i; return *this; }

	Property &operator &= (const value_type &v)
		{ *_i &= v; return *this; }

	Property &operator &= (const Property &v)
		{ *_i &= *v._i; return *this; }

	Property &operator |= (const value_type &v)
		{ *_i |= v; return *this; }

	Property &operator |= (const Property &v)
		{ *_i |= *v._i; return *this; }

	Property &operator ^= (const value_type &v)
		{ *_i ^= v; return *this; }

	Property &operator ^= (const Property &v)
		{ *_i ^= *v._i; return *this; }

	operator value_type () const
		{ return *_i; }
};

} // namespace LinBox

#endif // __LINBOX_UTIL_PROPERTY_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
