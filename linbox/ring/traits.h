/* linbox/ring/traits.h
 * Copyright 2004 Dan Roche
 * see COPYING for permissions etc.
 */

#ifndef __LINBOX_ring_traits_H
#define __LINBOX_ring_traits_H

#include <linbox/integer.h>

// Namespace in which all LinBox library code resides
namespace LinBox 
{

/*! \brief some basic information about each ring.  
  \ingroup ring

 * It will try
 * to take the information from the ring when possible, and use defaults
 * otherwise.
 * maxModulus returns the greatest modulus that is usable for a given ring, -1
 * for infinite. goodModulus returns takes an integer and returns true if that
 * integer is a valid modulus for the given ring.
 * maxExponent and goodExponent do the same for the prime power.
 */

class RingCategories {

public:

	//generic ring.
	struct GenericTag{};
	//If it is isomorphic to Z/mZ, for some m or its extensions.
	struct ModularTag : public virtual GenericTag{};
	//If it is isomorphic to Z
	struct IntegerTag : public virtual GenericTag{};
	//If it is isomorphic to Q
	struct RationalTag : public virtual GenericTag{};
};

template <class Ring>
struct ClassifyRing {
typedef	RingCategories::GenericTag categoryTag;
};

template <class Ring>
struct RingTraits
{
	typedef typename ClassifyRing<Ring>::categoryTag categoryTag;
	
	static integer& maxModulus( integer& i ) {
		return i = static_cast<integer>(Ring::getMaxModulus());
	}
	static bool goodModulus( const integer& i ) {
		integer max;
		maxModulus( max );
		if( max == -1 ) return ( i >= 2 );
		else if( max == 0 ) return ( i == 0 );
		else return ( i >= 2 && i <= max );
	}

	static integer& maxExponent( integer& i ) { return i = 1; }
	static bool goodExponent( const integer& i ) {
		integer max;
                maxExponent( max );
                if( max == -1 ) return ( i >= 1 );
                else return ( i >= 1 && i <= max );
	}
};

} // Namespace LinBox

#endif // __LINBOX_ring_traits_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
