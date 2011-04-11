/* linbox/field/gf2.inl
 * Copyright 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_gf2_INL
#define __LINBOX_field_gf2_INL

#include <iostream>
#include <time.h>

#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/hybrid.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/matrix/matrix-domain.h"

#include <cctype> //isdigit

namespace LinBox 
{ 

// Specialization of canonical vector types

template <>
class RawVector<bool>
{
    public:
	typedef BitVector<BigEndian<uint64> > Dense;
	typedef std::vector<size_t> Sparse;
	typedef HybridVector<BigEndian<uint64>, uint16, uint64> Hybrid;
};

// Vector traits for hybrid sparse-dense format
template <> 
struct GF2VectorTraits<HybridVector<BigEndian<uint64>, uint16, uint64> >
{ 
	typedef HybridVector<BigEndian<uint64>, uint16, uint64> VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

template <>
struct GF2VectorTraits<const HybridVector<BigEndian<uint64>, uint16, uint64> >
{ 
	typedef const HybridVector<BigEndian<uint64>, uint16, uint64> VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

template <> 
struct GF2VectorTraits<HybridVector<LittleEndian<uint64>, uint16, uint64> >
{ 
	typedef HybridVector<LittleEndian<uint64>, uint16, uint64> VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

template <>
struct GF2VectorTraits<const HybridVector<LittleEndian<uint64>, uint16, uint64> >
{ 
	typedef const HybridVector<LittleEndian<uint64>, uint16, uint64> VectorType;
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 
};

} // namespace LinBox

#include "linbox/vector/vector-domain-gf2.h"
#include "linbox/matrix/matrix-domain-gf2.h"

#ifdef __LINBOX_HAVE_M4RI
#  include "linbox/matrix/m4ri-matrix.h"
#endif

#endif // __LINBOX_field_gf2_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
