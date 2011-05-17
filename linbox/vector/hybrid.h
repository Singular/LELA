/* linbox/vector/hybrid.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_VECTOR_HYBRID_H
#define __LINBOX_VECTOR_HYBRID_H

#include "linbox/vector/sparse.h"

namespace LinBox
{

template <class _Endianness, class IndexType, class WordType>
class HybridVector : public SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> >
{
public:
	typedef VectorCategories::HybridZeroOneVectorTag VectorCategory; 

	typedef _Endianness Endianness;
	typedef IndexType index_type;
	typedef WordType word_type;

	HybridVector () {}
		
	template <class IV, class EV>
	HybridVector (IV &iv, EV &ev)
		: SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> > (iv, ev)
	{}

	template <class IIt, class EIt>
	HybridVector (IIt idx_begin, IIt idx_end, EIt elt_begin)
		: SparseVector<WordType, std::vector<IndexType>, std::vector<WordType> > (idx_begin, idx_end, elt_begin)
	{}
};

} // namespace LinBox

#endif // __LINBOX_VECTOR_HYBRID_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
