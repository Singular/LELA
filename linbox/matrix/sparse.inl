/* linbox/matrix/sparse.inl
 * Copyright 2001-2002 Bradford Hovinen
 *           1999-2001 William J Turner,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Based on sparse-base.h by William J Turner <wjturner@math.ncsu.edu>
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.inl to matrix/sparse.inl
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_sparse_INL
#define __LINBOX_matrix_sparse_INL

#include "linbox/linbox-config.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>

#include "linbox/matrix/sparse.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"
#include "linbox/vector/sparse-subvector.h"
#include "linbox/matrix/submatrix.h"
#include <linbox/util/commentator.h>

namespace LinBox
{

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseVectorTag>
	::setEntry (size_t i, size_t j, const Element &value) 
{
        typedef typename Row::value_type value_type;
	Row &v = _A[i];
	typename Row::iterator iter;
        
	if (v.size () == 0) {
		v.push_back (value_type (j, value));
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, value_type (j, value));
                else
                    	iter->second = value;
 	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseVectorTag>
	::eraseEntry (size_t i, size_t j) 
{
	Row &v = _A[i];
	typename Row::iterator iter;
	
	iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries ());

	if (iter != v.end () && iter->first == j)
		v.erase (iter);
}

template <class Element, class Row>
bool SparseMatrix<Element, Row, VectorCategories::SparseVectorTag>
	::getEntry (Element &x, size_t i, size_t j) const
{
	const Row &v = _A[i];
	typename Row::const_iterator iter;

	if (v.size () == 0)
		return false;
	else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries ());

		if (iter == v.end () || iter->first != j)
			return false;
		else {
			x = iter->second;
			return true;
		}
	}
}

template <class Element, class Row>
template <class Vector>
Vector &SparseMatrix<Element, Row, VectorCategories::SparseVectorTag>
	::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			++v[j->first];
	}

	return v;
}

template <class Element, class Row>
SparseMatrix<Element, Row, VectorCategories::SparseVectorTag>
	&SparseMatrix<Element, Row, VectorCategories::SparseVectorTag>::transpose (SparseMatrix &AT) const
{
	typename Row::const_iterator j;

	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row)
		for (j = i->begin (); j != i->end (); ++j)
			AT._A[j->first].push_back (std::pair<size_t, Element> (row, j->second));

	return AT;
}

template <class Element, class Row, class Trait>
class SubvectorFactory<SparseMatrix<Element, Row, Trait> >
{
    public:
	typedef SparseSubvector<typename SparseMatrix<Element, Row, Trait>::Row, Trait> RowSubvector;
	typedef SparseSubvector<typename SparseMatrix<Element, Row, Trait>::ConstRow, Trait> ConstRowSubvector;

	RowSubvector MakeRowSubvector (Submatrix<SparseMatrix<Element, Row, Trait> > &M, typename SparseMatrix<Element, Row, Trait>::RowIterator &pos)
		{ return RowSubvector (*pos, M.startCol (), M.startCol () + M.coldim ()); }

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<SparseMatrix<Element, Row, Trait> > &M, typename SparseMatrix<Element, Row, Trait>::ConstRowIterator &pos)
		{ return ConstRowSubvector (*pos, M.startCol (), M.startCol () + M.coldim ()); }
};

template <class Element, class Row, class Trait>
class SubvectorFactory<const SparseMatrix<Element, Row, Trait> >
{
    public:
	typedef SparseSubvector<typename SparseMatrix<Element, Row, Trait>::ConstRow, Trait> RowSubvector;
	typedef SparseSubvector<typename SparseMatrix<Element, Row, Trait>::ConstRow, Trait> ConstRowSubvector;

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<const SparseMatrix<Element, Row, Trait> > &M, typename SparseMatrix<Element, Row, Trait>::ConstRowIterator &pos)
		{ return ConstRowSubvector (*pos, M.startCol (), M.startCol () + M.coldim ());; }
};

template <class Row>
class HybridSubvectorFactory
{
    public:
	typedef HybridSubvectorFactory<Row> Self_t;

	typedef SparseSubvector<typename SparseMatrix<bool, Row, VectorCategories::HybridZeroOneVectorTag>::ConstRow, HybridSubvectorWordAlignedTag> RowSubvector;
	typedef RowSubvector ConstRowSubvector;

	ConstRowSubvector MakeConstRowSubvector (const Submatrix<const SparseMatrix<bool, Row, VectorCategories::HybridZeroOneVectorTag>, Self_t> &M,
						 typename SparseMatrix<bool, Row, VectorCategories::HybridZeroOneVectorTag>::ConstRowIterator &pos)
		{ return ConstRowSubvector (*pos, M.startCol () >> WordTraits<typename Row::word_type>::logof_size,
					    (M.startCol () + M.coldim ()) >> WordTraits<typename Row::word_type>::logof_size); }
};

} // namespace LinBox

#endif // __LINBOX_matrix_sparse_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
