/* lela/matrix/sparse.tcc
 * Copyright 2001-2002 Bradford Hovinen
 *           1999-2001 William J Turner,
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 * Based on sparse-base.h by William J Turner <wjturner@math.ncsu.edu>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_SPARSE_TCC
#define __LELA_MATRIX_SPARSE_TCC

#include "lela/lela-config.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>

#include "lela/util/debug.h"
#include "lela/matrix/sparse.h"
#include "lela/vector/traits.h"
#include "lela/matrix/submatrix.h"
#include "lela/util/commentator.h"

namespace LELA
{

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse>
	::setEntry (size_t i, size_t j, const Element &value) 
{
        typedef typename Row::value_type value_type;
	Row &v = _A[i];
	typename Row::iterator iter;
        
	if (v.size () == 0) {
		v.push_back (value_type (j, value));
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorUtils::FindSparseEntryLB ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, value_type (j, value));
                else
                    	iter->second = value;
 	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse>
	::eraseEntry (size_t i, size_t j) 
{
	Row &v = _A[i];
	typename Row::iterator iter;
	
	iter = std::lower_bound (v.begin (), v.end (), j, VectorUtils::FindSparseEntryLB ());

	if (iter != v.end () && iter->first == j)
		v.erase (iter);
}

template <class Element, class Row>
bool SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse>
	::getEntry (Element &x, size_t i, size_t j) const
{
	const Row &v = _A[i];
	typename Row::const_iterator iter;

	if (v.size () == 0)
		return false;
	else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorUtils::FindSparseEntryLB ());

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
Vector &SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse>
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
SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse>
	&SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse>::transpose (SparseMatrix &AT) const
{
	typename Row::const_iterator j;

	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row)
		for (j = i->begin (); j != i->end (); ++j)
			AT._A[j->first].push_back (std::pair<size_t, Element> (row, j->second));

	return AT;
}

} // namespace LELA

#endif // __LELA_MATRIX_SPARSE_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
