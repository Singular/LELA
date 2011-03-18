/* linbox/matrix/sparse.inl
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
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

template <class Element, class Row, class Tag>
template <class Field>
SparseMatrix<Element,Row,Tag>
	::SparseMatrix( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
template <class Field>
SparseMatrix<Element,Row,VectorCategories::SparseSequenceVectorTag>
	::SparseMatrix( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
template <class Field>
SparseMatrix<Element,Row,VectorCategories::SparseAssociativeVectorTag>
	::SparseMatrix( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
template <class Field>
SparseMatrix<Element,Row,VectorCategories::SparseParallelVectorTag>
	::SparseMatrix( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
        typedef typename Row::value_type value_type;
	Row &v = _A[i];
	typename Row::iterator iter;
        
	if (v.size () == 0) {
		v.push_back ( value_type(j, value));                
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, value_type(j, value));
                else
                    	iter->second = value;
 	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >
	::eraseEntry (size_t i, size_t j) 
{
	Row &v = _A[i];
	typename Row::iterator iter;
	
	iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

	if (iter != v.end () && iter->first == j)
		v.erase (iter);
}

template <class Element, class Row>
Element &SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >
	::refEntry (size_t i, size_t j) 
{
	static Element zero;

	Row &v = _A[i];
	typename Row::iterator iter;

	if (v.size () == 0) {
		v.push_back (std::pair <size_t, Element> (j, zero));
		return v.front ().second;
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, std::pair <size_t, Element> (j, zero));

		return iter->second;
	}
}

template <class Element, class Row>
const Element &SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >
	::getEntry (size_t i, size_t j) const
{
	static Element zero;

	const Row &v = _A[i];
	typename Row::const_iterator iter;

	if (v.size () == 0)
		return zero;
	else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			return zero;
		else
			return iter->second;
	}
}

template <class Element, class Row>
const Element &SparseMatrix<Element, Row, VectorCategories::SparseAssociativeVectorTag >
	::getEntry (size_t i, size_t j) const
{
	static Element zero;

	const Row &v = _A[i];
	typename Row::const_iterator iter;

	if (v.size () == 0)
		return zero;
	else {
		iter = v.find (j);

		if (iter == v.end () || iter->first != j)
			return zero;
		else
			return iter->second;
	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseParallelVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	while (_A.size() < i + 1) _A.push_back(Row());
	_m = _A.size(); 

	Row &v = _A[i];
	typename Row::first_type::iterator iter;

	if (v.first.size () == 0) {
		v.first.push_back (j);
		v.second.push_back (value);
	} else {
		iter = std::lower_bound (v.first.begin (), v.first.end (), j);

		if (iter == v.first.end () || *iter != j) {
			iter = v.first.insert (iter, j);
			v.second.insert (v.second.begin () + (iter - v.first.begin ()), value);
		} else
                    	*(v.second.begin () + (iter - v.first.begin ())) = value;
	}
}

template <class Element, class Row>
Element &SparseMatrix<Element, Row, VectorCategories::SparseParallelVectorTag >
	::refEntry (size_t i, size_t j) 
{
	static Element zero;

	Row &v = _A[i];
	typename Row::first_type::iterator iter;
	typename Row::second_type::iterator iter_elt;

	if (v.first.size () == 0) {
		v.first.push_back (j);
		v.second.push_back (zero);
		return v.second.front ();
	} else {
		iter = std::lower_bound (v.first.begin (), v.first.end (), j);

		if (iter == v.first.end () || *iter != j) {
			iter = v.first.insert (iter, j);
			iter_elt = v.second.insert (v.second.begin () + (iter - v.first.begin ()), zero);
		}
		else
			iter_elt = v.second.begin () + (iter - v.first.begin ());

		return *iter_elt;
	}
}

template <class Element, class Row>
const Element &SparseMatrix<Element, Row, VectorCategories::SparseParallelVectorTag >
	::getEntry (size_t i, size_t j) const
{
	static Element zero;

	const Row &v = _A[i];
	typename Row::first_type::const_iterator iter;

	if (v.first.size () == 0)
		return zero;
	else {
		iter = std::lower_bound (v.first.begin (), v.first.end (), j);

		if (iter == v.first.end () || *iter != j)
			return zero;
		else
			return *(v.second.begin () + (iter - v.first.begin ()));
	}
}

template <class Element, class Row>
template <class Vector>
Vector &SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >::columnDensity (Vector &v) const
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
template <class Vector>
Vector &SparseMatrix<Element, Row, VectorCategories::SparseParallelVectorTag >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::first_type::const_iterator j_idx = i->first.begin ();

		for (; j_idx != i->first.end (); ++j_idx)
			++v[*j_idx];
	}

	return v;
}

template <class Element, class Row>
template <class Vector>
Vector &SparseMatrix<Element, Row, VectorCategories::SparseAssociativeVectorTag >::columnDensity (Vector &v) const
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
SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >
	&SparseMatrix<Element, Row, VectorCategories::SparseSequenceVectorTag >::transpose (SparseMatrix &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			AT._A[j->first].push_back (std::pair<size_t, Element> (row, j->second));
	}

	return AT;
}

template <class Element, class Row>
SparseMatrix<Element, Row, VectorCategories::SparseAssociativeVectorTag >
	&SparseMatrix<Element, Row, VectorCategories::SparseAssociativeVectorTag >::transpose (SparseMatrix &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			AT._A[j->first][row] = j->second;
	}

	return AT;
}

template <class Element, class Row>
SparseMatrix<Element, Row, VectorCategories::SparseParallelVectorTag >
	&SparseMatrix<Element, Row, VectorCategories::SparseParallelVectorTag >::transpose (SparseMatrix &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::first_type::const_iterator j_idx = i->first.begin ();
		typename Row::second_type::const_iterator j_elt = i->second.begin ();

		for (; j_idx != i->first.end (); ++j_idx, ++j_elt) {
			AT._A[*j_idx].first.push_back (row);
			AT._A[*j_idx].second.push_back (*j_elt);
		}
	}

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

} // namespace LinBox

#endif // __LINBOX_matrix_sparse_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
