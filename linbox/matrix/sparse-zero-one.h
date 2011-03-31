/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/sparse-zero-one.inl
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of SparseMatrix and helpers for 0-1 matrices
 * 
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_sparse_zero_one_INL
#define __LINBOX_matrix_sparse_zero_one_INL

#include "linbox/matrix/sparse.h"
#include "linbox/vector/bit-vector.h"

namespace LinBox {

/* Specialization for sparse zero-one vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::SparseZeroOneVectorTag >
{
public:
	
	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrix<typename _Tp1::Element, _R1, VectorCategories::SparseZeroOneVectorTag> other; };

	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix (const SparseMatrix<Element, Row> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
    
    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType> &A)
		: _A(A._m), _m (A._m), _n (A._n) {
			typename Rep::iterator i;
			typename SparseMatrix<Element, VectorType>::ConstRowIterator i_A;

			for (i = _A.begin (), i_A = A.rowBegin (); i != _A.end (); ++i, ++i_A) {
				i->resize (i_A->size ());
				std::copy (i_A->begin (), i_A->end (), i->begin ());
			}
		}

	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

 	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const
	{ 
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s += LinBox::RawVector<_Element>::size (*it);

		return s;
        }

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void setEntry (size_t i, size_t j, const Element &value);
	void eraseEntry (size_t i, size_t j);
	bool getEntry (Element &x, size_t i, size_t j) const;

	ConstRowIterator rowBegin () const 
		{ return _A.begin (); }
	ConstRowIterator rowEnd () const
		{ return _A.end (); }
	RowIterator rowBegin ()
		{ return _A.begin (); }
	RowIterator rowEnd ()
		{ return _A.end (); }

	Row &getRow (size_t i) { return _A[i]; }
	Row &operator [] (size_t i) { return _A[i]; }
	ConstRow &operator [] (size_t i) const { return _A[i]; }

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::SparseZeroOneVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::SparseZeroOneVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrix &transpose (SparseMatrix &AT) const;

    protected:

	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrix;
};

/* Specialization for hybrid zero-one vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::HybridZeroOneVectorTag >
{
public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
	struct rebind
	{ typedef SparseMatrix<typename _Tp1::Element, _R1, VectorCategories::HybridZeroOneVectorTag> other; };

	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix (const SparseMatrix<Element, Row> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
    
	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType> &A)
		: _A(A._m), _m (A._m), _n (A._n)
	{
		typename Rep::iterator meit = this->_A.begin ();
		typename SparseMatrix<Element, VectorType>::Rep::const_iterator copit = A._A.begin ();

		for( ; meit != this->_A.end(); ++meit, ++copit)
			LinBox::RawVector<Element>::convert(*meit, *copit);
	}

	SparseMatrix (VectorStream<Row> &vs)
		: _A (vs.size ()), _m (vs.size ()), _n (vs.dim ())
	{
		for (RowIterator i = rowBegin (); i != rowEnd (); ++i)
			vs >> *i;
	}

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	size_t size () const
	{
		size_t s = 0;

		for (typename Rep::const_iterator it = _A.begin (); it != _A.end (); ++it)
			s+= LinBox::RawVector<_Element>::size (*it);

		return s;
	}

	void resize (size_t m, size_t n)
	{
		_m = m; _n = n;
		_A.resize (m);
	}

	void setEntry (size_t i, size_t j, const Element &value);
	void eraseEntry (size_t i, size_t j) {}
	bool getEntry (Element &x, size_t i, size_t j) const;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	ConstRowIterator rowBegin () const 
		{ return _A.begin (); }
	ConstRowIterator rowEnd () const
		{ return _A.end (); }
	RowIterator rowBegin ()
		{ return _A.begin (); }
	RowIterator rowEnd ()
		{ return _A.end (); }

	Row &getRow (size_t i) { return _A[i]; }
	Row &operator [] (size_t i) { return _A[i]; }
	ConstRow &operator [] (size_t i) const { return _A[i]; }

	typedef MatrixRawIterator<ConstRowIterator, VectorCategories::HybridZeroOneVectorTag> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorCategories::HybridZeroOneVectorTag, false> RawIndexedIterator;
	typedef RawIndexedIterator ConstRawIndexedIterator;

	ConstRawIndexedIterator rawIndexedBegin() const { return ConstRawIndexedIterator (rowBegin (), 0, rowEnd (), coldim ()); }
        ConstRawIndexedIterator rawIndexedEnd() const   { return ConstRawIndexedIterator (rowEnd (), rowdim (), rowEnd (), coldim ()); }

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrix &transpose (SparseMatrix &AT) const;

protected:

	Rep               _A;
	size_t            _m;
	size_t            _n;

	template<class F, class R, class T> friend class SparseMatrix;
};

template <class Element, class Row>
bool SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag>::getEntry (Element &x, size_t i, size_t j) const
{
	const Row &v = (*this)[i];

	typename Row::const_iterator idx = std::lower_bound (v.begin (), v.end (), j);

	if (idx != v.end () && *idx == j) {
		x = true;
		return true;
	} else
		return false;
}

template <class Element, class Row>
bool SparseMatrix<Element, Row, VectorCategories::HybridZeroOneVectorTag>::getEntry (Element &x, size_t i, size_t j) const
{
	const Row &v = (*this)[i];

	typename Row::const_iterator idx;

	idx = std::lower_bound (v.begin (), v.end (), j >> WordTraits<typename Row::word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

	if (idx != v.end () && idx->first == j >> WordTraits<typename Row::word_type>::logof_size) {
		x = idx->second & Row::Endianness::e_j (j & WordTraits<typename Row::word_type>::pos_mask);
		return true;
	} else
		return false;
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag>::eraseEntry (size_t i, size_t j)
{
	Row &v = (*this)[i];

	typename Row::iterator idx = std::lower_bound (v.begin (), v.end (), j);

	if (idx != v.end () && *idx == j)
		v.erase (idx);
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	Row &v = _A[i];
	typename Row::iterator iter;

	if (value && v.size () == 0) {
		v.push_back (j);
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j);

		if (value && (iter == v.end () || *iter != j))
			iter = v.insert (iter, j);
		else if (!value && iter != v.end () && *iter == j)
			v.erase (iter);
	}
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::HybridZeroOneVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	Row &v = _A[i];
	typename Row::iterator it;

	typename Row::word_type m = Row::Endianness::e_j (j & WordTraits<typename Row::word_type>::pos_mask);

	if (value && v.empty ()) {
		v.push_back (typename Row::value_type (j >> WordTraits<typename Row::word_type>::logof_size, m));
	} else {
		it = std::lower_bound (v.begin (), v.end (), j >> WordTraits<typename Row::word_type>::logof_size, VectorWrapper::CompareSparseEntries ());

		if (it == v.end () || it->first != (j >> WordTraits<typename Row::word_type>::logof_size)) {
			if (value)
				it = v.insert (it, typename Row::value_type (j >> WordTraits<typename Row::word_type>::logof_size, m));
		}
		else {
			if (value)
				it->second |= m;
			else {
				it->second &= ~m;

				if (!it->second)
					v.erase (it);
			}
		}
	}
}

} // namespace LinBox

#endif // __LINBOX_matrix_sparse_zero_one_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
