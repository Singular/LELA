/* lela/matrix/sparse-zero-one.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of SparseMatrix and helpers for 0-1 matrices
 * 
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_SPARSE_ZERO_ONE_H
#define __LELA_MATRIX_SPARSE_ZERO_ONE_H

#include <vector>

#include "lela/matrix/sparse.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/sparse-subvector.h"

namespace LELA {

/* Specialization for sparse zero-one vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorRepresentationTypes::Sparse01>
{
public:
	
	typedef _Element Element;
	typedef _Row Row;
	typedef SparseMatrix<Element, Row, VectorRepresentationTypes::Sparse01> Self_t;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;
	typedef MatrixIteratorTypes::Row IteratorType; 
	typedef MatrixStorageTypes::Rows StorageType;

	typedef Submatrix<Self_t> SubmatrixType;
	typedef Submatrix<const Self_t> ConstSubmatrixType;
	typedef Submatrix<Self_t> AlignedSubmatrixType;
	typedef Submatrix<const Self_t> ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = 1;

	typedef Self_t ContainerType;

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

        SparseMatrix () : _m (0), _n (0) {}
	SparseMatrix (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix (const SparseMatrix<Element, Row> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
    
    	template<class VectorType>
	SparseMatrix (const SparseMatrix<Element, VectorType> &A)
		: _A (A._m), _m (A._m), _n (A._n)
	{
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
			s += LELA::RawVector<_Element>::size (*it);

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

	typedef MatrixRawIterator<ConstRowIterator, VectorRepresentationTypes::Sparse01> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorRepresentationTypes::Sparse01, false> RawIndexedIterator;
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
class SparseMatrix<_Element, _Row, VectorRepresentationTypes::Hybrid01>
{
public:

	typedef _Element Element;
	typedef _Row Row;
	typedef SparseMatrix<Element, Row, VectorRepresentationTypes::Hybrid01> Self_t;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;
	typedef MatrixIteratorTypes::Row IteratorType; 
	typedef MatrixStorageTypes::Rows StorageType;

	typedef Submatrix<Self_t> SubmatrixType;
	typedef Submatrix<const Self_t> ConstSubmatrixType;
	typedef Submatrix<Self_t, AlignedSubmatrixTag> AlignedSubmatrixType;
	typedef Submatrix<const Self_t, AlignedSubmatrixTag> ConstAlignedSubmatrixType;

	static const size_t rowAlign = 1;
	static const size_t colAlign = WordTraits<typename Row::word_type>::bits;

	typedef Self_t ContainerType;

        SparseMatrix () : _m (0), _n (0) {}
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
			LELA::RawVector<Element>::convert(*meit, *copit);
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
			s+= LELA::RawVector<_Element>::size (*it);

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

	typedef MatrixRawIterator<ConstRowIterator, VectorRepresentationTypes::Hybrid01> RawIterator;
	typedef RawIterator ConstRawIterator;
    
	ConstRawIterator rawBegin () const { return ConstRawIterator (rowBegin (), 0, rowEnd (), coldim ()); }
	ConstRawIterator rawEnd () const   { return ConstRawIterator (rowEnd (), 0, rowEnd (), coldim ()); }

	typedef MatrixRawIndexedIterator<ConstRowIterator, VectorRepresentationTypes::Hybrid01, false> RawIndexedIterator;
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

} // namespace LELA

#include "lela/matrix/sparse-zero-one.tcc"

#endif // __LELA_MATRIX_SPARSE_ZERO_ONE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
