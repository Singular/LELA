/* linbox/matrix/sparse-zero-one.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Specialisation of SparseMatrix and helpers for 0-1 matrices
 * 
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_sparse_zero_one_H
#define __LINBOX_matrix_sparse_zero_one_H

#include <vector>

#include "linbox/matrix/sparse.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/sparse-subvector.h"

namespace LinBox {

/* Specialization for sparse zero-one vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::SparseZeroOneVectorTag>
{
public:
	
	typedef _Element Element;
	typedef _Row Row;
	typedef SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag> Self_t;
	typedef const Row ConstRow;
	typedef std::vector<Row> Rep;
	typedef MatrixCategories::RowMatrixTag MatrixCategory; 
	typedef SparseMatrixTag<bool, Row, VectorCategories::SparseZeroOneVectorTag> Tag;

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
class SparseMatrix<_Element, _Row, VectorCategories::HybridZeroOneVectorTag>
{
public:

	typedef _Element Element;
	typedef _Row Row;
	typedef SparseMatrix<Element, Row, VectorCategories::HybridZeroOneVectorTag> Self_t;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;
	typedef MatrixCategories::RowMatrixTag MatrixCategory; 
	typedef SparseMatrixTag<bool, Row, VectorCategories::HybridZeroOneVectorTag> Tag;

	typedef Submatrix<Self_t> SubmatrixType;
	typedef Submatrix<const Self_t> ConstSubmatrixType;
	typedef Submatrix<Self_t, HybridSubvectorWordAlignedTag> AlignedSubmatrixType;
	typedef Submatrix<const Self_t, HybridSubvectorWordAlignedTag> ConstAlignedSubmatrixType;

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

template <class Row, class Trait, class Submatrix>
class SubvectorFactory<SparseMatrixTag<bool, Row, Trait>, Submatrix, typename SparseMatrix<bool, Row, Trait>::RowIterator, HybridSubvectorWordAlignedTag>
{
    public:
	typedef typename SparseMatrix<bool, Row, Trait>::RowIterator IteratorType;
	typedef SparseSubvector<typename SparseMatrix<bool, Row, Trait>::Row, HybridSubvectorWordAlignedTag> Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (*pos, M.startCol (), M.startCol () + M.coldim ()); }
};

template <class Row, class Trait, class Submatrix>
class SubvectorFactory<SparseMatrixTag<bool, Row, Trait>, Submatrix, typename SparseMatrix<bool, Row, Trait>::ConstRowIterator, HybridSubvectorWordAlignedTag>
{
    public:
	typedef typename SparseMatrix<bool, Row, Trait>::ConstRowIterator IteratorType;
	typedef SparseSubvector<typename SparseMatrix<bool, Row, Trait>::ConstRow, HybridSubvectorWordAlignedTag> Subvector;

	Subvector MakeSubvector (Submatrix &M, IteratorType &pos)
		{ return Subvector (*pos, M.startCol (), M.startCol () + M.coldim ()); }
};

} // namespace LinBox

#include "linbox/matrix/sparse-zero-one.tcc"

#endif // __LINBOX_matrix_sparse_zero_one_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
