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

template <class _Element, class Row>
class SparseMatrixWriteHelper<_Element, Row, VectorCategories::SparseZeroOneVectorTag >
{
    public:
	typedef _Element Element;

	// Dummy class to avoid code duplication
	class NoField 
	{
	    public:
		typedef _Element Element;

		std::istream &read (std::istream &stream, Element &elt) const
			{ return stream >> elt; }
		std::ostream &write (std::ostream &stream, const Element &elt) const
			{ return stream << elt; }
	};

	template <class Field, class Trait>
	static std::ostream &write (const SparseMatrix<Element, Row, Trait> &A, std::ostream &os, const Field &F, FileFormatTag format);
};

template <class _Element, class Row>
class SparseMatrixWriteHelper<_Element, Row, VectorCategories::HybridZeroOneVectorTag >
{
    public:
	typedef _Element Element;

	// Dummy class to avoid code duplication
	class NoField 
	{
	    public:
		typedef _Element Element;

		std::istream &read (std::istream &stream, Element &elt) const
			{ return stream >> elt; }
		std::ostream &write (std::ostream &stream, const Element &elt) const
			{ return stream << elt; }
	};

	template <class Field, class Trait>
	static std::ostream &write (const SparseMatrix<Element, Row, Trait> &A, std::ostream &os, const Field &F, FileFormatTag format);
};

/* Specialization for sparse zero-one vectors */

template <class _Element, class _Row>
class SparseMatrix<_Element, _Row, VectorCategories::SparseZeroOneVectorTag >
{
public:
	
	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

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

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrix ( MatrixStream<Field>& ms );

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const { 
            size_t s(0);
            for(typename Rep::const_iterator it = _A.begin(); it != _A.end(); ++it)
                s+= LinBox::RawVector<_Element>::size(*it);
            return s;
        }

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::read
			  (*this, is, F, format); }
	std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::read
			  (*this, is, typename SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::NoField (),
			   format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::write
			  (*this, os, F, format); }
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::write
			  (*this, os, SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::NoField (),
			   format); }

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
			{ return x = getEntry (i, j); }

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

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrix &transpose (SparseMatrix &AT) const;

    protected:

	friend class SparseMatrixWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag >;
	friend class SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag >;

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
		: _A(A._m), _m (A._m), _n (A._n) {
		typename Rep::iterator meit = this->_A.begin();
		typename SparseMatrix<Element, VectorType>::Rep::const_iterator copit = A._A.begin();
		for( ; meit != this->_A.end(); ++meit, ++copit)
			LinBox::RawVector<Element>::convert(*meit, *copit);
	}

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrix ( MatrixStream<Field>& ms );

	~SparseMatrix () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const { 
		size_t s(0);
		for(typename Rep::const_iterator it = _A.begin(); it != _A.end(); ++it)
			s+= LinBox::RawVector<_Element>::size(*it);
		return s;
	}

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag>::read
				(*this, is, F, format); }
	std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag>::read
				(*this, is, typename SparseMatrixReadWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag>::NoField (),
				 format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag>::write
				(*this, os, F, format); }
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag>::write
				(*this, os, typename SparseMatrixReadWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag>::NoField (),
				 format); }

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
		{ return x = getEntry (i, j); }

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

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrix &transpose (SparseMatrix &AT) const;

protected:

	friend class SparseMatrixWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag >;
	friend class SparseMatrixReadWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag >;

	Rep               _A;
	size_t            _m;
	size_t            _n;

	template<class F, class R, class T> friend class SparseMatrix;
};

template <class Element, class Row>
template <class Field, class Trait>
std::ostream &SparseMatrixWriteHelper<Element, Row, VectorCategories::SparseZeroOneVectorTag >
	::write (const SparseMatrix<Element, Row, Trait> &A, std::ostream &os, const Field &F, 
		FileFormatTag format)
{
	typename SparseMatrix<Element, Row, Trait>::Rep::const_iterator i;
	typename Row::const_iterator j_idx;
	typename Field::Element zero;
	size_t i_idx, j_idx_1, col_idx;
	//int col_width;
	integer c;
        bool firstrow;
        
	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	typename Field::Element one;
	F.init (one, 1);

	switch (format) {
	    case FORMAT_DETECT:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "format != FORMAT_DETECT");
		break;

	    case FORMAT_TURNER:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->begin ();
			     j_idx != i->end ();
			     ++j_idx)
			{
				os << i_idx << ' ' << *j_idx << ' ';
				F.write (os, one) << std::endl;
			}
		}

		break;

	    case FORMAT_ONE_BASED:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->begin (); j_idx != i->end (); ++j_idx) {
				os << i_idx + 1 << ' ' << *j_idx + 1 << ' ';
				F.write (os, one) << std::endl;
			}
		}

		break;

	    case FORMAT_GUILLAUME:
		os << A._m << ' ' << A._n << " M" << std::endl;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->begin (); j_idx != i->end (); ++j_idx) {
				os << i_idx + 1 << ' ' << *j_idx + 1 << ' ';
				F.write (os, one) << std::endl;
			}
		}

		os << "0 0 0" << std::endl;

		break;

	    case FORMAT_MAPLE:
		F.init (zero, 0);
                firstrow=true;

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			if (firstrow) {
                            os << "[";
                            firstrow =false;
                        } else 
                             os << ", [";
                           
			j_idx = i->begin ();

			for (j_idx_1 = 0; j_idx_1 < A._n; j_idx_1++) {
				if (j_idx == i->end () || j_idx_1 != *j_idx)
					F.write (os, zero);
				else {
					F.write (os, one);
					++j_idx;
				}

				if (j_idx_1 < A._n - 1)
					os << ", ";
			}

			os << "]";
		}

		os << "]";

		break;

	    case FORMAT_MATLAB:
		F.init (zero, 0);

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			j_idx = i->begin ();

			for (j_idx_1 = 0; j_idx_1 < A._n; j_idx_1++) {
				if (j_idx == i->end () || j_idx_1 != *j_idx)
					F.write (os, zero);
				else {
					F.write (os, one);
					++j_idx;
				}

				if (j_idx_1 < A._n - 1)
					os << ", ";
			}

			os << "; ";
		}

		os << "]";

		break;

	    case FORMAT_PRETTY:
		//F.characteristic (c);
		//col_width = (int) ceil (log ((double) c) / M_LN10);
		F.init (zero, 0);

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			os << "  [ ";

			j_idx = i->begin ();

			for (col_idx = 0; col_idx < A._n; col_idx++) {
				//os.width (col_width);

				if (j_idx == i->end () || col_idx != *j_idx)
//					F.write (os, zero);
					os << '.';
				else {
					F.write (os, one);
					++j_idx;
				}

				os << ' ';
			}

			os << ']' << std::endl;
		}

		break;
	    case FORMAT_MAGMACPT:
		os << "sparse matrix written in MagmaCpt form is not implemented" << std::endl;
		break;
	}

	return os;
}

template <class Element, class Row>
template <class Field, class Trait>
std::ostream &SparseMatrixWriteHelper<Element, Row, VectorCategories::HybridZeroOneVectorTag >
	::write (const SparseMatrix<Element, Row, Trait> &A, std::ostream &os, const Field &F, 
		 FileFormatTag format)
{
	typedef typename std::iterator_traits<typename Row::second_type::const_word_iterator>::value_type word_type;

	typename SparseMatrix<Element, Row, Trait>::Rep::const_iterator i;
	typename Row::first_type::const_iterator j_idx;
	typename Row::second_type::const_word_iterator j_elt;
	typename Field::Element zero;
	typedef typename Row::second_type::Endianness Endianness;
	size_t i_idx, j_idx_1, col_idx;
	//int col_width;
	integer c;
        bool firstrow;

	word_type mask;
	size_t t;
        
	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	typename Field::Element one;
	F.init (one, 1);

	switch (format) {
	    case FORMAT_DETECT:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "format != FORMAT_DETECT");
		break;

	    case FORMAT_TURNER:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->first.begin (), j_elt = i->second.wordBegin ();
			     j_idx != i->first.end ();
			     ++j_idx, ++j_elt)
			{
				for (t = 0, mask = Endianness::e_0; mask != 0; mask = Endianness::shift_right (mask, 1), ++t) {
					if (*j_elt & mask) {
						os << i_idx << ' ' << (static_cast<size_t> (*j_idx) << WordTraits<word_type>::logof_size) + t << ' ';
						F.write (os, one) << std::endl;
					}
				}
			}
		}

		break;

	    case FORMAT_ONE_BASED:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->first.begin (), j_elt = i->second.wordBegin ();
			     j_idx != i->first.end ();
			     ++j_idx, ++j_elt)
			{
				for (t = 0, mask = Endianness::e_0; mask != 0; mask = Endianness::shift_right (mask, 1), ++t) {
					if (*j_elt & mask) {
						os << i_idx + 1 << ' ' << (static_cast<size_t> (*j_idx) << WordTraits<word_type>::logof_size) + t + 1 << ' ';
						F.write (os, one) << std::endl;
					}
				}
			}
		}

		break;

	    case FORMAT_GUILLAUME:
		os << A._m << ' ' << A._n << " M" << std::endl;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->first.begin (), j_elt = i->second.wordBegin ();
			     j_idx != i->first.end ();
			     ++j_idx, ++j_elt)
			{
				for (t = 0, mask = Endianness::e_0; mask != 0; mask = Endianness::shift_right (mask, 1), ++t) {
					if (*j_elt & mask) {
						os << i_idx + 1 << ' ' << (static_cast<size_t> (*j_idx) << WordTraits<word_type>::logof_size) + t + 1 << ' ';
						F.write (os, one) << std::endl;
					}
				}
			}
		}

		os << "0 0 0" << std::endl;

		break;

	    case FORMAT_MAPLE:
		F.init (zero, 0);
                firstrow=true;

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			if (firstrow) {
                            os << "[";
                            firstrow =false;
                        } else 
                             os << ", [";
                           
			j_idx = i->first.begin ();
			j_elt = i->second.wordBegin ();

			mask = Endianness::e_0;
			t = 0;

			for (j_idx_1 = 0; j_idx_1 < A._n; j_idx_1++) {
				if (mask == 0 && j_idx != i->first.end ()) {
					mask = Endianness::e_0;
					++t;

					if (*j_idx < t) {
						++j_idx;
						++j_elt;
					}
				}

				if (j_idx == i->first.end () || t != *j_idx || !(*j_elt & mask))
					F.write (os, zero);
				else
					F.write (os, one);

				if (j_idx_1 < A._n - 1)
					os << ", ";

				mask = Endianness::shift_right (mask, 1);
			}

			os << "]";
		}

		os << "]";

		break;

	    case FORMAT_MATLAB:
		F.init (zero, 0);

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			j_idx = i->first.begin ();
			j_elt = i->second.wordBegin ();

			mask = Endianness::e_0;
			t = 0;

			for (j_idx_1 = 0; j_idx_1 < A._n; j_idx_1++) {
				if (mask == 0 && j_idx != i->first.end ()) {
					mask = Endianness::e_0;
					++t;

					if (*j_idx < t) {
						++j_idx;
						++j_elt;
					}
				}

				if (j_idx == i->first.end () || t != *j_idx || !(*j_elt & mask))
					F.write (os, zero);
				else
					F.write (os, one);

				if (j_idx_1 < A._n - 1)
					os << ", ";

				mask = Endianness::shift_right (mask, 1);
			}

			os << "; ";
		}

		os << "]";

		break;

	    case FORMAT_PRETTY:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			os << "  [ ";

			j_idx = i->first.begin ();
			j_elt = i->second.wordBegin ();

			mask = Endianness::e_0;
			t = 0;

			for (col_idx = 0; col_idx < A._n; col_idx++) {
				if (mask == 0 && j_idx != i->first.end ()) {
					mask = Endianness::e_0;
					++t;

					if (*j_idx < t) {
						++j_idx;
						++j_elt;
					}
				}

				if (j_idx == i->first.end () || t != *j_idx || !(*j_elt & mask))
					os << '.';
				else
					F.write (os, one);

				os << ' ';

				mask = Endianness::shift_right (mask, 1);
			}

			os << ']' << std::endl;
		}

		break;

	    case FORMAT_SAGE:
	        os << "matrix(R,[" << std::endl;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			os << "  [ ";

			j_idx = i->first.begin ();
			j_elt = i->second.wordBegin ();

			mask = Endianness::e_0;
			t = 0;

			for (col_idx = 0; col_idx < A._n; col_idx++) {
				if (mask == 0 && j_idx != i->first.end ()) {
					mask = Endianness::e_0;
					++t;

					if (*j_idx < t) {
						++j_idx;
						++j_elt;
					}
				}

				if (j_idx == i->first.end () || t != *j_idx || !(*j_elt & mask))
					os << '.';
				else
					F.write (os, one);

				if (col_idx < A._n - 1)
					os << ", ";

				mask = Endianness::shift_right (mask, 1);
			}

			if (i_idx < A._m - 1)
				os << "]," << std::endl;
			else
				os << "] ])" << std::endl;
		}

		break;
	    case FORMAT_MAGMACPT:
		os << "sparse matrix written in MagmaCpt form is not implemented" << std::endl;
		break;
	}

	return os;
}

template <class Element, class Row>
void SparseMatrix<Element, Row, VectorCategories::SparseZeroOneVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	while (_A.size() < i + 1)
		_A.push_back (Row());

	_m = _A.size(); 

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
	while (_A.size() < i + 1)
		_A.push_back (Row());

	_m = _A.size(); 

	Row &v = _A[i];
	typename Row::first_type::iterator it_idx;
	typename Row::second_type::word_iterator it_elt;
	typedef typename std::iterator_traits<typename Row::second_type::word_iterator>::value_type word_type;

	typename Row::second_type::word_iterator::value_type m = 1ULL << (j & WordTraits<word_type>::pos_mask);

	if (value && v.first.size () == 0) {
		v.first.push_back (j & ~WordTraits<word_type>::pos_mask);
		v.second.push_word_back (m);
	} else {
		it_idx = std::lower_bound (v.first.begin (), v.first.end (), j & ~WordTraits<word_type>::pos_mask);

		if (it_idx == v.first.end () || *it_idx != j & ~WordTraits<word_type>::pos_mask) {
			if (value) {
				it_idx = v.first.insert (it_idx, j & ~WordTraits<word_type>::pos_mask);
				v.second.insertWord (v.second.wordBegin () + (it_idx - v.first.begin ()), m);
			}
		}
		else {
			if (value)
				*(v.second.wordBegin () + (it_idx - v.first.begin ())) |= m;
			else {
				it_elt = v.second.wordBegin () + (it_idx - v.first.begin ());
				*it_elt &= ~m;

				if (!*it_elt) {
					v.first.erase (it_idx);
					v.second.eraseWord (it_elt);
				}
			}
		}
	}
}

} // namespace LinBox

#endif // __LINBOX_matrix_sparse_zero_one_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
