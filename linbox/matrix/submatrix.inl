/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/submatrix.inl
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>,
 *
 * See COPYING for license information
 */

#ifndef __MATRIX_SUBMATRIX_INL
#define __MATRIX_SUBMATRIX_INL

#include "linbox/matrix/submatrix.h"

namespace std {
	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixBaseRowIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::RowSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::RowSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::RowSubvector value_type;
		typedef long difference_type;
	};

	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixBaseConstRowIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstRowSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstRowSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstRowSubvector value_type;
		typedef long difference_type;
	};

	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixBaseColIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ColSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ColSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ColSubvector value_type;
		typedef long difference_type;
	};

	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixBaseConstColIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstColSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstColSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstColSubvector value_type;
		typedef long difference_type;
	};
}

namespace LinBox {

template <class Matrix, class Trait>
class SubmatrixBaseConstRowIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixBaseConstRowIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixBaseConstRowIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixBaseConstRowIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixBaseConstRowIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixBaseConstRowIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixBaseConstRowIterator () {}
	SubmatrixBaseConstRowIterator (const SubmatrixBase<Matrix, Trait> *M, typename Matrix::ConstRowIterator pos)
		: _M (M), _pos (pos), _row_valid (false) {}
	SubmatrixBaseConstRowIterator (const SubmatrixBaseConstRowIterator &i) : _M (i._M), _pos (i._pos), _row (i._row), _row_valid (i._row_valid) {}

	SubmatrixBaseConstRowIterator &operator = (const SubmatrixBaseConstRowIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_row = i._row;
		_row_valid = i._row_valid;
		return *this;
	}

	SubmatrixBaseConstRowIterator &operator ++ () 
	{
		++_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixBaseConstRowIterator operator ++ (int) 
	{
		SubmatrixBaseConstRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixBaseConstRowIterator operator + (difference_type i) const
	{
		return SubmatrixBaseConstRowIterator (_M, _pos + i);
	}

	SubmatrixBaseConstRowIterator &operator += (difference_type i) 
	{
		_pos += i;
		_row_valid = false;
		return *this;
	}

	SubmatrixBaseConstRowIterator &operator -- () 
	{
		--_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixBaseConstRowIterator operator -- (int) 
	{
		SubmatrixBaseConstRowIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixBaseConstRowIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixBaseConstRowIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixBaseConstRowIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i); }

	const reference operator * ()
		{ update_row (); return _row; }

	const pointer operator -> ()
		{ update_row (); return &_row; }

	bool operator == (const SubmatrixBaseConstRowIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixBaseConstRowIterator &c) const 
		{ return (_pos != c._pos); }

    private:
	friend class SubmatrixBaseRowIterator<Matrix, Trait>;

	const SubmatrixBase<Matrix, Trait> *_M;
	typename Matrix::ConstRowIterator _pos;
	typename SubvectorFactory<Matrix, Trait>::ConstRowSubvector _row;
	bool _row_valid;

	inline void update_row () {
		if (!_row_valid) {
			_row = (SubvectorFactory<Matrix, Trait> ()).MakeConstRowSubvector (*_M, _pos);
			_row_valid = true;
		}
	}
};

template <class Matrix, class Trait>
class SubmatrixBaseRowIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixBaseRowIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixBaseRowIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixBaseRowIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixBaseRowIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixBaseRowIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixBaseRowIterator () {}
	SubmatrixBaseRowIterator (SubmatrixBase<Matrix, Trait> *M, typename Matrix::RowIterator pos)
		: _M (M), _pos (pos), _row_valid (false) {}
	SubmatrixBaseRowIterator (const SubmatrixBaseRowIterator &i)
		: _M (i._M), _pos (i._pos), _row (i._row), _row_valid (i._row_valid) {}

	SubmatrixBaseRowIterator &operator = (const SubmatrixBaseRowIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_row = i._row;
		_row_valid = i._row_valid;
		return *this;
	}

	SubmatrixBaseRowIterator &operator ++ () 
	{
		++_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixBaseRowIterator operator ++ (int) 
	{
		SubmatrixBaseRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixBaseRowIterator operator + (difference_type i) const
	{
		return SubmatrixBaseRowIterator (_M, _pos + i);
	}

	SubmatrixBaseRowIterator &operator += (difference_type i) 
	{
		_pos += i;
		_row_valid = false;
		return *this;
	}

	SubmatrixBaseRowIterator &operator -- () 
	{
		--_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixBaseRowIterator operator -- (int) 
	{
		SubmatrixBaseRowIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixBaseRowIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixBaseRowIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixBaseRowIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ update_row (); return _row; }

	pointer operator -> () 
		{ update_row (); return &_row; }

	bool operator == (const SubmatrixBaseRowIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixBaseRowIterator &c) const 
		{ return (_pos != c._pos); }

	operator SubmatrixBaseConstRowIterator<Matrix, Trait> ()
		{ return SubmatrixBaseConstRowIterator<Matrix, Trait> (_M, _pos); }

    private:
	friend class SubmatrixBaseConstRowIterator<Matrix, Trait>;

	SubmatrixBase<Matrix, Trait> *_M;
	typename Matrix::RowIterator _pos;
	typename SubvectorFactory<Matrix, Trait>::RowSubvector _row;
	bool _row_valid;

	inline void update_row () {
		if (!_row_valid) {
			_row = (SubvectorFactory<Matrix, Trait> ()).MakeRowSubvector (*_M, _pos);
			_row_valid = true;
		}
	}
};

template <class Matrix, class Trait>
class SubmatrixBaseConstColIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixBaseConstColIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixBaseConstColIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixBaseConstColIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixBaseConstColIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixBaseConstColIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixBaseConstColIterator () {}
	SubmatrixBaseConstColIterator (const SubmatrixBase<Matrix, Trait> *M, typename Matrix::ConstColIterator pos)
		: _M (M), _pos (pos), _col_valid (false) {}
	SubmatrixBaseConstColIterator (const SubmatrixBaseConstColIterator &i) : _M (i._M), _pos (i._pos), _col (i._col), _col_valid (i._col_valid) {}

	SubmatrixBaseConstColIterator &operator = (const SubmatrixBaseConstColIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_col = i._col;
		_col_valid = i._col_valid;
		return *this;
	}

	SubmatrixBaseConstColIterator &operator ++ () 
	{
		++_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixBaseConstColIterator operator ++ (int) 
	{
		SubmatrixBaseConstColIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixBaseConstColIterator operator + (difference_type i) const
	{
		return SubmatrixBaseConstColIterator (_M, _pos + i);
	}

	SubmatrixBaseConstColIterator &operator += (difference_type i) 
	{
		_pos += i;
		_col_valid = false;
		return *this;
	}

	SubmatrixBaseConstColIterator &operator -- () 
	{
		--_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixBaseConstColIterator operator -- (int) 
	{
		SubmatrixBaseConstColIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixBaseConstColIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixBaseConstColIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixBaseConstColIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i); }

	const reference operator * () 
		{ update_col (); return _col; }

	const pointer operator -> ()
		{ update_col (); return &_col; }

	bool operator == (const SubmatrixBaseConstColIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixBaseConstColIterator &c) const 
		{ return (_pos != c._pos); }

    private:
	friend class SubmatrixBaseColIterator<Matrix, Trait>;

	const SubmatrixBase<Matrix, Trait> *_M;
	typename Matrix::ConstColIterator _pos;
	typename SubvectorFactory<Matrix, Trait>::ConstColSubvector _col;
	bool _col_valid;

	inline void update_col () {
		if (!_col_valid) {
			_col = (SubvectorFactory<Matrix, Trait> ()).MakeConstColSubvector (*_M, _pos);
			_col_valid = true;
		}
	}
};

template <class Matrix, class Trait>
class SubmatrixBaseColIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixBaseColIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixBaseColIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixBaseColIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixBaseColIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixBaseColIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixBaseColIterator () {}
	SubmatrixBaseColIterator (SubmatrixBase<Matrix, Trait> *M, typename Matrix::ColIterator pos)
		: _M (M), _pos (pos), _col_valid (false) {}
	SubmatrixBaseColIterator (const SubmatrixBaseColIterator &i)
		: _M (i._M), _pos (i._pos), _col (i._col), _col_valid (i._col_valid) {}

	SubmatrixBaseColIterator &operator = (const SubmatrixBaseColIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_col = i._col;
		_col_valid = i._col_valid;
		return *this;
	}

	SubmatrixBaseColIterator &operator ++ () 
	{
		++_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixBaseColIterator operator ++ (int) 
	{
		SubmatrixBaseColIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixBaseColIterator operator + (difference_type i) const
	{
		return SubmatrixBaseColIterator (_M, _pos + i);
	}

	SubmatrixBaseColIterator &operator += (difference_type i) 
	{
		_pos += i;
		_col_valid = false;
		return *this;
	}

	SubmatrixBaseColIterator &operator -- () 
	{
		--_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixBaseColIterator operator -- (int) 
	{
		SubmatrixBaseColIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixBaseColIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixBaseColIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixBaseColIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ update_col (); return _col; }

	pointer operator -> () 
		{ update_col (); return &_col; }

	bool operator == (const SubmatrixBaseColIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixBaseColIterator &c) const 
		{ return (_pos != c._pos); }

	operator SubmatrixBaseConstColIterator<Matrix, Trait> ()
		{ return SubmatrixBaseConstColIterator<Matrix, Trait> (_M, _pos); }

    private:
	friend class SubmatrixBaseConstColIterator<Matrix, Trait>;

	SubmatrixBase<Matrix, Trait> *_M;
	typename Matrix::ColIterator _pos;
	typename SubvectorFactory<Matrix, Trait>::ColSubvector _col;
	bool _col_valid;

	inline void update_col () {
		if (!_col_valid) {
			_col = (SubvectorFactory<Matrix, Trait> ()).MakeColSubvector (*_M, _pos);
			_col_valid = true;
		}
	}
};

template <class Matrix, class Trait>
template <class Field>
std::istream& SubmatrixBase<Matrix, Trait>::read (std::istream &file, const Field& field)
{
	size_t i, j;

	for (i = _beg_row; i < _end_row; ++i) {
		for (j = _beg_col; j < _end_col; ++j) {
			file.ignore (1);
			field.read (file, refEntry (i, j));
		}
	}

	return file;
}

template <class Matrix, class Trait>
template <class Field>
std::ostream &SubmatrixBase<Matrix, Trait>::write (std::ostream &os, const Field& field, bool mapleFormat) const
{
	ConstRowIterator p;

	typename ConstRow::const_iterator pe;

	if (mapleFormat) os << "[";

	for (p = rowBegin (); p != rowEnd (); ++p) {
		if (mapleFormat && (p != rowBegin()))
			os << ',';
		if (mapleFormat) os << "[";

		for (pe = (*p).begin (); pe != (*p).end (); ++pe) {
			if (mapleFormat && (pe != (*p).begin()))
				os << ',';

			if (field.isZero (*pe))
				os << '.';
			else
				field.write(os, *pe);

			os << " ";
		}

		if (!mapleFormat)
			os << std::endl;
		else os << ']';
	}

	if (mapleFormat) os << ']';
	os << std::endl;

	return os;
}

template <class Matrix>
template <class Field>
std::istream& SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::read (std::istream &file, const Field& field)
{
	size_t i, j;

	for (i = _beg_row; i < _end_row; ++i) {
		for (j = _beg_col; j < _end_col; ++j) {
			file.ignore (1);
			field.read (file, refEntry (i, j));
		}
	}

	return file;
}

template <class Matrix>
template <class Field>
std::ostream &SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::write (std::ostream &os, const Field& field, bool mapleFormat) const
{
	ConstRowIterator p;

	typename ConstRow::const_iterator pe;

	if (mapleFormat) os << "[";

	for (p = rowBegin (); p != rowEnd (); ++p) {
		if (mapleFormat && (p != rowBegin()))
			os << ',';
		if (mapleFormat) os << "[";

		for (pe = (*p).begin (); pe != (*p).end (); ++pe) {
			if (mapleFormat && (pe != (*p).begin()))
				os << ',';

			if (field.isZero (*pe))
				os << '.';
			else
				field.write(os, *pe);

			os << " ";
		}

		if (!mapleFormat)
			os << std::endl;
		else os << ']';
	}

	if (mapleFormat) os << ']';
	os << std::endl;

	return os;
}

template <class Matrix, class Trait>
inline typename SubmatrixBase<Matrix, Trait>::RowIterator SubmatrixBase<Matrix, Trait>::rowBegin ()
{
	return RowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class Trait>
inline typename SubmatrixBase<Matrix, Trait>::RowIterator SubmatrixBase<Matrix, Trait>::rowEnd ()
{
	return RowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix, class Trait>
inline typename SubmatrixBase<Matrix, Trait>::ConstRowIterator SubmatrixBase<Matrix, Trait>::rowBegin () const
{
	return ConstRowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class Trait>
inline typename SubmatrixBase<Matrix, Trait>::ConstRowIterator SubmatrixBase<Matrix, Trait>::rowEnd () const
{
	return ConstRowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix>
inline typename SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::RowIterator SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::rowBegin ()
{
	return RowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix>
inline typename SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::RowIterator SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::rowEnd ()
{
	return RowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix>
inline typename SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::ConstRowIterator SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::rowBegin () const
{
	return ConstRowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix>
inline typename SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::ConstRowIterator SubmatrixBase<Matrix, MatrixCategories::RowMatrixTag>::rowEnd () const
{
	return ConstRowIterator (this, _M->rowBegin () + _end_row);
}

}

#endif // __MATRIX_SUBMATRIX_INL

