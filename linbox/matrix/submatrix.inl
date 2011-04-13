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

namespace LinBox {

template <class Matrix, class SubvectorFactory, class Trait>
class SubmatrixConstRowIterator {
    public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename SubvectorFactory::ConstRowSubvector &reference;
	typedef typename SubvectorFactory::ConstRowSubvector *pointer;
	typedef typename SubvectorFactory::ConstRowSubvector value_type;
	typedef ptrdiff_t difference_type;

	SubmatrixConstRowIterator () {}

	template <class It>
	SubmatrixConstRowIterator (const Submatrix<Matrix, SubvectorFactory, Trait> *M, It pos)
		: _M (M), _pos (pos), _row_valid (false) {}

	SubmatrixConstRowIterator (const SubmatrixConstRowIterator &i) : _M (i._M), _pos (i._pos), _row (i._row), _row_valid (i._row_valid) {}

	SubmatrixConstRowIterator &operator = (const SubmatrixConstRowIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_row = i._row;
		_row_valid = i._row_valid;
		return *this;
	}

	SubmatrixConstRowIterator &operator ++ () 
	{
		++_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixConstRowIterator operator ++ (int) 
	{
		SubmatrixConstRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixConstRowIterator operator + (difference_type i) const
	{
		return SubmatrixConstRowIterator (_M, _pos + i);
	}

	SubmatrixConstRowIterator &operator += (difference_type i) 
	{
		_pos += i;
		_row_valid = false;
		return *this;
	}

	SubmatrixConstRowIterator &operator -- () 
	{
		--_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixConstRowIterator operator -- (int) 
	{
		SubmatrixConstRowIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixConstRowIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixConstRowIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixConstRowIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i); }

	const reference operator * ()
		{ update_row (); return _row; }

	const pointer operator -> ()
		{ update_row (); return &_row; }

	bool operator == (const SubmatrixConstRowIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixConstRowIterator &c) const 
		{ return (_pos != c._pos); }

    private:
	friend class SubmatrixRowIterator<Matrix, SubvectorFactory, Trait>;

	const Submatrix<Matrix, SubvectorFactory, Trait> *_M;
	typename Matrix::ConstRowIterator _pos;
	typename SubvectorFactory::ConstRowSubvector _row;
	bool _row_valid;

	inline void update_row () {
		if (!_row_valid) {
			_row = (SubvectorFactory ()).MakeConstRowSubvector (*_M, _pos);
			_row_valid = true;
		}
	}
};

template <class Matrix, class SubvectorFactory, class Trait>
class SubmatrixRowIterator {
    public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename SubvectorFactory::RowSubvector &reference;
	typedef typename SubvectorFactory::RowSubvector *pointer;
	typedef typename SubvectorFactory::RowSubvector value_type;
	typedef ptrdiff_t difference_type;

	SubmatrixRowIterator () {}
	SubmatrixRowIterator (Submatrix<Matrix, SubvectorFactory, Trait> *M, typename Matrix::RowIterator pos)
		: _M (M), _pos (pos), _row_valid (false) {}
	SubmatrixRowIterator (const SubmatrixRowIterator &i)
		: _M (i._M), _pos (i._pos), _row (i._row), _row_valid (i._row_valid) {}

	SubmatrixRowIterator &operator = (const SubmatrixRowIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_row = i._row;
		_row_valid = i._row_valid;
		return *this;
	}

	SubmatrixRowIterator &operator ++ () 
	{
		++_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixRowIterator operator ++ (int) 
	{
		SubmatrixRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixRowIterator operator + (difference_type i)
	{
		return SubmatrixRowIterator (_M, _pos + i);
	}

	SubmatrixRowIterator &operator += (difference_type i) 
	{
		_pos += i;
		_row_valid = false;
		return *this;
	}

	SubmatrixRowIterator &operator -- () 
	{
		--_pos;
		_row_valid = false;
		return *this;
	}

	SubmatrixRowIterator operator -- (int) 
	{
		SubmatrixRowIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixRowIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixRowIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixRowIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ update_row (); return _row; }

	pointer operator -> () 
		{ update_row (); return &_row; }

	bool operator == (const SubmatrixRowIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixRowIterator &c) const 
		{ return (_pos != c._pos); }

	operator SubmatrixConstRowIterator<Matrix, SubvectorFactory, Trait> ()
		{ return SubmatrixConstRowIterator<Matrix, SubvectorFactory, Trait> (_M, _pos); }

    private:
	friend class SubmatrixConstRowIterator<Matrix, SubvectorFactory, Trait>;

	Submatrix<Matrix, SubvectorFactory, Trait> *_M;
	typename Matrix::RowIterator _pos;
	typename SubvectorFactory::RowSubvector _row;
	bool _row_valid;

	inline void update_row () {
		if (!_row_valid) {
			_row = (SubvectorFactory ()).MakeRowSubvector (*_M, _pos);
			_row_valid = true;
		}
	}
};

template <class Matrix, class SubvectorFactory, class Trait>
class SubmatrixConstColIterator {
    public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename SubvectorFactory::ColSubvector &reference;
	typedef typename SubvectorFactory::ColSubvector *pointer;
	typedef typename SubvectorFactory::ColSubvector value_type;
	typedef ptrdiff_t difference_type;

	SubmatrixConstColIterator () {}
	SubmatrixConstColIterator (const Submatrix<Matrix, SubvectorFactory, Trait> *M, typename Matrix::ConstColIterator pos)
		: _M (M), _pos (pos), _col_valid (false) {}
	SubmatrixConstColIterator (const SubmatrixConstColIterator &i) : _M (i._M), _pos (i._pos), _col (i._col), _col_valid (i._col_valid) {}

	SubmatrixConstColIterator &operator = (const SubmatrixConstColIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_col = i._col;
		_col_valid = i._col_valid;
		return *this;
	}

	SubmatrixConstColIterator &operator ++ () 
	{
		++_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixConstColIterator operator ++ (int) 
	{
		SubmatrixConstColIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixConstColIterator operator + (difference_type i) const
	{
		return SubmatrixConstColIterator (_M, _pos + i);
	}

	SubmatrixConstColIterator &operator += (difference_type i) 
	{
		_pos += i;
		_col_valid = false;
		return *this;
	}

	SubmatrixConstColIterator &operator -- () 
	{
		--_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixConstColIterator operator -- (int) 
	{
		SubmatrixConstColIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixConstColIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixConstColIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixConstColIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) const
		{ return *(*this + i); }

	const reference operator * () 
		{ update_col (); return _col; }

	const pointer operator -> ()
		{ update_col (); return &_col; }

	bool operator == (const SubmatrixConstColIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixConstColIterator &c) const 
		{ return (_pos != c._pos); }

    private:
	friend class SubmatrixColIterator<Matrix, SubvectorFactory, Trait>;

	const Submatrix<Matrix, SubvectorFactory, Trait> *_M;
	typename Matrix::ConstColIterator _pos;
	typename SubvectorFactory::ConstColSubvector _col;
	bool _col_valid;

	inline void update_col () {
		if (!_col_valid) {
			_col = (SubvectorFactory ()).MakeConstColSubvector (*_M, _pos);
			_col_valid = true;
		}
	}
};

template <class Matrix, class SubvectorFactory, class Trait>
class SubmatrixColIterator {
    public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename SubvectorFactory::ColSubvector &reference;
	typedef typename SubvectorFactory::ColSubvector *pointer;
	typedef typename SubvectorFactory::ColSubvector value_type;
	typedef ptrdiff_t difference_type;

	SubmatrixColIterator () {}
	SubmatrixColIterator (Submatrix<Matrix, SubvectorFactory, Trait> *M, typename Matrix::ColIterator pos)
		: _M (M), _pos (pos), _col_valid (false) {}
	SubmatrixColIterator (const SubmatrixColIterator &i)
		: _M (i._M), _pos (i._pos), _col (i._col), _col_valid (i._col_valid) {}

	SubmatrixColIterator &operator = (const SubmatrixColIterator &i) {
		_M = i._M;
		_pos = i._pos;
		_col = i._col;
		_col_valid = i._col_valid;
		return *this;
	}

	SubmatrixColIterator &operator ++ () 
	{
		++_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixColIterator operator ++ (int) 
	{
		SubmatrixColIterator tmp (*this);
		++*this;
		return tmp;
	}

	SubmatrixColIterator operator + (difference_type i) const
	{
		return SubmatrixColIterator (_M, _pos + i);
	}

	SubmatrixColIterator &operator += (difference_type i) 
	{
		_pos += i;
		_col_valid = false;
		return *this;
	}

	SubmatrixColIterator &operator -- () 
	{
		--_pos;
		_col_valid = false;
		return *this;
	}

	SubmatrixColIterator operator -- (int) 
	{
		SubmatrixColIterator tmp (*this);
		--*this;
		return tmp;
	}

	SubmatrixColIterator operator - (difference_type i) const
		{ return *this + -i; }

	SubmatrixColIterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (SubmatrixColIterator &i) const 
		{ return _pos - i._pos; }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ update_col (); return _col; }

	pointer operator -> () 
		{ update_col (); return &_col; }

	bool operator == (const SubmatrixColIterator &c) const 
		{ return (_pos == c._pos); }

	bool operator != (const SubmatrixColIterator &c) const 
		{ return (_pos != c._pos); }

	operator SubmatrixConstColIterator<Matrix, SubvectorFactory, Trait> ()
		{ return SubmatrixConstColIterator<Matrix, SubvectorFactory, Trait> (_M, _pos); }

    private:
	friend class SubmatrixConstColIterator<Matrix, SubvectorFactory, Trait>;

	Submatrix<Matrix, SubvectorFactory, Trait> *_M;
	typename Matrix::ColIterator _pos;
	typename SubvectorFactory::ColSubvector _col;
	bool _col_valid;

	inline void update_col () {
		if (!_col_valid) {
			_col = (SubvectorFactory ()).MakeColSubvector (*_M, _pos);
			_col_valid = true;
		}
	}
};

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::RowIterator Submatrix<Matrix, SubvectorFactory, Trait>::rowBegin ()
{
	return RowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::RowIterator Submatrix<Matrix, SubvectorFactory, Trait>::rowEnd ()
{
	return RowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::ConstRowIterator Submatrix<Matrix, SubvectorFactory, Trait>::rowBegin () const
{
	return ConstRowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::ConstRowIterator Submatrix<Matrix, SubvectorFactory, Trait>::rowEnd () const
{
	return ConstRowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::ColIterator Submatrix<Matrix, SubvectorFactory, Trait>::colBegin ()
{
	return ColIterator (this, _M->colBegin () + _beg_col);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::ColIterator Submatrix<Matrix, SubvectorFactory, Trait>::colEnd ()
{
	return ColIterator (this, _M->colBegin () + _end_col);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::ConstColIterator Submatrix<Matrix, SubvectorFactory, Trait>::colBegin () const
{
	return ConstColIterator (this, _M->colBegin () + _beg_col);
}

template <class Matrix, class SubvectorFactory, class Trait>
inline typename Submatrix<Matrix, SubvectorFactory, Trait>::ConstColIterator Submatrix<Matrix, SubvectorFactory, Trait>::colEnd () const
{
	return ConstColIterator (this, _M->colBegin () + _end_col);
}

template <class Matrix, class SubvectorFactory>
inline typename Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::RowIterator Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::rowBegin ()
{
	return RowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class SubvectorFactory>
inline typename Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::RowIterator Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::rowEnd ()
{
	return RowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix, class SubvectorFactory>
inline typename Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::ConstRowIterator Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::rowBegin () const
{
	return ConstRowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class SubvectorFactory>
inline typename Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::ConstRowIterator Submatrix<Matrix, SubvectorFactory, MatrixCategories::RowMatrixTag>::rowEnd () const
{
	return ConstRowIterator (this, _M->rowBegin () + _end_row);
}

}

#endif // __MATRIX_SUBMATRIX_INL

