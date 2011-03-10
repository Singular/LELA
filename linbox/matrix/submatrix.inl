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
	struct iterator_traits<LinBox::SubmatrixRowIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::RowSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::RowSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::RowSubvector value_type;
		typedef long difference_type;
	};

	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixConstRowIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstRowSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstRowSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ConstRowSubvector value_type;
		typedef long difference_type;
	};

	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixColIterator<Matrix, Trait> >
	{
		typedef random_access_iterator_tag iterator_category;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ColSubvector &reference;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ColSubvector *pointer;
		typedef typename LinBox::SubvectorFactory<Matrix, Trait>::ColSubvector value_type;
		typedef long difference_type;
	};

	template <class Matrix, class Trait>
	struct iterator_traits<LinBox::SubmatrixConstColIterator<Matrix, Trait> >
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
class SubmatrixConstRowIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixConstRowIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixConstRowIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixConstRowIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixConstRowIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixConstRowIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixConstRowIterator () {}
	SubmatrixConstRowIterator (const Submatrix<Matrix, Trait> *M, typename Matrix::ConstRowIterator pos)
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
	friend class SubmatrixRowIterator<Matrix, Trait>;

	const Submatrix<Matrix, Trait> *_M;
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
class SubmatrixRowIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixRowIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixRowIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixRowIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixRowIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixRowIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixRowIterator () {}
	SubmatrixRowIterator (Submatrix<Matrix, Trait> *M, typename Matrix::RowIterator pos)
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

	SubmatrixRowIterator operator + (difference_type i) const
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

	operator SubmatrixConstRowIterator<Matrix, Trait> ()
		{ return SubmatrixConstRowIterator<Matrix, Trait> (_M, _pos); }

    private:
	friend class SubmatrixConstRowIterator<Matrix, Trait>;

	Submatrix<Matrix, Trait> *_M;
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
class SubmatrixConstColIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixConstColIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixConstColIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixConstColIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixConstColIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixConstColIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixConstColIterator () {}
	SubmatrixConstColIterator (const Submatrix<Matrix, Trait> *M, typename Matrix::ConstColIterator pos)
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
	friend class SubmatrixColIterator<Matrix, Trait>;

	const Submatrix<Matrix, Trait> *_M;
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
class SubmatrixColIterator {
    public:
	typedef typename std::iterator_traits<SubmatrixColIterator<Matrix, Trait> >::iterator_category iterator_category;
	typedef typename std::iterator_traits<SubmatrixColIterator<Matrix, Trait> >::reference reference;
	typedef typename std::iterator_traits<SubmatrixColIterator<Matrix, Trait> >::pointer pointer;
	typedef typename std::iterator_traits<SubmatrixColIterator<Matrix, Trait> >::value_type value_type;
	typedef typename std::iterator_traits<SubmatrixColIterator<Matrix, Trait> >::difference_type difference_type;

	SubmatrixColIterator () {}
	SubmatrixColIterator (Submatrix<Matrix, Trait> *M, typename Matrix::ColIterator pos)
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

	operator SubmatrixConstColIterator<Matrix, Trait> ()
		{ return SubmatrixConstColIterator<Matrix, Trait> (_M, _pos); }

    private:
	friend class SubmatrixConstColIterator<Matrix, Trait>;

	Submatrix<Matrix, Trait> *_M;
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
inline typename Submatrix<Matrix, Trait>::RowIterator Submatrix<Matrix, Trait>::rowBegin ()
{
	return RowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class Trait>
inline typename Submatrix<Matrix, Trait>::RowIterator Submatrix<Matrix, Trait>::rowEnd ()
{
	return RowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix, class Trait>
inline typename Submatrix<Matrix, Trait>::ConstRowIterator Submatrix<Matrix, Trait>::rowBegin () const
{
	return ConstRowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix, class Trait>
inline typename Submatrix<Matrix, Trait>::ConstRowIterator Submatrix<Matrix, Trait>::rowEnd () const
{
	return ConstRowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix>
inline typename Submatrix<Matrix, MatrixCategories::RowMatrixTag>::RowIterator Submatrix<Matrix, MatrixCategories::RowMatrixTag>::rowBegin ()
{
	return RowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix>
inline typename Submatrix<Matrix, MatrixCategories::RowMatrixTag>::RowIterator Submatrix<Matrix, MatrixCategories::RowMatrixTag>::rowEnd ()
{
	return RowIterator (this, _M->rowBegin () + _end_row);
}

template <class Matrix>
inline typename Submatrix<Matrix, MatrixCategories::RowMatrixTag>::ConstRowIterator Submatrix<Matrix, MatrixCategories::RowMatrixTag>::rowBegin () const
{
	return ConstRowIterator (this, _M->rowBegin () + _beg_row);
}

template <class Matrix>
inline typename Submatrix<Matrix, MatrixCategories::RowMatrixTag>::ConstRowIterator Submatrix<Matrix, MatrixCategories::RowMatrixTag>::rowEnd () const
{
	return ConstRowIterator (this, _M->rowBegin () + _end_row);
}

}

#endif // __MATRIX_SUBMATRIX_INL

