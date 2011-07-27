/* lela/matrix/raw-iterator.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>,
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_MATRIX_RAW_ITERATOR_H
#define __LELA_MATRIX_RAW_ITERATOR_H

#include <utility>
#include <vector>

namespace LELA
{

/** Generic raw iterator
 *
 * Forward iterator which iterates through all entries in the
 * matrix. Values are entries.
 *
 * \ingroup matrix
 */
template <class Iterator, class VectorTrait>
class MatrixRawIterator;

/** Generic raw indexed iterator
 *
 * Forward iterator which iterates through all entries in the
 * matrix. Values are STL-pairs of row- and column-index of
 * corresponding entry.
 *
 * \ingroup matrix
 */
template <class Iterator, class VectorTrait, bool Direction>
class MatrixRawIndexedIterator;

template <class Iterator>
class MatrixRawIterator<Iterator, VectorRepresentationTypes::Dense>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef typename std::iterator_traits<typename Vector::iterator>::reference reference;
	typedef const reference const_reference;
	typedef typename std::iterator_traits<typename Vector::iterator>::value_type value_type;
	typedef typename std::iterator_traits<typename Vector::iterator>::difference_type difference_type;

	MatrixRawIterator (Iterator rowcol, size_t pos, Iterator rowcol_end, size_t rowcol_len)
		: _rowcol (rowcol), _pos (rowcol->begin ()) {}

	MatrixRawIterator () {}

	MatrixRawIterator &operator = (const MatrixRawIterator &p)
	{
		_rowcol = p._rowcol;
		_pos = p._pos;
		return *this;
	}
    
	MatrixRawIterator &operator ++ ()
	{
		++_pos;

		if (_pos == _rowcol->end ()) {
			++_rowcol;
			_pos = _rowcol->begin ();
		}

		return *this;
	}
    
	MatrixRawIterator operator ++ (int)
	{
		MatrixRawIterator tmp (*this);
		++*this;
		return tmp;
	}

	reference operator * ()
		{ return *_pos; }
 
	bool operator == (const MatrixRawIterator &c) const
		{ return (_rowcol == c._rowcol) && (_pos == c._pos); }

	bool operator != (const MatrixRawIterator &c) const
		{ return (_rowcol != c._rowcol) || (_pos != c._pos); }

    private:
	Iterator _rowcol;
	typename Vector::iterator _pos;
};

template <class Iterator>
class MatrixRawIterator<Iterator, VectorRepresentationTypes::Sparse>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type::second_type value_type;
	typedef typename std::iterator_traits<typename Vector::const_iterator>::difference_type difference_type;
	typedef value_type reference;
	typedef const reference const_reference;

	MatrixRawIterator (Iterator rowcol, size_t pos, Iterator rowcol_end, size_t rowcol_len)
		: _rowcol (rowcol), _rowcol_end (rowcol_end), _iter_valid (false)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ())
			++_rowcol;
	}

	MatrixRawIterator () {}

	MatrixRawIterator &operator = (const MatrixRawIterator &p)
	{
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_pos = p._pos;
		_iter_valid = p._iter_valid;
		return *this;
	}
	
	MatrixRawIterator &operator ++ ()
	{
		++_pos;

		if (_pos == _rowcol->end ()) {
			do {
				++_rowcol;
			} while (_rowcol != _rowcol_end && _rowcol->empty ());

			_iter_valid = false;
		}

		return *this;
	}
    
	MatrixRawIterator operator ++ (int)
	{
		MatrixRawIterator tmp (*this);
		++*this;
		return tmp;
	}

	const_reference operator * ()
		{ update_iter (); return (const_reference) _pos->second; }
 
	bool operator == (const MatrixRawIterator &c) const
		{ return (_rowcol == c._rowcol) && (!(_iter_valid || c._iter_valid) || ((_iter_valid && c._iter_valid) && (_pos == c._pos))); }

	bool operator != (const MatrixRawIterator &c) const
		{ return (_rowcol != c._rowcol) || ((_iter_valid || c._iter_valid) && (!(_iter_valid && c._iter_valid) || (_pos != c._pos))); }

    private:
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _pos;
	bool _iter_valid;

	void update_iter ()
	{
		if (!_iter_valid) {
			_pos = _rowcol->begin ();
			_iter_valid = true;
		}
	}
};

template <class Iterator>
class MatrixRawIterator<Iterator, VectorRepresentationTypes::Dense01>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef typename std::iterator_traits<typename Vector::const_iterator>::reference reference;
	typedef typename Vector::const_iterator::const_reference const_reference;
	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type value_type;
	typedef typename std::iterator_traits<typename Vector::const_iterator>::difference_type difference_type;

	MatrixRawIterator (Iterator rowcol, size_t pos, Iterator rowcol_end, size_t rowcol_len)
		: _rowcol (rowcol), _pos (rowcol->begin ()) {}

	MatrixRawIterator () {}

	MatrixRawIterator &operator = (const MatrixRawIterator &p)
	{
		_rowcol = p._rowcol;
		_pos = p._pos;
		return *this;
	}
    
	MatrixRawIterator &operator ++ ()
	{
		++_pos;
		return *this;
	}
    
	MatrixRawIterator operator ++ (int)
	{
		MatrixRawIterator tmp (*this);
		++*this;
		return tmp;
	}

	reference operator * ()
	{
		if (_pos == _rowcol->end ()) {
			++_rowcol;
			_pos = _rowcol->begin ();
		}

		return *_pos;
	}
 
	bool operator == (const MatrixRawIterator &c) const
		{ return (_rowcol == c._rowcol) && (_pos == c._pos); }

	bool operator != (const MatrixRawIterator &c) const
		{ return (_rowcol != c._rowcol) || (_pos != c._pos); }

    private:
	Iterator _rowcol;
	typename Vector::const_iterator _pos;
};

template <class Iterator>
class MatrixRawIterator<Iterator, VectorRepresentationTypes::Sparse01>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef bool const_reference;
	typedef const_reference reference;
	typedef bool value_type;
	typedef typename std::iterator_traits<typename Vector::const_iterator>::difference_type difference_type;

	MatrixRawIterator (Iterator rowcol, size_t pos, Iterator rowcol_end, size_t rowcol_len)
		: _rowcol (rowcol), _rowcol_end (rowcol_end), _iter_valid (false)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ())
			++_rowcol;
	}

	MatrixRawIterator () {}

	MatrixRawIterator &operator = (const MatrixRawIterator &p)
	{
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_pos = p._pos;
		_iter_valid = p._iter_valid;
		return *this;
	}
    
	MatrixRawIterator &operator ++ ()
	{
		++_pos;

		if (_pos == _rowcol->end ()) {
			do {
				++_rowcol;
			} while (_rowcol != _rowcol_end && _rowcol->empty ());

			_iter_valid = false;
		}

		return *this;
	}
    
	MatrixRawIterator operator ++ (int)
	{
		MatrixRawIterator tmp (*this);
		++*this;
		return tmp;
	}

	const_reference operator * ()
		{ update_iter (); return true; }
 
	bool operator == (const MatrixRawIterator &c) const
		{ return (_rowcol == c._rowcol) && (!(_iter_valid || c._iter_valid) || ((_iter_valid && c._iter_valid) && (_pos == c._pos))); }

	bool operator != (const MatrixRawIterator &c) const
		{ return (_rowcol != c._rowcol) || ((_iter_valid || c._iter_valid) && (!(_iter_valid && c._iter_valid) || (_pos != c._pos))); }

    private:
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _pos;
	bool _iter_valid;

	void update_iter ()
	{
		if (!_iter_valid) {
			_pos = _rowcol->begin ();
			_iter_valid = true;
		}
	}
};

template <class Iterator>
class MatrixRawIterator<Iterator, VectorRepresentationTypes::Hybrid01>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef typename Vector::Endianness Endianness;
	typedef bool value_type;
	typedef value_type reference;
	typedef value_type const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _t (Endianness::e_0), _rowcol_len (rowcol_len)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.first;
		}

		if (_rowcol != _rowcol_end) {
			_iter = _rowcol->begin ();
			_pos.second = _iter->first << WordTraits<typename Vector::word_type>::logof_size;
		}
	}

	MatrixRawIterator () {}

	MatrixRawIterator &operator = (const MatrixRawIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_t = p._t;
		_rowcol_len = p._rowcol_len;
		return *this;
	}
    
	MatrixRawIterator &operator ++ ()
	{
		_t = Endianness::shift_right (_t, 1);
		++_pos.second;

		if (_pos.second >= _rowcol_len) {
			_t = Endianness::e_0;
			advance_rowcol ();
		}
		else if (_t == 0) {
			_t = Endianness::e_0;
			++_iter;

			if (_iter == _rowcol->end ())
				advance_rowcol ();
			else
				_pos.second = _iter->first << WordTraits<typename Vector::word_type>::logof_size;
		}

		return *this;
	}
    
	MatrixRawIterator operator ++ (int)
	{
		MatrixRawIterator tmp (*this);
		++*this;
		return tmp;
	}

	value_type operator * ()
		{ return _iter->second & _t; }
 
	bool operator == (const MatrixRawIterator &c) const
		{ return (_rowcol == c._rowcol) && ((_rowcol == _rowcol_end) || ((_iter == c._iter) && (_t == c._t))); }

	bool operator != (const MatrixRawIterator &c)
		{ return (_rowcol != c._rowcol) || ((_rowcol != _rowcol_end) && ((_iter != c._iter) || (_t != c._t))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	bool _iter_valid;
	typename Vector::word_type _t;
	size_t _rowcol_len;

	void advance_rowcol ()
	{
		do {
			++_rowcol;
			++_pos.first;
		} while (_rowcol != _rowcol_end && _rowcol->empty ());

		if (_rowcol != _rowcol_end) {
			_iter = _rowcol->begin ();
			_pos.second = _iter->first << WordTraits<typename Vector::word_type>::logof_size;
		}
	}
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Dense, false>
{
    public:
	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol) {}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator &operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		return *this;
	}
    
	MatrixRawIndexedIterator &operator ++ ()
	{
		++_pos.second;

		if (_pos.second == _rowcol->size ()) {
			++_pos.first;
			_pos.second = 0;
		}

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> () const
		{ return &_pos; }

	value_type operator * () const
		{ return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c) const
		{ return (_pos == c._pos); }

	bool operator != (const MatrixRawIndexedIterator &c) const
		{ return (_pos != c._pos); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol;
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Sparse, false>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _iter_valid (false)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.first;
		}
	}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator& operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_iter_valid = p._iter_valid;
		return *this;
	}
    
	MatrixRawIndexedIterator& operator ++ ()
	{
		++_iter;

		if (_iter == _rowcol->end ()) {
			do {
				++_rowcol;
				++_pos.first;
			} while (_rowcol != _rowcol_end && _rowcol->empty ());

			_iter_valid = false;
		} else
			_pos.second = _iter->first;

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> ()
		{ update_iter (); return &_pos; }

	value_type operator * () const
		{ update_iter (); return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol == c._rowcol) && (!(_iter_valid || c._iter_valid) || ((_iter_valid && c._iter_valid) && (_iter == c._iter))); }

	bool operator != (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol != c._rowcol) || ((_iter_valid || c._iter_valid) && (!(_iter_valid && c._iter_valid) || (_iter != c._iter))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	bool _iter_valid;

	void update_iter ()
	{
		if (!_iter_valid) {
			_iter = _rowcol->begin ();
			_pos.second = _iter->first;
			_iter_valid = true;
		}
	}
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Sparse, true>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _iter_valid (false)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.second;
		}
	}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator& operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_iter_valid = p._iter_valid;
		return *this;
	}
    
	MatrixRawIndexedIterator& operator ++ ()
	{
		++_iter;

		if (_iter == _rowcol->end ()) {
			do {
				++_rowcol;
				++_pos.second;
			} while (_rowcol != _rowcol_end && _rowcol->empty ());

			_iter_valid = false;
		} else
			_pos.first = _iter->first;

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> ()
		{ update_iter (); return &_pos; }

	value_type operator * () const
		{ update_iter (); return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol == c._rowcol) && (!(_iter_valid || c._iter_valid) || ((_iter_valid && c._iter_valid) && (_iter == c._iter))); }

	bool operator != (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol != c._rowcol) || ((_iter_valid || c._iter_valid) && (!(_iter_valid && c._iter_valid) || (_iter != c._iter))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	bool _iter_valid;

	void update_iter ()
	{
		if (!_iter_valid) {
			_iter = _rowcol->begin ();
			_pos.first = _iter->first;
			_iter_valid = true;
		}
	}
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Dense01, false>
{
    public:
	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol) {}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator& operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		return *this;
	}
    
	MatrixRawIndexedIterator& operator ++ ()
	{
		++_pos.second;

		if (_pos.second == _rowcol->size ()) {
			++_pos.first;
			_pos.second = 0;
		}

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> () const
		{ return &_pos; }

	value_type operator * () const
		{ return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c) const
		{ return (_pos == c._pos); }

	bool operator != (const MatrixRawIndexedIterator &c) const
		{ return (_pos != c._pos); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol;
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Sparse01, false>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _iter_valid (false)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.first;
		}
	}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator& operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_iter_valid = p._iter_valid;
		return *this;
	}
    
	MatrixRawIndexedIterator& operator ++ ()
	{
		++_iter;

		if (_iter == _rowcol->end ()) {
			do {
				++_rowcol;
				++_pos.first;
			} while (_rowcol != _rowcol_end && _rowcol->empty ());

			_iter_valid = false;
		} else
			_pos.second = *_iter;

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> ()
		{ update_iter (); return &_pos; }

	value_type operator * () const
		{ update_iter (); return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol == c._rowcol) && (!(_iter_valid || c._iter_valid) || ((_iter_valid && c._iter_valid) && (_iter == c._iter))); }

	bool operator != (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol != c._rowcol) || ((_iter_valid || c._iter_valid) && (!(_iter_valid && c._iter_valid) || (_iter != c._iter))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	bool _iter_valid;

	void update_iter ()
	{
		if (!_iter_valid) {
			_iter = _rowcol->begin ();
			_pos.second = *_iter;
			_iter_valid = true;
		}
	}
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Sparse01, true>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _iter_valid (false)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.second;
		}
	}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator& operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_iter_valid = p._iter_valid;
		return *this;
	}
    
	MatrixRawIndexedIterator& operator ++ ()
	{
		++_iter;

		if (_iter == _rowcol->end ()) {
			do {
				++_rowcol;
				++_pos.second;
			} while (_rowcol != _rowcol_end && _rowcol->empty ());

			_iter_valid = false;
		} else
			_pos.first = *_iter;

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> ()
		{ update_iter (); return &_pos; }

	value_type operator * () const
		{ update_iter (); return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol == c._rowcol) && (!(_iter_valid || c._iter_valid) || ((_iter_valid && c._iter_valid) && (_iter == c._iter))); }

	bool operator != (const MatrixRawIndexedIterator &c) const
		{ return (_rowcol != c._rowcol) || ((_iter_valid || c._iter_valid) && (!(_iter_valid && c._iter_valid) || (_iter != c._iter))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	bool _iter_valid;

	void update_iter ()
	{
		if (!_iter_valid) {
			_iter = _rowcol->begin ();
			_pos.first = *_iter;
			_iter_valid = true;
		}
	}
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Hybrid01, false>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef typename Vector::Endianness Endianness;
	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _t (Endianness::e_0), _rowcol_len (rowcol_len)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.first;
		}

		if (_rowcol != _rowcol_end) {
			_iter = _rowcol->begin ();
			_pos.second = (_iter->first << WordTraits<typename Vector::word_type>::logof_size);
		}
	}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator &operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_t = p._t;
		_rowcol_len = p._rowcol_len;
		return *this;
	}
    
	MatrixRawIndexedIterator &operator ++ ()
	{
		_t = Endianness::shift_right (_t, 1);
		++_pos.second;

		if (_pos.second >= _rowcol_len) {
			_t = Endianness::e_0;
			advance_rowcol ();
		}
		else if (_t == 0) {
			_t = Endianness::e_0;
			++_iter;

			if (_iter == _rowcol->end ()) {
				advance_rowcol ();
			} else
				_pos.second = (_iter->first << WordTraits<typename Vector::word_type>::logof_size);
		}

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> ()
		{ return &_pos; }

	value_type operator * () const
		{ return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c)
		{ return (_rowcol == c._rowcol) && ((_rowcol == _rowcol_end) || ((_iter == c._iter) && (_t == c._t))); }

	bool operator != (const MatrixRawIndexedIterator &c)
		{ return (_rowcol != c._rowcol) || ((_rowcol != _rowcol_end) && ((_iter != c._iter) || (_t != c._t))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	typename Vector::word_type _t;
	size_t _rowcol_len;

	void advance_rowcol ()
	{
		++_rowcol;
		++_pos.first;

		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.first;
		}

		if (_rowcol != _rowcol_end) {
			_iter = _rowcol->begin ();
			_pos.second = (_iter->first << WordTraits<typename Vector::word_type>::logof_size);
		}
	}
};

template <class Iterator>
class MatrixRawIndexedIterator<Iterator, VectorRepresentationTypes::Hybrid01, true>
{
    public:
	typedef typename Iterator::value_type Vector;

	typedef typename Vector::Endianness Endianness;
	typedef std::pair<size_t, size_t> value_type;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef typename Iterator::difference_type difference_type;

	MatrixRawIndexedIterator (Iterator rowcol, size_t i, Iterator rowcol_end, size_t rowcol_len)
		: _pos (i, 0), _rowcol (rowcol), _rowcol_end (rowcol_end), _t (Endianness::e_0), _rowcol_len (rowcol_len)
	{
		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.second;
		}

		if (_rowcol != _rowcol_end) {
			_iter = _rowcol->begin ();
			_pos.first = (_iter->first << WordTraits<typename Vector::word_type>::logof_size);
		}
	}

	MatrixRawIndexedIterator () {}

	MatrixRawIndexedIterator &operator = (const MatrixRawIndexedIterator &p)
	{
		_pos = p._pos;
		_rowcol = p._rowcol;
		_rowcol_end = p._rowcol_end;
		_iter = p._iter;
		_t = p._t;
		_rowcol_len = p._rowcol_len;
		return *this;
	}
    
	MatrixRawIndexedIterator &operator ++ ()
	{
		_t = Endianness::shift_right (_t, 1);
		++_pos.first;

		if (_pos.first >= _rowcol_len) {
			_t = Endianness::e_0;
			advance_rowcol ();
		}
		else if (_t == 0) {
			_t = Endianness::e_0;
			++_iter;

			if (_iter == _rowcol->end ()) {
				advance_rowcol ();
			} else
				_pos.first = (_iter->first << WordTraits<typename Vector::word_type>::logof_size);
		}

		return *this;
	}
    
	MatrixRawIndexedIterator operator ++ (int)
	{
		MatrixRawIndexedIterator tmp (*this);
		++*this;
		return tmp;
	}

	const value_type *operator -> ()
		{ return &_pos; }

	value_type operator * () const
		{ return _pos; }
 
	bool operator == (const MatrixRawIndexedIterator &c)
		{ return (_rowcol == c._rowcol) && ((_rowcol == _rowcol_end) || ((_iter == c._iter) && (_t == c._t))); }

	bool operator != (const MatrixRawIndexedIterator &c)
		{ return (_rowcol != c._rowcol) || ((_rowcol != _rowcol_end) && ((_iter != c._iter) || (_t != c._t))); }

    private:
	std::pair<size_t, size_t> _pos;
	Iterator _rowcol, _rowcol_end;
	typename Vector::const_iterator _iter;
	typename Vector::word_type _t;
	size_t _rowcol_len;

	void advance_rowcol ()
	{
		++_rowcol;
		++_pos.second;

		while (_rowcol != _rowcol_end && _rowcol->empty ()) {
			++_rowcol;
			++_pos.second;
		}

		if (_rowcol != _rowcol_end) {
			_iter = _rowcol->begin ();
			_pos.first = (_iter->first << WordTraits<typename Vector::word_type>::logof_size);
		}
	}
};

} // namespace LELA

namespace std {
	template <class Iterator, class VectorTrait>
	struct iterator_traits<LELA::MatrixRawIterator<Iterator, VectorTrait> >
	{
		typedef forward_iterator_tag iterator_category;
		typedef typename LELA::MatrixRawIterator<Iterator, VectorTrait>::reference &reference;
		typedef typename LELA::MatrixRawIterator<Iterator, VectorTrait>::value_type *pointer;
		typedef typename LELA::MatrixRawIterator<Iterator, VectorTrait>::value_type value_type;
		typedef typename LELA::MatrixRawIterator<Iterator, VectorTrait>::difference_type difference_type;
	};

	template <class Iterator, class VectorTrait, bool Direction>
	struct iterator_traits<LELA::MatrixRawIndexedIterator<Iterator, VectorTrait, Direction> >
	{
		typedef forward_iterator_tag iterator_category;
		typedef typename LELA::MatrixRawIndexedIterator<Iterator, VectorTrait, Direction>::reference &reference;
		typedef typename LELA::MatrixRawIndexedIterator<Iterator, VectorTrait, Direction>::value_type *pointer;
		typedef typename LELA::MatrixRawIndexedIterator<Iterator, VectorTrait, Direction>::value_type value_type;
		typedef typename LELA::MatrixRawIndexedIterator<Iterator, VectorTrait, Direction>::difference_type difference_type;
	};
}

#endif // __LELA_MATRIX_RAW_ITERATOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
