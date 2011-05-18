/* linbox/vector/bit-subvector.h
 * Copyright 2010 Bradford Hovinen
 *
 * Evolved from bit-vector.h
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __VECTOR_BIT_SUBVECTOR_H
#define __VECTOR_BIT_SUBVECTOR_H

#include <iterator>
#include <vector>
#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"

namespace LinBox
{

/** A subvector of a BitVector. The interface is identical to that of \ref BitVector.
 \ingroup vector
 */

template <class Iterator, class ConstIterator = Iterator>
class BitSubvector
{
    public:
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
	typedef typename std::iterator_traits<Iterator>::pointer	 pointer;
	typedef typename std::iterator_traits<Iterator>::reference	 reference;
	typedef Iterator iterator;
	typedef ConstIterator const_iterator;
	typedef VectorCategories::DenseZeroOneVectorTag VectorCategory; 

	typedef typename Iterator::size_type	 size_type;
	typedef typename Iterator::const_reference const_reference;

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	class word_reference;
	class word_iterator;
	class const_word_iterator;

	typedef std::reverse_iterator<word_iterator> reverse_word_iterator;
	typedef std::reverse_iterator<const_word_iterator> const_reverse_word_iterator;

	friend class word_reference;
	friend class word_iterator;
	friend class const_word_iterator;

	typedef typename Iterator::Endianness Endianness;
	typedef typename std::iterator_traits<word_iterator>::value_type word_type;

	BitSubvector () {}

	template <class Endianness>
	BitSubvector (BitVector<Endianness> &v, size_t start, size_t end)
		: _begin (v.begin () + start), _end (v.begin () + end) { set_end_word (); }

	BitSubvector (BitSubvector &v, size_t start, size_t end)
		: _begin (v._begin + start), _end (v._begin + end) { set_end_word (); }

	BitSubvector (iterator begin, iterator end)
		: _begin (begin), _end (end) { set_end_word (); }
		
	BitSubvector (const BitSubvector &v) 
		: _begin (v._begin), _end (v._end), _end_word (v._end_word), _end_marker (v._end_marker) {}

	~BitSubvector () {}

	inline iterator                    begin      (void)       { return _begin; }
	inline const_iterator              begin      (void) const { return const_iterator (_begin); }
	inline iterator                    end        (void)       { return _end; }
	inline const_iterator              end        (void) const { return const_iterator (_end); }

	inline reverse_iterator            rbegin     (void)       { return reverse_iterator (_end); }
	inline const_reverse_iterator      rbegin     (void) const { return const_reverse_iterator (_end); }
	inline reverse_iterator            rend       (void)       { return reverse_iterator (_begin); }
	inline const_reverse_iterator      rend       (void) const { return const_reverse_iterator (_begin); }

	inline word_iterator               wordBegin  (void)       { return word_iterator (*this, _begin.word ()); }
	inline const_word_iterator         wordBegin  (void) const { return const_word_iterator (*this, _begin.word ()); }
	inline word_iterator               wordEnd    (void)       { return word_iterator (*this, _end_word); }
	inline const_word_iterator         wordEnd    (void) const { return const_word_iterator (*this, _end_word); }

	inline reverse_word_iterator       wordRbegin (void)       { return reverse_word_iterator (wordEnd ()); }
	inline const_reverse_word_iterator wordRbegin (void) const { return const_reverse_word_iterator (wordEnd ()); }
	inline reverse_word_iterator       wordRend   (void)       { return reverse_word_iterator (wordBegin ()); }
	inline const_reverse_word_iterator wordRend   (void) const { return const_reverse_word_iterator (wordBegin ()); }

	// Element access

	inline reference       operator[] (size_type n)       { return *(_begin + n); }
	inline const_reference operator[] (size_type n) const { return *(_begin + n); }

	inline reference at (size_type n)
	{   
		iterator p = _begin + n;
		if ( _begin <= p && p < _end ) 
			return *p;
		else 
			throw std::out_of_range (std::string ("n"));
	}

	inline const_reference at (size_type n) const
	{
		const_iterator p = _begin + n;
		if ( _begin <= p && p < _end)
			return *p;
		else 
			throw std::out_of_range (std::string ("n"));
	}

	BitSubvector &operator = (const BitSubvector &sub)
		{ _begin = sub._begin; _end = sub._end; _end_word = sub._end_word; _end_marker = sub._end_marker; return *this; }

	inline reference       front (void)       { return *_begin; }
	inline const_reference front (void) const { return *_begin; }
	inline reference       back  (void)       { return *(_end - 1); }
	inline const_reference back  (void) const { return *(_end - 1); }

	inline size_type size      (void) const { return _end - _begin;    }
	inline bool      empty     (void) const { return _end == _begin;   }
	inline size_type max_size  (void) const { return _end - _begin; }

	inline bool operator == (const BitSubvector &v) const
		{ return (_begin == v._begin) && (_end == v._end); }
	inline bool operator != (const BitSubvector &v) const
		{ return (_begin != v._begin) || (_end != v._end); }

    protected:

	Iterator     _begin;
	Iterator     _end;

	typename Iterator::word_iterator _end_word, _end_marker;

	void set_end_word ()
	{
		_end_marker = _end_word = _end.word ();

		if (_end.pos () > _begin.pos ())
			++_end_word;
		if (_end.pos () > 0)
			++_end_marker;
	}

    public:

	class word_reference
	{
	public:
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type word;
		typedef typename Iterator::Endianness Endianness;

		word_reference () {}

		word_reference (typename Iterator::word_iterator position, uint8 shift)
			: _pos (position), _shift (shift)
			{}

		~word_reference () {}

		word_reference &operator = (word_reference &a)
		{
			*this = (word) a;
			return a;
		}

		word_reference &operator = (word v)
		{
			v &= _mask;

			*_pos = (*_pos & (Endianness::shift_right (~_mask, _shift) | Endianness::mask_left (_shift))) | Endianness::shift_right (v, _shift);

			if (_shift != 0 && !_just_this_word)
				_pos[1] = (_pos[1] & (~Endianness::shift_left (_mask, WordTraits<word>::bits - _shift) | Endianness::mask_right (_shift))) |
					Endianness::shift_left (v, WordTraits<word>::bits - _shift);

			return *this;
		}

		word_reference &operator &= (word_reference &a)
			{ return *this = (word) *this & (word) a; }

		word_reference &operator &= (word v) 
			{ return *this = (word) *this & v; }

		word_reference &operator |= (word_reference &a) 
			{ return *this = (word) *this | (word) a; }

		word_reference &operator |= (word v) 
			{ return *this = (word) *this | v; }

		word_reference &operator ^= (word_reference &a) 
			{ return *this = ((word) *this) ^ ((word) a); }

		word_reference &operator ^= (word v) 
			{ return *this = ((word) *this) ^ v; }

		operator word (void) const
			{ return (Endianness::shift_left (*_pos, _shift) | ((_just_this_word || _shift == 0) ? 0 : Endianness::shift_right (_pos[1], WordTraits<word>::bits - _shift))) & _mask; }

	private:
		friend class word_iterator;
		friend class const_word_iterator;
		friend class const_word_reference;

		typename Iterator::word_iterator _pos;
		uint8 _shift;

		word _mask;
		bool _just_this_word;
	};

	class word_iterator
	{
	public:
		typedef std::random_access_iterator_tag iterator_category;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference reference;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference *pointer;
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::value_type value_type;
		typedef typename std::iterator_traits<typename Iterator::word_iterator>::difference_type difference_type;
		typedef typename Iterator::Endianness Endianness;

		word_iterator () {}
		word_iterator (BitSubvector<Iterator, ConstIterator> &v, typename Iterator::word_iterator pos)
			: _v (&v), _ref (pos, v._begin.pos ()) { init_mask (); }
		word_iterator (const word_iterator &i)
			: _v (i._v), _ref (i._ref._pos, i._ref._shift) { init_mask (); }

		word_iterator &operator = (const word_iterator &i)
		{
			_v = i._v;
			_ref._pos = i._ref._pos;
			_ref._shift = i._ref._shift;
			_ref._mask = i._ref._mask;
			_ref._just_this_word = i._ref._just_this_word;
			return *this;
		}

		word_iterator &operator ++ () 
		{
			++_ref._pos;
			update_mask ();
			return *this;
		}

		word_iterator operator ++ (int) 
		{
			word_iterator tmp (*this);
			++*this;
			return tmp;
		}

		word_iterator operator + (difference_type i) const
			{ return word_iterator (*const_cast<BitSubvector<Iterator, ConstIterator> *> (_v), _ref._pos + i); }

		word_iterator &operator += (difference_type i) 
		{
			_ref._pos += i;
			init_mask ();
			return *this;
		}

		word_iterator &operator -- () 
		{
			--_ref._pos;
			init_mask ();
			return *this;
		}

		word_iterator operator -- (int) 
		{
			word_iterator tmp (*this);
			--*this;
			return tmp;
		}

		word_iterator operator - (difference_type i) const
			{ return *this + -i; }

		word_iterator &operator -= (difference_type i) 
			{ return *this += -i; }

		difference_type operator - (word_iterator &i) const 
			{ return _ref._pos - i._ref._pos; }

		typename BitSubvector<Iterator, ConstIterator>::word_reference operator [] (long i) 
			{ return *(*this + i); }

		typename BitSubvector<Iterator, ConstIterator>::word_reference operator * () 
			{ return _ref; }

		typename BitSubvector<Iterator, ConstIterator>::word_reference *operator -> () 
			{ return &_ref; }

		value_type operator * () const 
			{ return _ref; }

		bool operator == (const word_iterator &c) const 
			{ return (_ref._pos == c._ref._pos); }

		bool operator != (const word_iterator &c) const 
			{ return (_ref._pos != c._ref._pos); }

	private:
		friend class const_word_iterator;

		BitSubvector<Iterator, ConstIterator> *_v;
		word_reference _ref;

		inline void init_mask ()
		{
			_ref._mask = (value_type) 0 - (value_type) 1;
			update_mask ();
		}

		inline void update_mask ()
		{
			if (((_ref._pos - _v->_begin.word () + 1) << WordTraits<value_type>::logof_size) > _v->_end - _v->_begin)
				_ref._mask = Endianness::mask_left ((_v->_end - _v->_begin) & WordTraits<value_type>::pos_mask);

			if (_ref._pos + 1 == _v->_end_marker)
				_ref._just_this_word = true;
			else
				_ref._just_this_word = false;
		}
	};

	class const_word_iterator
	{
	public:
		typedef std::random_access_iterator_tag iterator_category;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference reference;
		typedef typename BitSubvector<Iterator, ConstIterator>::word_reference *pointer;
		typedef typename std::iterator_traits<typename Iterator::const_word_iterator>::value_type value_type;
		typedef typename std::iterator_traits<typename Iterator::const_word_iterator>::difference_type difference_type;
		typedef typename Iterator::Endianness Endianness;

		const_word_iterator () {}
		const_word_iterator (const BitSubvector<Iterator, ConstIterator> &v, typename Iterator::const_word_iterator pos)
			: _v (&v), _pos (pos), _ref_valid (false) {}
		const_word_iterator (const const_word_iterator  &i)
			: _v (i._v), _pos (i._pos), _word (i._word), _ref_valid (i._ref_valid) {}
		const_word_iterator (const word_iterator  &i)
			: _v (i._v), _pos (i._ref._pos), _ref_valid (false) {}

		const_word_iterator &operator = (const const_word_iterator  &i)
		{
			_v = i._v;
			_pos = i._pos;
			_word = i._word;
			_ref_valid = i._ref_valid;
			return *this;
		}

		const_word_iterator &operator = (const word_iterator  &i)
		{
			_v = i._v;
			_pos = i._ref._pos;
			_ref_valid = false;
			return *this;
		}

		const_word_iterator &operator ++ () 
		{
			++_pos;
			_ref_valid = false;
			return *this;
		}

		const_word_iterator operator ++ (int) 
		{
			const_word_iterator tmp (*this);
			++*this;
			return tmp;
		}

		const_word_iterator operator + (difference_type i) const
			{ return const_word_iterator (*_v, _pos + i); }

		const_word_iterator &operator += (difference_type i) 
		{
			_pos += i;
			_ref_valid = false;
			return *this;
		}

		const_word_iterator &operator -- () 
		{
			--_pos;
			_ref_valid = false;
			return *this;
		}

		const_word_iterator operator -- (int) 
		{
			const_word_iterator tmp (*this);
			--*this;
			return tmp;
		}

		const_word_iterator operator - (difference_type i) const
			{ return *this + -i; }

		const_word_iterator &operator -= (difference_type i) 
			{ return *this += -i; }

		difference_type operator - (const_word_iterator &i) const 
			{ return _pos - i._pos; }

		value_type operator [] (long i) 
			{ return *(*this + i); }

		value_type operator * ()
			{ update_ref (); return _word; }

		bool operator == (const word_iterator &c) const 
			{ return (_pos == c._ref._pos); }

		bool operator == (const const_word_iterator &c) const 
			{ return (_pos == c._pos); }

		bool operator != (const word_iterator &c) const 
			{ return (_pos != c._ref._pos); }

		bool operator != (const const_word_iterator &c) const 
			{ return (_pos != c._pos); }

	private:
		friend class word_iterator;

		const BitSubvector<Iterator, ConstIterator> *_v;
		typename Iterator::const_word_iterator _pos;
		value_type _word;
		bool _ref_valid;

		inline void update_ref ()
		{
			if (!_ref_valid) {
				if (_v->_begin.pos () == 0)
					_word = *_pos;
				else if (_pos + 1 == (_v->_end.pos () ? (_v->_end.word () + 1) : _v->_end.word ()))
					_word = Endianness::shift_left (*_pos, _v->_begin.pos ());
				else
					_word = Endianness::shift_left (*_pos, _v->_begin.pos ()) | Endianness::shift_right (_pos[1], WordTraits<value_type>::bits - _v->_begin.pos ());

				if (((_pos - _v->_begin.word () + 1) << WordTraits<value_type>::logof_size) > _v->_end - _v->_begin)
					_word &= Endianness::mask_left ((_v->_end - _v->_begin) & WordTraits<value_type>::pos_mask);

				_ref_valid = true;
			}
		}
	};

}; // template <class Iterator> class BitSubvector

} // namespace LinBox

#endif // __VECTOR_BIT_SUBVECTOR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
