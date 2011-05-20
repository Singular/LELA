/* linbox/vector/subvector.h
 * Copyright 2002 William J. Turner
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 *
 * Written by William J. Turner <wjturner@acm.org>
 * Mods by -bds 
 * Maintained by: -bds <saunders@udel.edu> 
 * (where there is missing or buggy function, please contact me rather than workaround)
 */

#ifndef __LINBOX_subvector_H
#define __LINBOX_subvector_H

#include <iterator>
#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/subiterator.h"

namespace LinBox
{

//wrapper Iterator to get a const Iterator

/** \brief Dense subvector 
\ingroup vector

 * This class provides a statically sized subvector of a 
 * random access container (such as std::vector, deque).
 * It does not work on sparse linbox vectors.
 * It implements all of the types and methods of a std::vector
 * except for those that invalidate iterators, i.e.,
 * those (potentially) involving vector resizing, such as
 * push_back(), insert(), resize().
 */
template <typename Iterator, typename ConstIterator = Iterator> 
class Subvector //: public Vector // for types
{
    public:
	// Types
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	typedef typename VectorCategories::DenseVectorTag VectorCategory;
    
	typedef size_t                                              size_type;
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
	typedef typename std::iterator_traits<Iterator>::pointer	    pointer;
	typedef typename std::iterator_traits<Iterator>::reference	    reference;
	typedef const reference	                                    const_reference; 
	typedef Iterator                                            iterator;
	//typedef typename ConstIteratorType<Iterator>::const_iterator         const_iterator;

	typedef ConstIterator                                            const_iterator;

	typedef std::reverse_iterator<iterator>	                    reverse_iterator;
	typedef std::reverse_iterator<const_iterator>               const_reverse_iterator;
   	
	Subvector () : _begin (Iterator ()), _end (Iterator ()) {}
   	
	template<class Vector>
	Subvector (Vector &v, size_type start, size_type stride, size_type length)
		: _begin (iterator (v.begin() + start, stride)),
		  _end   (iterator (v.begin() + start + (stride * length), stride)) {}
	
	Subvector (iterator begin, iterator end)
		: _begin (begin), _end (end) {}
	
	Subvector (iterator begin, size_type length)
		: _begin (begin), _end (begin + length) {}
	
	//copy constructor
	Subvector (const Subvector &x)
		: _begin (x._begin), _end (x._end) {}

	template <class It, class CIt>
	Subvector (const Subvector<It, CIt> &x)
		: _begin (x._begin), _end (x._end) {}
	
	~Subvector () {}
	
	// Iterators
	
	iterator               begin  ()       { return _begin; }
	const_iterator         begin  () const { return _begin; }
	iterator               end    ()       { return _end; }
	const_iterator         end    () const { return _end; }
	
	reverse_iterator       rbegin ()       { return reverse_iterator (_end); }
	const_reverse_iterator rbegin () const { return reverse_iterator (_end); }
	reverse_iterator       rend   ()       { return reverse_iterator (_begin); }
	const_reverse_iterator rend   () const { return reverse_iterator (_begin); }
	
	// Element access
	
	reference       operator[] (size_type n)       { return _begin[n]; }
	const_reference operator[] (size_type n) const { return _begin[n]; }
	
	// the method "at" does appear to be implemented 
	// in the gnu implementation of the STL
	reference at(size_type n)  // validity is relative to valid _begin, _end
	{   
		iterator p = _begin + n;
		if ( _begin <= p && p < _end ) 
			return *p;
		else 
			throw std::out_of_range("out of range"); //out of range error message.
	}
	
	const_reference at(size_type n) const 
	{
		const_iterator p = _begin + n;
		if ( _begin <= p && p < _end)
			return *p;
		else 
			throw std::out_of_range("out of range"); //out of range error message
	}
	
	reference       front ()       { return *_begin; }
	const_reference front () const { return *_begin; }
	reference       back  ()       { return *( _end - 1 ); }
	const_reference back  () const { return *( _end - 1 ); }
	
	template<class Container>
	/** assign the elements of Container one by one to *this.
	 *  Container must be at least as long as this.
	 */
	Subvector& operator= (const Container& x)
	{
		typename Container::const_iterator q = x.begin ();
		
		for (iterator p = _begin (); p != _end (); ++p, ++q)
			*p = *q;
		
		return *this;
	}
	
	template<class Iterator2, class ConstIterator2>
	Subvector& operator = (const Subvector<Iterator2, ConstIterator2>& sub)
	{
		_begin=sub.begin();
		_end=sub.end();
		return *this;
	}
	
	//		template <class In> void assign(In first, In last);
	//		void assign(size_type n, const T& val);
	
	// Stack operations:  
	// 	not implemented because they invalidate iterators
	
	// List operations:
	// 	not implemented because they invalidate iterators
	
	// Capacity
	// 	resize, reserve: not implemented because they 
	// 		invalidate iterators
	
	//copy assignment
	Subvector& operator=(const Subvector& sub)
	{ _begin = sub._begin; _end = sub._end; return *this; }
	
	size_type size      () const { return _end - _begin; }
	bool      empty     () const { return _end == _begin; }
	size_type max_size  () const { return _end - _begin; }
	size_type capacity  () const { return _end - _begin; }
	
    protected:
	template <class It, class CIt>
	friend class Subvector;
	
	iterator _begin; // a subiterator of wrapped vector
	iterator _end;	 // a subiterator of wrapped vector
	
}; // template <class Vector> class Subvector

/** \brief Mutable subvector 
 *
 * This is a mutable version of @ref Subvector, providing operations
 * which invalidate iterators. It is not available for const vectors.
 */
template <typename Vector> 
class MutableSubvector
{
    public:
	// Types
	typedef typename VectorCategories::DenseVectorTag VectorCategory;

	typedef size_t                                                   size_type;
	typedef typename Vector::iterator                                iterator;
	typedef typename Vector::const_iterator                          const_iterator;
	typedef typename std::iterator_traits<iterator>::value_type      value_type;
	typedef typename std::iterator_traits<iterator>::difference_type difference_type;
	typedef typename std::iterator_traits<iterator>::pointer	 pointer;
	typedef typename std::iterator_traits<iterator>::reference	 reference;
	typedef const reference	                                         const_reference; 

	typedef std::reverse_iterator<iterator>                          reverse_iterator;
	typedef std::reverse_iterator<const_iterator>                    const_reverse_iterator;

	MutableSubvector () {}

	MutableSubvector (Vector &v, typename Vector::iterator begin, typename Vector::iterator end)
		: _v (&v), _idx_begin (begin - v.begin ()), _idx_end (end - v.begin ()) {}

	MutableSubvector (const MutableSubvector &v)
		: _v (v._v), _idx_begin (v._idx_begin), _idx_end (v._idx_end) {}

	~MutableSubvector () {}

	MutableSubvector &operator = (const MutableSubvector &v)
	{
		_v = v._v;
		_idx_begin = v._idx_begin;
		_idx_end = v._idx_end;
		return *this;
	}
	
	iterator               begin  ()       { return _v->begin () + _idx_begin; }
	const_iterator         begin  () const { return _v->begin () + _idx_begin; }
	iterator               end    ()       { return _v->begin () + _idx_end; }
	const_iterator         end    () const { return _v->begin () + _idx_end; }
	
	reverse_iterator       rbegin ()       { return reverse_iterator (_v->begin () + _idx_end); }
	const_reverse_iterator rbegin () const { return reverse_iterator (_v->begin () + _idx_end); }
	reverse_iterator       rend   ()       { return reverse_iterator (_v->begin () + _idx_begin); }
	const_reverse_iterator rend   () const { return reverse_iterator (_v->begin () + _idx_begin); }
	
	// Element access
	
	reference       operator[] (size_type n)       { return _v[_idx_begin + n]; }
	const_reference operator[] (size_type n) const { return _v[_idx_begin + n]; }
	
	reference at (size_type n)
	{   
		if (n < _idx_end - _idx_begin) 
			return _v[_idx_begin + n];
		else 
			throw std::out_of_range ("n");
	}
	
	const_reference at (size_type n) const 
	{
		if (n < _idx_end - _idx_begin) 
			return _v[_idx_begin + n];
		else 
			throw std::out_of_range ("n");
	}
	
	reference       front ()       { return _v[_idx_begin]; }
	const_reference front () const { return _v[_idx_begin]; }
	reference       back  ()       { return _v[_idx_end - 1]; }
	const_reference back  () const { return _v[_idx_end - 1]; }
	
	template<class Container>
	MutableSubvector &operator = (const Container &x)
	{
		typename Container::const_iterator q = x.begin ();
		iterator p_end = _v->begin () + _idx_end;
		
		for (iterator p = _v->begin () + _idx_begin; p != p_end; ++p, ++q)
			*p = *q;
		
		return *this;
	}

	size_type size     () const { return _idx_end - _idx_begin; }
	bool      empty    () const { return _idx_end == _idx_begin; }
	size_type max_size () const { return _v->max_size () - (_v->size () - _idx_end) - _idx_begin; }
	size_type capacity () const { return _v->capacity () - (_v->size () - _idx_end) - _idx_begin; }

	iterator insert (iterator position, const value_type &x)
	{
		position = _v->insert (position, x);
		++_idx_end;
		return position;
	}

	void insert (iterator position, size_type n, const value_type &x)
	{
		_v->insert (position, n, x);
		_idx_end += n;
	}	

	template <class InputIterator>
	void insert (iterator position, InputIterator first, InputIterator last)
	{
		size_type old_size = _v->size ();
		_v->insert (position, first, last);
		_idx_end += _v->size () - old_size;
	}

	void erase (iterator position)
	{
		position = _v->erase (position);
		--_idx_end;
		return position;
	}

	void erase (iterator first, iterator last)
	{
		_idx_end -= last - first;
		first = _v->erase (first, last);
		return first;
	}

	void push_back (const value_type &x)
		{ insert (end (), x); }

	void pop_back ()
		{ erase (end () - 1); }

	void resize (size_type n)
	{
		if (n > size ()) {
			size_type old_size = _v->size ();
			_v->resize (n + _idx_begin + _v->size () - _idx_end);
			std::copy (_v->begin () + _idx_end, _v->begin () + old_size, _v->begin () + (_idx_begin + n));
		} else if (n < size ()) {
			std::copy (_v->begin () + _idx_end, _v->end (), _v->begin () + (_idx_begin + n));
			_v->resize (n + _idx_begin + _v->size () - _idx_end);
		}

		_idx_end = _idx_begin + n;
	}

	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
		{ assign_spec (first, last, typename std::iterator_traits<InputIterator>::iterator_category ()); }

	void clear () { resize (0); }

	void reserve (size_type n)
		{ _v->reserve (n + _idx_begin + _v->size () - _idx_end); }

	// Direct access to internals -- only for use by SparseSubvector
	Vector &parent     () const { return *_v; }
	void set_parent    (Vector &v)                 { _v         = &v; }
	void set_idx_begin (difference_type idx_begin) { _idx_begin = idx_begin; }
	void set_idx_end   (difference_type idx_end)   { _idx_end   = idx_end; }
	
    protected:
	template <class V2>
	friend class MutableSubvector;

	template <class V, class T>
	class SparseSubvector;

	template <class InputIterator>
	void assign_spec (InputIterator first, InputIterator last, std::input_iterator_tag)
	{
		clear ();

		while (first != last)
			insert (_v->begin () + _idx_end, *first++);
	}

	template <class InputIterator>
	void assign_spec (InputIterator first, InputIterator last, std::forward_iterator_tag)
	{
		resize (last - first);
		std::copy (first, last, _v->begin () + _idx_begin);
	}
	
	Vector *_v;
	difference_type _idx_begin, _idx_end;

}; // template <class Vector> class MutableSubvector
  
} // namespace LinBox

namespace std {

template<class Iterator, class ConstIterator>
void swap (LinBox::Subvector<Iterator, ConstIterator> &x, LinBox::Subvector<Iterator, ConstIterator> &y)
{
	linbox_check (x.size () == y.size ());

	typename LinBox::Subvector<Iterator, ConstIterator>::iterator i_x, i_y;

	for (i_x = x.begin (), i_y = y.begin (); i_x != x.end (); ++i_x, ++i_y)
		std::swap (*i_x, *i_y);
}

template<class Vector>
void swap (LinBox::MutableSubvector<Vector> &x, LinBox::MutableSubvector<Vector> &y)
{
	linbox_check (x.size () == y.size ());

	typename LinBox::MutableSubvector<Vector>::iterator i_x, i_y;

	for (i_x = x.begin (), i_y = y.begin (); i_x != x.end (); ++i_x, ++i_y)
		std::swap (*i_x, *i_y);
}

}
	
#endif //__LINBOX_subvector_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

