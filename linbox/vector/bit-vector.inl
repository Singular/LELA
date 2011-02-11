/* linbox/vector/bit-vector.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_bit_vector_INL
#define __LINBOX_bit_vector_INL

#include <stdexcept>
#include <vector>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"

namespace LinBox
{

template <class Endianness>
inline typename BitVector<Endianness>::iterator BitVector<Endianness>::begin (void)
	{ return iterator (_v.begin (), 0UL); }

template <class Endianness>
inline typename BitVector<Endianness>::const_iterator BitVector<Endianness>::begin (void) const
	{ return const_iterator (_v.begin (), 0UL); }

template <class Endianness>
inline typename BitVector<Endianness>::iterator BitVector<Endianness>::end (void)
{
	if ((_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return iterator (_v.end (), 0UL);
	else
		return iterator (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

template <class Endianness>
inline typename BitVector<Endianness>::const_iterator BitVector<Endianness>::end (void) const
{
	if ((_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return const_iterator (_v.end (), 0UL);
	else
		return const_iterator (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

template <class Endianness>
inline typename BitVector<Endianness>::reverse_iterator BitVector<Endianness>::rbegin (void)
	{ return reverse_iterator (end () - 1UL); }

template <class Endianness>
inline typename BitVector<Endianness>::const_reverse_iterator BitVector<Endianness>::rbegin (void) const
	{ return const_reverse_iterator (end () - 1UL); }

template <class Endianness>
inline typename BitVector<Endianness>::reverse_iterator BitVector<Endianness>::rend (void)
	{ return reverse_iterator (begin () - 1UL); }

template <class Endianness>
inline typename BitVector<Endianness>::const_reverse_iterator BitVector<Endianness>::rend (void) const
	{ return const_reverse_iterator (begin () - 1UL); }

template <class Endianness>
inline typename BitVector<Endianness>::reference BitVector<Endianness>::operator[] (BitVector::size_type n)
	{ return *(begin () + n); }

template <class Endianness>
inline typename BitVector<Endianness>::const_reference BitVector<Endianness>::operator[] (BitVector::size_type n) const
	{ return *(begin () + n); }

template <class Endianness>
inline typename BitVector<Endianness>::reference BitVector<Endianness>::at (BitVector<Endianness>::size_type n)
{
	if (n >= _size)
		throw std::out_of_range ("LinBox::BitVector");
	else
		return (*this)[n];
}

template <class Endianness>
inline typename BitVector<Endianness>::const_reference BitVector<Endianness>::at (BitVector<Endianness>::size_type n) const
{
	if (n >= _size)
		throw std::out_of_range ("LinBox::BitVector");
	else
		return (*this)[n];
}

template <class Endianness>
inline typename BitVector<Endianness>::reference BitVector<Endianness>::front (void)
	{ return reference (_v.begin (), 0UL); }

template <class Endianness>
inline typename BitVector<Endianness>::const_reference BitVector<Endianness>::front (void) const
	{ return const_reference (_v.begin (), 0UL); }

template <class Endianness>
inline typename BitVector<Endianness>::reference BitVector<Endianness>::back (void)
{
	if ( (_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return reference (_v.end (), 0UL);
	else
		return reference (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

template <class Endianness>
inline typename BitVector<Endianness>::const_reference BitVector<Endianness>::back (void) const
{
	if ( (_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return const_reference (_v.end (), 0UL);
	else
		return const_reference (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

template <class Endianness>
inline void BitVector<Endianness>::push_back (bool v)
{
	if (v) {
		if ((_size & __LINBOX_POS_ALL_ONES) == 0UL)
			push_word_back (1UL);
		else
			_v.back () |= 1 << (_size & __LINBOX_POS_ALL_ONES);
	} else {
		if ((_size & __LINBOX_POS_ALL_ONES) == 0UL)
			push_word_back (0UL);
		else
			_v.back () &= (1 << (_size & __LINBOX_POS_ALL_ONES)) - 1;
	}

	++_size;
}

template <class Endianness>
template<class Container>
inline BitVector<Endianness> &BitVector<Endianness>::operator = (const Container &v)
{
	typename Container::const_iterator i;
	typename Container::const_iterator i_end = v.begin () + (v.size () >> __LINBOX_LOGOF_SIZE);
	std::vector<__LINBOX_BITVECTOR_WORD_TYPE>::iterator j;
	__LINBOX_BITVECTOR_WORD_TYPE idx;

	_v.resize ((v.size () >> __LINBOX_LOGOF_SIZE) + ((v.size () & __LINBOX_POS_ALL_ONES) ? 1UL : 0UL));

	for (j = _v.begin (); i != i_end; ++j) {
		*j = 0UL;
		for (idx = 0UL; idx < __LINBOX_BITSOF_LONG; ++idx, ++i) {
			*j <<= 1UL;
			*j |= *i & 1UL;
		}
	}

	if (v.size () & __LINBOX_POS_ALL_ONES) {
		*j = 0UL;

		for (idx = 0UL; idx < (v.size () & __LINBOX_POS_ALL_ONES); ++idx) {
			*j <<= 1UL;
			*j |= *i & 1UL;
		}
	}

	_size = v.size ();

	return *this;
}

template <class Endianness>
inline void BitVector<Endianness>::resize (BitVector<Endianness>::size_type new_size, bool val)
	{ _v.resize ((new_size >> __LINBOX_LOGOF_SIZE) + ((new_size & __LINBOX_POS_ALL_ONES) ? 1UL : 0UL), val ? __LINBOX_ALL_ONES : 0UL); _size = new_size; }

template <class Endianness>
inline bool BitVector<Endianness>::operator == (const BitVector<Endianness> &v) const
{
	const_word_iterator i, j;
	__LINBOX_BITVECTOR_WORD_TYPE mask;

	if (_size != v._size) return false;

	for (i = wordBegin (), j = v.wordBegin (); i != wordEnd () - 1UL; ++i, ++j)
		if (*i != *j) return false;

	mask = (1UL << (_size & (8 * sizeof (__LINBOX_BITVECTOR_WORD_TYPE) - 1UL))) - 1UL;
	if (mask == 0UL) mask = (__LINBOX_BITVECTOR_WORD_TYPE) -1UL;

	if ((*i & mask) == (*j & mask))
		return true;
	else
		return false;
}

} // namespace LinBox

#endif // __LINBOX_bit_vector_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
