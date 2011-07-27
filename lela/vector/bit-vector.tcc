/* lela/vector/bit-vector.tcc
 * Copyright 2003 Bradford Hovinen
 *
 * Dense vectors over GF(2)
 *
 * -------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_BIT_VECTOR_TCC
#define __LELA_BIT_VECTOR_TCC

#include <stdexcept>
#include <vector>

#include "lela/vector/bit-vector.h"

namespace LELA
{

template <class Endianness>
inline typename BitVector<Endianness>::iterator BitVector<Endianness>::end (void)
{
	if ((_size & WordTraits<word_type>::pos_mask) == 0UL)
		return iterator (_v.end (), 0UL);
	else
		return iterator (_v.end () - 1UL, _size & WordTraits<word_type>::pos_mask);
}

template <class Endianness>
inline typename BitVector<Endianness>::const_iterator BitVector<Endianness>::end (void) const
{
	if ((_size & WordTraits<word_type>::pos_mask) == 0UL)
		return const_iterator (_v.end (), 0UL);
	else
		return const_iterator (_v.end () - 1UL, _size & WordTraits<word_type>::pos_mask);
}

template <class Endianness>
inline typename BitVector<Endianness>::reference BitVector<Endianness>::at (BitVector<Endianness>::size_type n)
{
	if (n >= _size)
		throw std::out_of_range ("LELA::BitVector");
	else
		return (*this)[n];
}

template <class Endianness>
inline typename BitVector<Endianness>::const_reference BitVector<Endianness>::at (BitVector<Endianness>::size_type n) const
{
	if (n >= _size)
		throw std::out_of_range ("LELA::BitVector");
	else
		return (*this)[n];
}

template <class Endianness>
inline typename BitVector<Endianness>::reference BitVector<Endianness>::back (void)
{
	if ( (_size & WordTraits<word_type>::pos_mask) == 0UL)
		return reference (_v.end (), 0UL);
	else
		return reference (_v.end () - 1UL, _size & WordTraits<word_type>::pos_mask);
}

template <class Endianness>
inline typename BitVector<Endianness>::const_reference BitVector<Endianness>::back (void) const
{
	if ( (_size & WordTraits<word_type>::pos_mask) == 0UL)
		return const_reference (_v.end (), 0UL);
	else
		return const_reference (_v.end () - 1UL, _size & WordTraits<word_type>::pos_mask);
}

template <class Endianness>
inline void BitVector<Endianness>::push_back (bool v)
{
	if (v) {
		if ((_size & WordTraits<word_type>::pos_mask) == 0UL)
			push_word_back (1UL);
		else
			_v.back () |= Endianness::e_j (_size & WordTraits<word_type>::pos_mask);
	} else {
		if ((_size & WordTraits<word_type>::pos_mask) == 0UL)
			push_word_back (0UL);
		else
			_v.back () &= Endianness::mask_left (_size & WordTraits<word_type>::pos_mask);
	}

	++_size;
}

template <class Endianness>
template<class Container>
inline BitVector<Endianness> &BitVector<Endianness>::operator = (const Container &v)
{
	typename Container::const_iterator i;
	typename Container::const_iterator i_end = v.begin () + (v.size () >> WordTraits<word_type>::logof_size);
	typename std::vector<word_type>::iterator j;
	word_type idx;

	_v.resize ((v.size () >> WordTraits<word_type>::logof_size) + ((v.size () & WordTraits<word_type>::pos_mask) ? 1UL : 0UL));

	for (j = _v.begin (); i != i_end; ++j) {
		*j = 0UL;
		for (idx = 0UL; idx < WordTraits<word_type>::bits; ++idx, ++i) {
			*j <<= 1UL;
			*j |= *i & 1UL;
		}
	}

	if (v.size () & WordTraits<word_type>::pos_mask) {
		*j = 0UL;

		for (idx = 0UL; idx < (v.size () & WordTraits<word_type>::pos_mask); ++idx) {
			*j <<= 1UL;
			*j |= *i & 1UL;
		}
	}

	_size = v.size ();

	return *this;
}

template <class Endianness>
inline bool BitVector<Endianness>::operator == (const BitVector<Endianness> &v) const
{
	const_word_iterator i, j;
	word_type mask;

	if (_size != v._size) return false;

	for (i = word_begin (), j = v.word_begin (); i != word_end (); ++i, ++j)
		if (*i != *j) return false;

	mask = (1UL << (_size & (8 * sizeof (word_type) - 1UL))) - 1UL;
	if (mask == 0UL) mask = (word_type) -1UL;

	if ((back_word () & mask) == (v.back_word () & mask))
		return true;
	else
		return false;
}

} // namespace LELA

#endif // __LELA_BIT_VECTOR_TCC

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
