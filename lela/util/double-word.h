/* lela/util/double-word.h
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ----------------------------
 *
 * Note: This header must not have dependencies, since
 * lela/integer.h depends on it to provide a 128-bit data-type when
 * one is not provided by the system.
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_DOUBLE_WORD_H
#define __LELA_DOUBLE_WORD_H

namespace LELA
{

/** Class which given an integral data-type, constructs a data-type which is twice as wide
 *
 * This class provides only bitwise binary operators &, |, and ^,
 * shift operators >> and <<, and their assignment-equivalents. As in
 * the C++-standard, shifting by greater than the size of this
 * datatype (which is twice the size of its parameter) is undefined.
 */
template <class word>
struct DoubleWord {
	word _low, _high;

	// Workaround for the fact that non-trivial constructors are
	// not allowed in unions -- this class appears in a union in
	// lela/vector/bit-iterator.h
	static DoubleWord make_double_word (word low, word high)
	{
		DoubleWord w;
		w._low = low;
		w._high = high;
		return w;
	}

	DoubleWord operator & (const DoubleWord &w) const
		{ return make_double_word (_low & w._low, _high & w._high); }

	DoubleWord operator | (const DoubleWord &w) const
		{ return make_double_word (_low | w._low, _high | w._high); }

	DoubleWord operator ^ (const DoubleWord &w) const
		{ return make_double_word (_low ^ w._low, _high ^ w._high); }

	DoubleWord operator << (unsigned int shift) const
	{
		if (shift == 0)
			return make_double_word (_low, _high);
		else if (shift >= sizeof (word) * 8)
			return make_double_word (0, _low << (shift - sizeof (word) * 8));
		else
			return make_double_word (_low << shift, (_low >> (sizeof (word) * 8 - shift)) | (_high << shift));
	}

	DoubleWord operator >> (unsigned int shift) const
	{
		if (shift == 0)
			return make_double_word (_low, _high);
		else if (shift >= sizeof (word) * 8)
			return make_double_word (_high >> (shift - sizeof (word) * 8), 0);
		else
			return make_double_word ((_low >> shift) | (_high << (sizeof (word) * 8 - shift)), _high >> shift);
	}

	DoubleWord &operator &= (const DoubleWord &w)
		{ _low &= w._low; _high &= w._high; return *this; }

	DoubleWord &operator |= (const DoubleWord &w)
		{ _low |= w._low; _high |= w._high; return *this; }

	DoubleWord &operator ^= (const DoubleWord &w)
		{ _low ^= w._low; _high ^= w._high; return *this; }

	DoubleWord &operator <<= (unsigned int shift)
	{
		if (shift >= sizeof (word) * 8) {
			_low = 0;
			_high = _low << (shift - sizeof (word) * 8);
		}
		else if (shift > 0) {
			_low <<= shift;
			_high = (_high << shift) | (_low >> (sizeof (word) * 8 - shift));
		}

		return *this;
	}

	DoubleWord &operator >>= (unsigned int shift)
	{
		if (shift >= sizeof (word) * 8) {
			_low = _high >> (shift - sizeof (word) * 8);
			_high = 0;
		}
		else if (shift > 0) {
			_low = (_low >> shift) | (_high << (sizeof (word) * 8 - shift));
			_high >>= shift;
		}

		return *this;
	}
};

} // namespace LELA

#endif // __LELA_DOUBLE_WORD_H
