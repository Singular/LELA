/* -*- mode: c++;-*-
 * test-bit-subvector.C
 * Tests for BitSubvector
 * Written by Bradford Hovinen <hovinen@gmail.com>
 * Copyright 2011 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector.h"

namespace F4Tests {

using namespace LinBox;

typedef __LINBOX_BITVECTOR_WORD_TYPE word;

static const word pattern[2] = { 0xf0f0f0f0f0f0f0f0ULL, 0xaaaaaaaaaaaaaaaaULL };

template <class Endianness>
word connect (word word1, word word2, int shift) 
{
	if (shift == 0)
		return word1;
	else
		return Endianness::shift_left (word1, shift) | Endianness::shift_right (word2, __LINBOX_BITSOF_LONG - shift);
}

// To be tested:
//  - Termination of word-iterator (gets right number of words at all offsets, all bit-lengths modulo __LINBOX_BITSOF_LONG)
//  - Bit-vector doesn't end at word-boundary - correct masking
//  - BIt-subvector doesn't have length which is multiple of word-length

void testConstIterator ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	const int n = 256;
	const int k = 128;

	BitVector<> v (n);

	BitVector<>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<>::const_iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<>::const_iterator>::const_word_iterator i;

		flip = ((offset / __LINBOX_BITSOF_LONG) % 2 == 0) ? 0 : 1;
				
		for (i = vp.wordBegin (); i != vp.wordEnd (); ++i, flip = 1 - flip) {
			word check = connect<BitVector<>::Endianness> (pattern[flip], pattern[1-flip], offset % __LINBOX_BITSOF_LONG);
			
			if (*i != check) {
				std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
				std::cerr << __FUNCTION__ << ": Pattern should be " << std::hex << check << std::endl;
				std::cerr << __FUNCTION__ << ": Detected " << std::hex << *i << std::endl;
				return;
			}
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;
}

void testIterator ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	const int n = 256;
	const int k = 128;

	BitVector<> v (n);

	BitVector<>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<>::iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<>::iterator>::word_iterator i;
		BitSubvector<BitVector<>::iterator>::const_word_iterator j;

		flip = ((offset / __LINBOX_BITSOF_LONG) % 2 == 0) ? 0 : 1;
				
		for (i = vp.wordBegin (), j = vp.wordBegin (); i != vp.wordEnd (); ++i, ++j, flip = 1 - flip) {
			if (*j != *i) {
				std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
				std::cerr << __FUNCTION__ << ": word_iterator and const_word_iterator don't agree" << std::endl;
				std::cerr << __FUNCTION__ << ": word_iterator: " << std::hex << *i << std::endl;
				std::cerr << __FUNCTION__ << ": const_word_iterator: " << std::hex << *j << std::endl;
				return;
			}

			word val = *j;
			
			*i ^= val;
			
			if (*j != 0) {
				std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
				std::cerr << __FUNCTION__ << ": Pattern should be 0 after clearing" << std::endl;
				std::cerr << __FUNCTION__ << ": Detected " << std::hex << *j << std::endl;
				return;
			}

			*i ^= val;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;
}

} // namespace F4Tests

int main (int argc, char **argv)
{
	F4Tests::testConstIterator ();
	F4Tests::testIterator ();
}
