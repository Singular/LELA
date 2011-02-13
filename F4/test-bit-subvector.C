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

int testConstIterator (int n, int k)
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

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
				return -1;
			}
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;

	return 0;
}

int testIterator (int n, int k)
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

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
				return -1;
			}

			word val = *j;
			
			*i ^= val;
			
			if (*j != 0) {
				std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
				std::cerr << __FUNCTION__ << ": Pattern should be 0 after clearing" << std::endl;
				std::cerr << __FUNCTION__ << ": Detected " << std::hex << *j << std::endl;
				return -1;
			}

			*i ^= val;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;

	return 0;
}

int testWordLength (int n, int k)
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	BitVector<> v (n);

	BitVector<>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;

	size_t correct_len = k / __LINBOX_BITSOF_LONG;

	if (k % __LINBOX_BITSOF_LONG != 0)
		++correct_len;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<>::iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<>::iterator>::word_iterator i;
		BitSubvector<BitVector<>::iterator>::const_word_iterator j;

		size_t len;

		for (i = vp.wordBegin (), len = 0; i != vp.wordEnd (); ++i, ++len);
		
		if (len != correct_len) {
			std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
			std::cerr << __FUNCTION__ << ": length for word_iterator should be " << correct_len << std::endl;
			std::cerr << __FUNCTION__ << ": but computed " << len << std::endl;
		}

		for (j = vp.wordBegin (), len = 0; j != vp.wordEnd (); ++j, ++len);
		
		if (len != correct_len) {
			std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
			std::cerr << __FUNCTION__ << ": length for const_word_iterator should be " << correct_len << std::endl;
			std::cerr << __FUNCTION__ << ": but computed " << len << std::endl;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;

	return 0;
}

void runTests () 
{
	int n, k;
	
	for (n = 256; n < 256 + __LINBOX_BITSOF_LONG; ++n) {
		for (k = 128; k < 128 + __LINBOX_BITSOF_LONG; ++k) {
			std::cout << __FUNCTION__ << ": Running tests for main vector-length " << n << ", subvector-length " << k << std::endl;
			if (testWordLength (n, k) == -1 || testConstIterator (n, k) == -1 || testIterator (n, k) == -1)
				return;
		}
	}
}

} // namespace F4Tests

int main (int argc, char **argv)
{
	F4Tests::runTests ();
}
