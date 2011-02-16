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

typedef BigEndian<word> Endianness;

static const word pattern[2] = { 0xf0f0f0f0f0f0f0f0ULL, 0xaaaaaaaaaaaaaaaaULL };

word connect (word word1, word word2, int shift) 
{
	if (shift == 0)
		return word1;
	else
		return Endianness::shift_left (word1, shift) | Endianness::shift_right (word2, WordTraits<word>::bits - shift);
}

// To be tested:
//  - Bits before subvector not touched when subvector modified

int testConstIterator (int n, int k)
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<Endianness>::const_iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<Endianness>::const_iterator>::const_word_iterator i;

		flip = ((offset / WordTraits<word>::bits) % 2 == 0) ? 0 : 1;

		size_t idx = 0;
				
		for (i = vp.wordBegin (); i != vp.wordEnd (); ++i, flip = 1 - flip, idx += WordTraits<word>::bits) {
			word check = connect (pattern[flip], pattern[1-flip], offset % WordTraits<word>::bits);

			if (idx + WordTraits<word>::bits > k)
				check &= Endianness::mask_left (k % WordTraits<word>::bits);
			
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

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<Endianness>::iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<Endianness>::iterator>::word_iterator i;
		BitSubvector<BitVector<Endianness>::iterator>::const_word_iterator j, k;

		flip = ((offset / WordTraits<word>::bits) % 2 == 0) ? 0 : 1;
				
		for (i = vp.wordBegin (), j = vp.wordBegin (), k = vp.wordBegin (); i != vp.wordEnd (); ++i, ++j, ++k, flip = 1 - flip) {
			if (*j != *i) {
				std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
				std::cerr << __FUNCTION__ << ": word_iterator and const_word_iterator don't agree" << std::endl;
				std::cerr << __FUNCTION__ << ": word_iterator: " << std::hex << *i << std::endl;
				std::cerr << __FUNCTION__ << ": const_word_iterator: " << std::hex << *j << std::endl;
				return -1;
			}

			word val = *j;
			
			*i ^= val;
			
			if (*k != 0) {
				std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
				std::cerr << __FUNCTION__ << ": Pattern should be 0 after clearing" << std::endl;
				std::cerr << __FUNCTION__ << ": Detected " << std::hex << *k << std::endl;
				return -1;
			}

			*i ^= val;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;

	return 0;
}

int testMask ()
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	static const int n = 3 * WordTraits<word>::bits;
	static const int k = 16;

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset < WordTraits<word>::bits; ++offset) {
		BitSubvector<BitVector<Endianness>::iterator> vp1 (v.begin () + offset, v.begin () + (offset + k));
		BitSubvector<BitVector<Endianness>::const_iterator> vp2 (v.begin () + offset + k, v.begin () + (offset + k + WordTraits<word>::bits));

		BitSubvector<BitVector<Endianness>::iterator>::word_iterator i = vp1.wordBegin ();
		BitSubvector<BitVector<Endianness>::const_iterator>::const_word_iterator j = vp2.wordBegin ();

		*i ^= pattern[0];

		if (k + offset >= WordTraits<word>::bits)
			flip = 1;
		else
			flip = 0;

		word check = connect (pattern[flip], pattern[1 - flip], (k + offset) % WordTraits<word>::bits);

		if (*j != check) {
			std::cerr << __FUNCTION__ << ": error at offset " << offset << std::endl;
			std::cerr << __FUNCTION__ << ": Pattern should be " << std::hex << check << std::endl;
			std::cerr << __FUNCTION__ << ": Detected " << std::hex << *i << std::endl;
			return -1;
		}
	}

	std::cout << __FUNCTION__ << ": done" << std::endl;

	return 0;
}

int testWordLength (int n, int k)
{
	std::cout << __FUNCTION__ << ": Enter" << std::endl;

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.wordBegin (); w != v.wordEnd (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;

	size_t correct_len = k / WordTraits<word>::bits;

	if (k % WordTraits<word>::bits != 0)
		++correct_len;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<Endianness>::iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<Endianness>::iterator>::word_iterator i;
		BitSubvector<BitVector<Endianness>::iterator>::const_word_iterator j;

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

	testMask ();
	
	for (n = 256; n < 256 + WordTraits<word>::bits; ++n) {
		for (k = 128; k < 128 + WordTraits<word>::bits; ++k) {
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
