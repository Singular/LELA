/* tests/test-bit-subvector.C
 * Copyright 2011 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Tests for BitSubvector
 *
 * ---------------------------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include <iostream>
#include <sstream>

#include "test-common.h"

#include "lela/util/commentator.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/bit-subvector.h"

using namespace LELA;

typedef BitVector<>::Endianness Endianness;

static const uint64 pattern[2] = { 0xf0f0f0f0f0f0f0f0ULL, 0xaaaaaaaaaaaaaaaaULL };

uint64 connect (uint64 word1, uint64 word2, int shift) 
{
	if (shift == 0)
		return word1;
	else
		return Endianness::shift_left (word1, shift) | Endianness::shift_right (word2, WordTraits<uint64>::bits - shift);
}

// To be tested:
//  - Bits before subvector not touched when subvector modified

bool testConstIterator (size_t n, size_t k)
{
	commentator.start ("Testing BitSubvector::const_word_iterator", __FUNCTION__);

	bool pass = true;

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.word_begin (); w != v.word_end (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	v.back_word () = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<Endianness>::const_iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<Endianness>::const_iterator>::const_word_iterator i;

		flip = ((offset / WordTraits<BitVector<Endianness>::word_type>::bits) % 2 == 0) ? 0 : 1;

		size_t idx = 0;
				
		for (i = vp.word_begin (); i != vp.word_end (); ++i, flip = 1 - flip, idx += WordTraits<BitVector<Endianness>::word_type>::bits) {
			BitVector<Endianness>::word_type check = connect (pattern[flip], pattern[1-flip], offset % WordTraits<BitVector<Endianness>::word_type>::bits);

			if (idx + WordTraits<BitVector<Endianness>::word_type>::bits > k)
				check &= Endianness::mask_left (k % WordTraits<BitVector<Endianness>::word_type>::bits);
			
			if (*i != check) {
				error << "ERROR: error at offset " << offset << std::endl;
				error << "ERROR: Pattern should be " << std::hex << check << std::endl;
				error << "ERROR: Detected " << std::hex << *i << std::endl;
				pass = false;
				goto end_test;
			}
		}
	}

end_test:
	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int testIterator (size_t n, size_t k)
{
	commentator.start ("Testing BitSubvector::word_iterator", __FUNCTION__);

	bool pass = true;

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.word_begin (); w != v.word_end (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	v.back_word () = pattern[flip];

	size_t offset;
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<Endianness>::iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<Endianness>::iterator>::word_iterator i;
		BitSubvector<BitVector<Endianness>::iterator>::const_word_iterator j1, j2;

		flip = ((offset / WordTraits<BitVector<Endianness>::word_type>::bits) % 2 == 0) ? 0 : 1;
				
		for (i = vp.word_begin (), j1 = vp.word_begin (), j2 = vp.word_begin (); i != vp.word_end (); ++i, ++j1, ++j2, flip = 1 - flip) {
			if (*j1 != *i) {
				error << "ERROR: error at offset " << offset << " (n = " << n << ", k = " << k << ")" << std::endl;
				error << "ERROR: word_iterator and const_word_iterator don't agree" << std::endl;
				error << "ERROR: word_iterator: " << std::hex << *i << std::endl;
				error << "ERROR: const_word_iterator: " << std::hex << *j1 << std::endl;
				pass = false;
				goto end_test;
			}

			BitVector<Endianness>::word_type val = *j1;
			
			*i ^= val;
			
			if (*j2 != 0) {
				error << "ERROR: error at offset " << offset << " (n = " << n << ", k = " << k << ")" << std::endl;
				error << "ERROR: Pattern should be 0 after clearing" << std::endl;
				error << "ERROR: Detected " << std::hex << *j2 << std::endl;
				pass = false;
				goto end_test;
			}

			*i ^= val;
		}
	}

end_test:
	commentator.stop (MSG_STATUS (pass));

	return pass;
}

int testBackWord ()
{
	commentator.start ("Testing back_word", __FUNCTION__);

	bool pass = true;

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	static const int n = 3 * WordTraits<BitVector<Endianness>::word_type>::bits;
	static const int k = 16;

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	
	size_t offset;
	unsigned int flip = 0;
	
	for (offset = 0; offset < WordTraits<BitVector<Endianness>::word_type>::bits; ++offset, flip = 1 - flip) {
		BitSubvector<BitVector<Endianness>::iterator> vp1 (v.begin () + offset, v.begin () + (offset + k));
		BitSubvector<BitVector<Endianness>::const_iterator> vp2 (v.begin () + offset, v.begin () + (offset + WordTraits<BitVector<Endianness>::word_type>::bits));

		vp1.back_word () = pattern[flip];

		BitVector<Endianness>::word_type check = pattern[flip] & Endianness::mask_left (k);

		if (vp1.back_word () != check) {
			error << "ERROR: error at offset " << std::dec << offset << " (rereading back_word)" << std::endl;
			error << "ERROR: Pattern should be " << std::hex << check << std::endl;
			error << "ERROR: Detected " << std::hex << vp1.back_word () << std::endl;
			pass = false;
		}

		if (vp2.front_word () != check) {
			error << "ERROR: error at offset " << std::dec << offset << " (checking against a const subvector)" << std::endl;
			error << "ERROR: Pattern should be " << std::hex << check << std::endl;
			error << "ERROR: Detected " << std::hex << vp2.front_word () << std::endl;
			pass = false;
		}
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

bool testWordLength (size_t n, size_t k)
{
	commentator.start ("Testing word-length", __FUNCTION__);

	bool pass = true;

	std::ostream &error = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	BitVector<Endianness> v (n);

	BitVector<Endianness>::word_iterator w;
	unsigned int flip = 0;
	
	for (w = v.word_begin (); w != v.word_end (); ++w, flip = 1 - flip)
		*w = pattern[flip];

	size_t offset;

	size_t correct_len = k / WordTraits<BitVector<Endianness>::word_type>::bits;

	if (k % WordTraits<BitVector<Endianness>::word_type>::bits != 0)
		++correct_len;

	--correct_len; // New semantics
	
	for (offset = 0; offset <= n - k; ++offset) {
		BitSubvector<BitVector<Endianness>::iterator> vp (v.begin () + offset, v.begin () + (offset + k));

		BitSubvector<BitVector<Endianness>::iterator>::word_iterator i;
		BitSubvector<BitVector<Endianness>::iterator>::const_word_iterator j;

		size_t len;

		for (i = vp.word_begin (), len = 0; i != vp.word_end (); ++i, ++len);
		
		if (len != correct_len) {
			error << "ERROR: error at offset " << offset << std::endl;
			error << "ERROR: length for word_iterator should be " << correct_len << std::endl;
			error << "ERROR: but computed " << len << std::endl;
			pass = false;
		}

		for (j = vp.word_begin (), len = 0; j != vp.word_end (); ++j, ++len);
		
		if (len != correct_len) {
			error << "ERROR: error at offset " << offset << std::endl;
			error << "ERROR: length for const_word_iterator should be " << correct_len << std::endl;
			error << "ERROR: but computed " << len << std::endl;
			pass = false;
		}
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

bool runTests () 
{
	bool pass = true;

	size_t n, k;

	std::ostringstream str;

	testBackWord ();
	
	for (n = 256; n < 256 + WordTraits<BitVector<Endianness>::word_type>::bits; ++n) {
		for (k = 128; k < 128 + WordTraits<BitVector<Endianness>::word_type>::bits; ++k) {
			str << "Running tests for main vector-length " << n << ", subvector-length " << k << std::ends;
			commentator.start (str.str ().c_str (), __FUNCTION__);
			pass = testWordLength (n, k);
			pass = testConstIterator (n, k) && pass;
			pass = testIterator (n, k) && pass;
			commentator.stop (MSG_STATUS (pass));

			if (!pass)
				goto end_test;
		}
	}

end_test:

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	commentator.start ("BitSubvector test-suite", "BitSubvector");

	pass = runTests () && pass;

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
