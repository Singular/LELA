/* -*- mode: c++;-*-
 * test-bit-subvector.C
 * Test for subvector of a bit-vector
 * Written by Bradford Hovinen <hovinen@gmail.com>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include <cmath>

#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector.h"

using namespace LinBox;

typedef BitVector::word_iterator::value_type Word;

void sieve (BitVector &vec, int offset) {
	int p, m;

	BitVector::iterator i;

	for (i = vec.begin (); i != vec.end (); ++i)
		*i = true;

	if (offset < 1)
		vec[-offset] = false;
	if (offset < 2)
		vec[1-offset] = false;

	for (p = 2; p < sqrt (vec.size () + offset) + 1; ++p) {
		if (p < offset || vec[p - offset]) {
			for (m = 2 * p; m < vec.size () + offset; m += p) {
				if (m >= offset)
					vec[m - offset] = false;
			}
		}
	}
}

bool testSubvectorWordIterator ()
{
	BitVector u, v;
	BitSubvector<BitVector::iterator> w;
	Word t;
	bool fail = false;

	int i;

	std::cout << "Testing subvector word-iterator..." << std::endl;

	u.resize (2 * 8 * sizeof (Word));
	v.resize (1 * 8 * sizeof (Word));

	for (i = 0; i < 1 * 8 * sizeof (Word); ++i) {
		std::cout << "Testing subvector offset " << i << "...";
		std::cout.flush ();

		sieve (u, 0);
		sieve (v, i);
		w = BitSubvector<BitVector::iterator> (u.begin () + i, u.begin () + i + v.size ());

		if (*(u.wordBegin ()) == 0ULL) {
			std::cout << std::endl << "Error: sieved vector u is zero" << std::endl;
			fail = true;
		}

		if (*(v.wordBegin ()) == 0ULL) {
			std::cout << std::endl << "Error: sieved vector v is zero" << std::endl;
			fail = true;
		}

		if (*(w.wordBegin ()) != *(v.wordBegin ())) {
			std::cout << std::endl << "Error: subvector word is " << std::hex << *(w.wordBegin ()) << " but should be " << std::hex << *(v.wordBegin ()) << std::endl;
			fail = true;
		}

		*(w.wordBegin ()) ^= *(v.wordBegin ());

		if (*(w.wordBegin ()) != 0) {
			std::cout << std::endl << "Error: after xor, subvector word is " << std::hex << *(w.wordBegin ()) << " but should be 0" << std::endl;
			fail = true;
		}

		if (!fail)
			std::cout << "passed" << std::endl;
	}

	std::cout << "Finished testing subvector word-iterator" << std::endl;
}

int main (int argc, char **argv)
{
	testSubvectorWordIterator ();
}
