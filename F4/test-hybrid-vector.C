/* -*- mode: c++;-*-
 * test-hybrid-vector.C
 * Test for hybrid sparse-dense vector-format
 * Written by Bradford Hovinen <hovinen@math.uni-hannover.de>
 * Copyright 2010 Bradford Hovinen
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2.0 or greater
 */

#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"

using namespace LinBox;

typedef GF2 Field;
typedef BigEndian<__LINBOX_BITVECTOR_WORD_TYPE> Endianness;
typedef std::pair<std::vector<uint16>, BitVector<Endianness> > HybridVector;
typedef HybridVector::second_type::word_iterator::value_type Word;

bool testAdd ()
{
	Field F (2);
	VectorDomain<Field> VD (F);

	HybridVector u, v, w, t;

	std::cout << "Testing VD.add..." << std::endl;

	// Test 1: Indices in same block
	u.first.clear ();
	u.second.clear ();
	v.first.clear ();
	v.second.clear ();
	t.first.clear ();
	t.second.clear ();

	u.first.push_back (0);
	v.first.push_back (0);

	u.second.push_word_back (1);
	v.second.push_word_back (1 | (1 << 1));

	t.first.push_back (0);
	t.second.push_word_back (1 << 1);

	VD.add (w, u, v);

	std::cout << "Test 1:" << std::endl;
	std::cout << "    u = ";
	VD.write (std::cout, u) << std::endl;
	std::cout << "    v = ";
	VD.write (std::cout, v) << std::endl;
	std::cout << "u + v = ";
	VD.write (std::cout, w) << std::endl;
	std::cout << "    t = ";
	VD.write (std::cout, t) << std::endl;

	if (VD.areEqual (w, t))
		std::cout << "Test 1 okay" << std::endl;
	else
		std::cout << "Test 1 NOT okay" << std::endl;

	// Test 2: Indices in different blocks, u before v
	u.first.clear ();
	u.second.clear ();
	v.first.clear ();
	v.second.clear ();
	t.first.clear ();
	t.second.clear ();

	u.first.push_back (0);
	v.first.push_back (8 * sizeof (Word));

	u.second.push_word_back (1);
	v.second.push_word_back (1);

	t.first.push_back (0);
	t.second.push_word_back (1);
	t.first.push_back (8 * sizeof (Word));
	t.second.push_word_back (1);

	VD.add (w, u, v);

	std::cout << "Test 2:" << std::endl;
	std::cout << "    u = ";
	VD.write (std::cout, u) << std::endl;
	std::cout << "    v = ";
	VD.write (std::cout, v) << std::endl;
	std::cout << "u + v = ";
	VD.write (std::cout, w) << std::endl;
	std::cout << "    t = ";
	VD.write (std::cout, t) << std::endl;

	if (VD.areEqual (w, t))
		std::cout << "Test 2 okay" << std::endl;
	else
		std::cout << "Test 2 NOT okay" << std::endl;

	// Test 3: Indices in different blocks, v before u
	u.first.clear ();
	u.second.clear ();
	v.first.clear ();
	v.second.clear ();
	t.first.clear ();
	t.second.clear ();

	u.first.push_back (8 * sizeof (Word));
	v.first.push_back (0);

	u.second.push_word_back (1);
	v.second.push_word_back (1);

	t.first.push_back (0);
	t.second.push_word_back (1);
	t.first.push_back (8 * sizeof (Word));
	t.second.push_word_back (1);

	VD.add (w, u, v);

	std::cout << "Test 3:" << std::endl;
	std::cout << "    u = ";
	VD.write (std::cout, u) << std::endl;
	std::cout << "    v = ";
	VD.write (std::cout, v) << std::endl;
	std::cout << "u + v = ";
	VD.write (std::cout, w) << std::endl;
	std::cout << "    t = ";
	VD.write (std::cout, t) << std::endl;

	if (VD.areEqual (w, t))
		std::cout << "Test 3 okay" << std::endl;
	else
		std::cout << "Test 3 NOT okay" << std::endl;

	// Test 4: u = 0
	u.first.clear ();
	u.second.clear ();
	v.first.clear ();
	v.second.clear ();
	t.first.clear ();
	t.second.clear ();

	v.first.push_back (0);
	v.second.push_word_back (1);

	t.first.push_back (0);
	t.second.push_word_back (1);

	VD.add (w, u, v);

	std::cout << "Test 4:" << std::endl;
	std::cout << "    u = ";
	VD.write (std::cout, u) << std::endl;
	std::cout << "    v = ";
	VD.write (std::cout, v) << std::endl;
	std::cout << "u + v = ";
	VD.write (std::cout, w) << std::endl;
	std::cout << "    t = ";
	VD.write (std::cout, t) << std::endl;

	if (VD.areEqual (w, t))
		std::cout << "Test 4 okay" << std::endl;
	else
		std::cout << "Test 4 NOT okay" << std::endl;

	// Test 5: v = 0
	u.first.clear ();
	u.second.clear ();
	v.first.clear ();
	v.second.clear ();
	t.first.clear ();
	t.second.clear ();

	u.first.push_back (0);
	u.second.push_word_back (1);

	t.first.push_back (0);
	t.second.push_word_back (1);

	VD.add (w, u, v);

	std::cout << "Test 5:" << std::endl;
	std::cout << "    u = ";
	VD.write (std::cout, u) << std::endl;
	std::cout << "    v = ";
	VD.write (std::cout, v) << std::endl;
	std::cout << "u + v = ";
	VD.write (std::cout, w) << std::endl;
	std::cout << "    t = ";
	VD.write (std::cout, t) << std::endl;

	if (VD.areEqual (w, t))
		std::cout << "Test 5 okay" << std::endl;
	else
		std::cout << "Test 5 NOT okay" << std::endl;

	// Test 6: u, v have lengths longer than 1
	u.first.clear ();
	u.second.clear ();
	v.first.clear ();
	v.second.clear ();
	t.first.clear ();
	t.second.clear ();

	u.first.push_back (0);
	u.second.push_word_back (2ULL);
	u.first.push_back (8 * sizeof (Word));
	u.second.push_word_back (1ULL);

	v.first.push_back (0);
	v.second.push_word_back (4ULL);
	v.first.push_back (8 * sizeof (Word));
	v.second.push_word_back (1ULL);

	t.first.push_back (0);
	t.second.push_word_back (2 | 4);
//	t.first.push_back (8 * sizeof (Word));
//	t.second.push_word_back (0);

	VD.add (w, u, v);

	std::cout << "Test 4:" << std::endl;
	std::cout << "    u = ";
	VD.write (std::cout, u) << std::endl;
	std::cout << "    v = ";
	VD.write (std::cout, v) << std::endl;
	std::cout << "u + v = ";
	VD.write (std::cout, w) << std::endl;
	std::cout << "    t = ";
	VD.write (std::cout, t) << std::endl;

	if (VD.areEqual (w, t))
		std::cout << "Test 6 okay" << std::endl;
	else
		std::cout << "Test 6 NOT okay" << std::endl;

	std::cout << "Finished testing VD.add" << std::endl;
}

void testFirstNonzeroEntry ()
{
	std::cout << __FUNCTION__ << ": Testing VD.firstNonzeroEntry" << std::endl;

	Field F (2);
	VectorDomain<Field> VD (F);

	bool a;

	typedef HybridVector::second_type::word_iterator::value_type Word;

	HybridVector v;

	int idx;

	v.first.push_back (2 * __LINBOX_BITSOF_LONG);
	v.second.push_word_back (Endianness::e_0);

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 2 * __LINBOX_BITSOF_LONG)
		std::cout << "Test 1 okay" << std::endl;
	else
		std::cout << "Test 1 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 2 * 8 * sizeof (Word) << ")" << std::endl;

	v.first.clear ();
	v.second.clear ();

	v.first.push_back (2 * __LINBOX_BITSOF_LONG);
	v.second.push_word_back (Endianness::e_j (16));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 2 * __LINBOX_BITSOF_LONG + 16)
		std::cout << "Test 2 okay" << std::endl;
	else
		std::cout << "Test 2 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 2 * 8 * sizeof (Word) + 16 << ")" << std::endl;

	v.first.clear ();
	v.second.clear ();

	v.first.push_back (3 * __LINBOX_BITSOF_LONG);
	v.second.push_word_back (Endianness::e_j (__LINBOX_BITSOF_LONG - 1));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 3 * __LINBOX_BITSOF_LONG + (__LINBOX_BITSOF_LONG - 1))
		std::cout << "Test 3 okay" << std::endl;
	else
		std::cout << "Test 3 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 3 * 8 * sizeof (Word) + (8 * sizeof (Word) - 1) << ")" << std::endl;

	std::cout << "Finished testing VD.firstNonzeroEntry" << std::endl;
}

int main (int argc, char **argv)
{
	testAdd ();
	testFirstNonzeroEntry ();
}
