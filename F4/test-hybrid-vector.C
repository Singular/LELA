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
#include "linbox/vector/sparse-subvector-hybrid.h"

using namespace LinBox;

typedef GF2 Field;

bool testAdd ()
{
	Field F (2);
	VectorDomain<Field> VD (F);

	Vector<GF2>::Hybrid u, v, w, t;

	std::cout << "Testing VD.add..." << std::endl;

	// Test 1: Indices in same block
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	v.push_back (Vector<GF2>::Hybrid::value_type (0, 1 | (1 << 1)));
	t.push_back (Vector<GF2>::Hybrid::value_type (0, 1 << 1));

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
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, 1));
	t.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	t.push_back (Vector<GF2>::Hybrid::value_type (1, 1));

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
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (Vector<GF2>::Hybrid::value_type (1, 1));
	v.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	t.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	t.push_back (Vector<GF2>::Hybrid::value_type (1, 1));

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
	u.clear ();
	v.clear ();
	t.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	t.push_back (Vector<GF2>::Hybrid::value_type (0, 1));

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
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (Vector<GF2>::Hybrid::value_type (0, 1));
	t.push_back (Vector<GF2>::Hybrid::value_type (0, 1));

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
	u.clear ();
	v.clear ();
	t.clear ();

	u.push_back (Vector<GF2>::Hybrid::value_type (0, 2ULL));
	u.push_back (Vector<GF2>::Hybrid::value_type (1, 1ULL));

	v.push_back (Vector<GF2>::Hybrid::value_type (0, 4ULL));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, 1ULL));

	t.push_back (Vector<GF2>::Hybrid::value_type (0, 2 | 4));

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

	Vector<GF2>::Hybrid v;

	int idx;

	v.push_back (Vector<GF2>::Hybrid::value_type (2, Vector<GF2>::Hybrid::Endianness::e_0));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 2 * WordTraits<Vector<GF2>::Hybrid::word_type>::bits)
		std::cout << "Test 1 okay" << std::endl;
	else
		std::cout << "Test 1 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 2 * 8 * sizeof (Vector<GF2>::Hybrid::word_type) << ")" << std::endl;

	v.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (2, Vector<GF2>::Hybrid::Endianness::e_j (WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 2)));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 5 * WordTraits<Vector<GF2>::Hybrid::word_type>::bits / 2)
		std::cout << "Test 2 okay" << std::endl;
	else
		std::cout << "Test 2 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 2 * 8 * sizeof (Vector<GF2>::Hybrid::word_type) + 16 << ")" << std::endl;

	v.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (3, Vector<GF2>::Hybrid::Endianness::e_j (WordTraits<Vector<GF2>::Hybrid::word_type>::bits - 1)));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 3 * WordTraits<Vector<GF2>::Hybrid::word_type>::bits + (WordTraits<Vector<GF2>::Hybrid::word_type>::bits - 1))
		std::cout << "Test 3 okay" << std::endl;
	else
		std::cout << "Test 3 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 3 * 8 * sizeof (Vector<GF2>::Hybrid::word_type) + (8 * sizeof (Vector<GF2>::Hybrid::word_type) - 1) << ")" << std::endl;

	std::cout << "Finished testing VD.firstNonzeroEntry" << std::endl;
}

Vector<GF2>::Hybrid::word_type connect (Vector<GF2>::Hybrid::word_type word1, Vector<GF2>::Hybrid::word_type word2, int shift) 
{
	if (shift == 0)
		return word1;
	else
		return Vector<GF2>::Hybrid::Endianness::shift_left (word1, shift) | Vector<GF2>::Hybrid::Endianness::shift_right (word2, WordTraits<Vector<GF2>::Hybrid::word_type>::bits - shift);
}

void testSparseSubvectorHybrid ()
{
	std::cout << __FUNCTION__ << ": Testing SparseSubvector<Hybrid>..." << std::endl;

	Vector<GF2>::Hybrid v;

	Vector<GF2>::Hybrid::word_type pattern[2] = { 0xf0f0f0f0f0f0f0f0ULL, 0xaaaaaaaaaaaaaaaaULL };
	Vector<GF2>::Hybrid::word_type check;

	uint16 offset = 5;
	uint16 len = WordTraits<Vector<GF2>::Hybrid::word_type>::bits + 4;

	std::cout << __FUNCTION__ << ": Test 1" << std::endl;

	v.push_back (Vector<GF2>::Hybrid::value_type (0, pattern[0]));
	v.push_back (Vector<GF2>::Hybrid::value_type (1, pattern[1]));

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v1 (v, offset, offset + len);

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag>::const_iterator i_v1 = v1.begin ();

	if (i_v1->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << std::dec << i_v1->first << std::endl;

	check = connect (pattern[0], pattern[1], offset);

	if (i_v1->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v1->second) << std::endl;

	++i_v1;

	if (i_v1->first != 1)
		std::cerr << __FUNCTION__ << ": ERROR: Second index should be 1, is " << std::dec << i_v1->first << std::endl;

	check = connect (pattern[1], 0ULL, offset) & Vector<GF2>::Hybrid::Endianness::mask_left (len % WordTraits<Vector<GF2>::Hybrid::word_type>::bits);

	if (i_v1->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Second word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v1->second) << std::endl;

	++i_v1;

	if (i_v1 != v1.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;

	v.clear ();

	std::cout << __FUNCTION__ << ": Test 2" << std::endl;

	v.push_back (Vector<GF2>::Hybrid::value_type (1, pattern[0]));

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v2 (v, offset, offset + len);

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag>::const_iterator i_v2 = v2.begin ();

	if (i_v2->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << std::dec << i_v2->first << std::endl;

	check = connect (0ULL, pattern[0], offset);

	if (i_v2->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v2->second) << std::endl;

	++i_v2;

	if (i_v2 == v2.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Iterator ended prematurely." << std::endl;

	if (i_v2->first != 1)
		std::cerr << __FUNCTION__ << ": ERROR: Second index should be 1, is " << std::dec << i_v2->first << std::endl;

	check = connect (pattern[0], 0ULL, offset) & Vector<GF2>::Hybrid::Endianness::mask_left (len % WordTraits<Vector<GF2>::Hybrid::word_type>::bits);

	if (i_v2->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Second word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v2->second) << std::endl;

	++i_v2;

	if (i_v2 != v2.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;

	std::cout << __FUNCTION__ << ": Test 3" << std::endl;

	v.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (0, pattern[0]));
	v.push_back (Vector<GF2>::Hybrid::value_type (2, pattern[1]));

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v3 (v, offset, WordTraits<Vector<GF2>::Hybrid::word_type>::bits + offset + len);

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag>::const_iterator i_v3 = v3.begin ();

	if (i_v3->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << std::dec << i_v3->first << std::endl;

	check = connect (pattern[0], 0ULL, offset);

	if (i_v3->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v3->second) << std::endl;

	++i_v3;

	if (i_v3 == v3.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Iterator ended prematurely." << std::endl;

	if (i_v3->first != 1)
		std::cerr << __FUNCTION__ << ": ERROR: Second index should be 1, is " << std::dec << i_v3->first << std::endl;

	check = connect (0ULL, pattern[1], offset);

	if (i_v3->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Second word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v3->second) << std::endl;

	++i_v3;

	if (i_v3 == v3.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Iterator ended prematurely." << std::endl;

	if (i_v3->first != 2)
		std::cerr << __FUNCTION__ << ": ERROR: Third index should be 2, is " << std::dec << i_v3->first << std::endl;

	check = connect (pattern[1], 0ULL, offset) & Vector<GF2>::Hybrid::Endianness::mask_left (len % WordTraits<Vector<GF2>::Hybrid::word_type>::bits);

	if (i_v3->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Third word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v3->second) << std::endl;

	++i_v3;

	if (i_v3 != v3.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;

	std::cout << __FUNCTION__ << ": Test 4" << std::endl;

	v.clear ();

	v.push_back (Vector<GF2>::Hybrid::value_type (0, pattern[0]));

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v4 (v, offset, WordTraits<Vector<GF2>::Hybrid::word_type>::bits + offset + len);

	SparseSubvector<Vector<GF2>::Hybrid, VectorCategories::HybridZeroOneVectorTag>::const_iterator i_v4 = v4.begin ();

	if (i_v4->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << std::dec << i_v4->first << std::endl;

	check = connect (pattern[0], 0ULL, offset);

	if (i_v4->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << static_cast<uint64> (check) << ", is " << std::hex << static_cast<uint64> (i_v4->second) << std::endl;

	++i_v4;

	if (i_v4 != v4.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;
	
	std::cout << __FUNCTION__ << ": done" << std::endl;
}

int main (int argc, char **argv)
{
	testAdd ();
	testFirstNonzeroEntry ();
	testSparseSubvectorHybrid ();
}
