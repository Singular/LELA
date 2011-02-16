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

	v.first.push_back (2 * WordTraits<Word>::bits);
	v.second.push_word_back (Endianness::e_0);

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 2 * WordTraits<Word>::bits)
		std::cout << "Test 1 okay" << std::endl;
	else
		std::cout << "Test 1 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 2 * 8 * sizeof (Word) << ")" << std::endl;

	v.first.clear ();
	v.second.clear ();

	v.first.push_back (2 * WordTraits<Word>::bits);
	v.second.push_word_back (Endianness::e_j (16));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 2 * WordTraits<Word>::bits + 16)
		std::cout << "Test 2 okay" << std::endl;
	else
		std::cout << "Test 2 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 2 * 8 * sizeof (Word) + 16 << ")" << std::endl;

	v.first.clear ();
	v.second.clear ();

	v.first.push_back (3 * WordTraits<Word>::bits);
	v.second.push_word_back (Endianness::e_j (WordTraits<Word>::bits - 1));

	idx = VD.firstNonzeroEntry (a, v);

	if (idx == 3 * WordTraits<Word>::bits + (WordTraits<Word>::bits - 1))
		std::cout << "Test 3 okay" << std::endl;
	else
		std::cout << "Test 3 not okay: VD.firstNonzeroEntry (a, v) = " << idx
			  << " (should be " << 3 * 8 * sizeof (Word) + (8 * sizeof (Word) - 1) << ")" << std::endl;

	std::cout << "Finished testing VD.firstNonzeroEntry" << std::endl;
}

word connect (word word1, word word2, int shift) 
{
	if (shift == 0)
		return word1;
	else
		return Endianness::shift_left (word1, shift) | Endianness::shift_right (word2, WordTraits<Word>::bits - shift);
}

void testSparseSubvectorHybrid ()
{
	std::cout << __FUNCTION__ << ": Testing SparseSubvector<Hybrid>..." << std::endl;

	HybridVector v;

	Word pattern[2] = { 0xf0f0f0f0f0f0f0f0ULL, 0xaaaaaaaaaaaaaaaaULL };
	Word check;

	uint16 offset = 10;
	uint16 len = WordTraits<Word>::bits + 8;

	std::cout << __FUNCTION__ << ": Test 1" << std::endl;

	v.first.push_back (0);
	v.second.push_word_back (pattern[0]);
	v.first.push_back (WordTraits<Word>::bits);
	v.second.push_word_back (pattern[1]);

	SparseSubvector<HybridVector> v1 (v, offset, offset + len);

	SparseSubvector<HybridVector>::const_iterator i_v1 = v1.begin ();

	if (i_v1->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << i_v1->first << std::endl;

	check = connect (pattern[0], pattern[1], offset);

	if (i_v1->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << check << ", is " << std::hex << i_v1->second << std::endl;

	++i_v1;

	if (i_v1->first != WordTraits<Word>::bits)
		std::cerr << __FUNCTION__ << ": ERROR: Second index should be " << WordTraits<Word>::bits << ", is " << i_v1->first << std::endl;

	check = connect (pattern[1], 0ULL, offset) & Endianness::mask_left (len % WordTraits<Word>::bits);

	if (i_v1->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Second word should be " << std::hex << check << ", is " << std::hex << i_v1->second << std::endl;

	++i_v1;

	if (i_v1 != v1.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;

	v.first.clear ();
	v.second.clear ();

	std::cout << __FUNCTION__ << ": Test 2" << std::endl;

	v.first.push_back (WordTraits<Word>::bits);
	v.second.push_word_back (pattern[0]);

	SparseSubvector<HybridVector> v2 (v, offset, offset + len);

	SparseSubvector<HybridVector>::const_iterator i_v2 = v2.begin ();

	if (i_v2->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << i_v2->first << std::endl;

	check = connect (0ULL, pattern[0], offset);

	if (i_v2->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << check << ", is " << std::hex << i_v2->second << std::endl;

	++i_v2;

	if (i_v2->first != WordTraits<Word>::bits)
		std::cerr << __FUNCTION__ << ": ERROR: Second index should be " << WordTraits<Word>::bits << ", is " << i_v2->first << std::endl;

	check = connect (pattern[0], 0ULL, offset) & Endianness::mask_left (len % WordTraits<Word>::bits);

	if (i_v2->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Second word should be " << std::hex << check << ", is " << std::hex << i_v2->second << std::endl;

	++i_v2;

	if (i_v2 != v2.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;

	std::cout << __FUNCTION__ << ": Test 3" << std::endl;

	v.first.clear ();
	v.second.clear ();

	v.first.push_back (0);
	v.second.push_word_back (pattern[0]);
	v.first.push_back (2 * WordTraits<Word>::bits);
	v.second.push_word_back (pattern[1]);

	SparseSubvector<HybridVector> v3 (v, offset, WordTraits<Word>::bits + offset + len);

	SparseSubvector<HybridVector>::const_iterator i_v3 = v3.begin ();

	if (i_v3->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << i_v3->first << std::endl;

	check = connect (pattern[0], 0ULL, offset);

	if (i_v3->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << check << ", is " << std::hex << i_v3->second << std::endl;

	++i_v3;

	if (i_v3->first != WordTraits<Word>::bits)
		std::cerr << __FUNCTION__ << ": ERROR: Second index should be " << WordTraits<Word>::bits << ", is " << i_v3->first << std::endl;

	check = connect (0ULL, pattern[1], offset);

	if (i_v3->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Second word should be " << std::hex << check << ", is " << std::hex << i_v3->second << std::endl;

	++i_v3;

	if (i_v3->first != 2 * WordTraits<Word>::bits)
		std::cerr << __FUNCTION__ << ": ERROR: Third index should be " << 2 * WordTraits<Word>::bits << ", is " << i_v3->first << std::endl;

	check = connect (pattern[1], 0ULL, offset) & Endianness::mask_left (len % WordTraits<Word>::bits);

	if (i_v3->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: Third word should be " << std::hex << check << ", is " << std::hex << i_v3->second << std::endl;

	++i_v3;

	if (i_v3 != v3.end ())
		std::cerr << __FUNCTION__ << ": ERROR: Did not hit end when expected." << std::endl;

	std::cout << __FUNCTION__ << ": Test 4" << std::endl;

	v.first.clear ();
	v.second.clear ();

	v.first.push_back (0);
	v.second.push_word_back (pattern[0]);

	SparseSubvector<HybridVector> v4 (v, offset, WordTraits<Word>::bits + offset + len);

	SparseSubvector<HybridVector>::const_iterator i_v4 = v4.begin ();

	if (i_v4->first != 0)
		std::cerr << __FUNCTION__ << ": ERROR: First index should be 0, is " << i_v4->first << std::endl;

	check = connect (pattern[0], 0ULL, offset);

	if (i_v4->second != check)
		std::cerr << __FUNCTION__ << ": ERROR: First word should be " << std::hex << check << ", is " << std::hex << i_v4->second << std::endl;

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
