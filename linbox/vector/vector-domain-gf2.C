/* linbox/vector/vector-domain-gf2.C
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include <iostream>

#include "linbox/field/gf2.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/bit-subvector.h"
#include "linbox/vector/bit-subvector-word-aligned.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/vector-domain-gf2.tcc"
#include "linbox/vector/sparse-subvector.h"

using namespace LinBox;

typedef GF2 Field;

template <class Vector1, class Vector2, class Vector3>
static void instantiate_ternary (const Field &F, Vector1 &v1, Vector2 &v2, Vector3 &v3)
{
	VectorDomain<Field> VD (F);

	Field::Element a;

	VD.add (v1, v2, v3);
	VD.sub (v1, v2, v3);
	VD.axpy (v1, a, v2, v3);
}

template <class Vector1, class Vector2>
static void add_two_and_instantiate (const Field &F, Vector1 &v, Vector2 &w)
{
	Vector<Field>::Dense v1;
	Vector<Field>::Sparse v2;
	Vector<Field>::Hybrid v3;

	BitSubvector<BitVector<>::iterator, BitVector<>::const_iterator> v4;
	BitSubvectorWordAligned<BitVector<>::word_iterator, BitVector<>::const_word_iterator, DefaultEndianness> v5;
//	SparseSubvector<Vector<Field>::Sparse, VectorCategories::SparseZeroOneVectorTag> v6;
//	SparseSubvector<Vector<Field>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v7;

	instantiate_ternary (F, v, w, v);
	instantiate_ternary (F, v, w, w);
	instantiate_ternary (F, v, w, v1);
	instantiate_ternary (F, v, w, v2);
	instantiate_ternary (F, v, w, v3);
	instantiate_ternary (F, v, w, v4);
	instantiate_ternary (F, v, w, v5);
//	instantiate_ternary (F, v, w, v6);
//	instantiate_ternary (F, v, w, v7);
}

template <class Vector1, class Vector2>
static void instantiate_binary (const Field &F, Vector1 &v1, Vector2 &v2)
{
	VectorDomain<Field> VD (F);

	Field::Element a;

	BitVectorReference<BitVector<>::word_iterator, BitVector<>::Endianness> b1 (BitVector<>::word_iterator (), 0);
	BitVectorReference<uint64 *, LittleEndian<uint64> > b2 (NULL, 0);
	BitVectorReference<uint64 *, BigEndian<uint64> > b3 (NULL, 0);

	VD.copy (v1, v2);
	VD.areEqual (v1, v2);
	VD.dot (a, v1, v2);
	VD.dot (b1, v1, v2);
	VD.dot (b2, v1, v2);
	VD.dot (b3, v1, v2);
	VD.addin (v1, v2);
	VD.subin (v1, v2);
	VD.mul (v1, v2, a);
	VD.neg (v1, v2);
	VD.axpyin (v1, a, v2);

	add_two_and_instantiate (F, v1, v2);
}

template <class Vector1, class T>
static void instantiate_permute (const Field &F, Vector1 &v)
{
	VectorDomain<Field> VD (F);

	std::vector<std::pair<T, T> > P1;
	const std::vector<std::pair<T, T> > P2;
	std::list<std::pair<T, T> > P3;
	const std::list<std::pair<T, T> > P4;
	std::pair<T, T> *P5;
	const std::pair<T, T> *P6;

	VD.permute (v, P1.begin (), P1.end ());
	VD.permute (v, P2.begin (), P2.end ());
	VD.permute (v, P3.begin (), P3.end ());
	VD.permute (v, P4.begin (), P4.end ());
	VD.permute (v, P5, P5);
	VD.permute (v, P6, P6);
}

template <class Vector1>
static void add_one_and_instantiate (const Field &F, Vector1 &v)
{
	Vector<Field>::Dense v1;
	Vector<Field>::Sparse v2;
	Vector<Field>::Hybrid v3;

	BitSubvector<BitVector<>::iterator, BitVector<>::const_iterator> v4;
	BitSubvectorWordAligned<BitVector<>::word_iterator, BitVector<>::const_word_iterator, DefaultEndianness> v5;
//	SparseSubvector<Vector<Field>::Sparse, VectorCategories::SparseZeroOneVectorTag> v6;
//	SparseSubvector<Vector<Field>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v7;

	instantiate_binary (F, v, v);
	instantiate_binary (F, v, v1);
	instantiate_binary (F, v, v2);
	instantiate_binary (F, v, v3);
	instantiate_binary (F, v, v4);
	instantiate_binary (F, v, v5);
//	instantiate_binary (F, v, v6);
//	instantiate_binary (F, v, v7);
}

template <class Vector1>
static void instantiate_unary (const Field &F, Vector1 &v, std::ostream &os, std::istream &is)
{
	VectorDomain<Field> VD (F);

	Field::Element a;

//	VD.read (is, v);
	VD.write (os, v);
	VD.firstNonzeroEntry (a, v);
	VD.isZero (v);
	VD.negin (v);
	VD.mulin (v, a);
	VD.swap (v, v);

	instantiate_permute<Vector1, size_t> (F, v);
	instantiate_permute<Vector1, uint32> (F, v);
	instantiate_permute<Vector1, uint16> (F, v);

	add_one_and_instantiate (F, v);
}

static void instantiate_VectorDomain_GF2 (const Field &F, std::ostream &os, std::istream &is)
{
	Vector<Field>::Dense v1;
	Vector<Field>::Sparse v2;
	Vector<Field>::Hybrid v3;

	BitSubvector<BitVector<>::iterator, BitVector<>::const_iterator> v4;
	BitSubvectorWordAligned<BitVector<>::word_iterator, BitVector<>::const_word_iterator, DefaultEndianness> v5;
//	SparseSubvector<Vector<Field>::Sparse, VectorCategories::SparseZeroOneVectorTag> v6;
//	SparseSubvector<Vector<Field>::Hybrid, VectorCategories::HybridZeroOneVectorTag> v7;

	instantiate_unary (F, v1, os, is);
	instantiate_unary (F, v2, os, is);
	instantiate_unary (F, v3, os, is);
	instantiate_unary (F, v4, os, is);
	instantiate_unary (F, v5, os, is);
//	instantiate_unary (F, v6, os, is);
//	instantiate_unary (F, v7, os, is);
}
