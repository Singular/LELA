/* linbox/vector/vector-domain-gf2.tcc
 * Copyright 2003, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_VECTOR_DOMAIN_GF2_INL
#define __LINBOX_VECTOR_DOMAIN_GF2_INL

namespace LinBox 
{ 

template <class Vector1, class Vector2>
inline bool &DotProductDomain<GF2>::dotSpecializedDD
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	linbox_check (v1.size () == v2.size ());

	typename Vector1::word_type t = 0;
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_word_iterator j = v2.wordBegin ();

	while (i != v1.wordEnd () - 1)
		t ^= *i++ & *j++;
        
        const size_t zeroing = WordTraits<typename Vector1::word_type>::bits - (v1.size() % WordTraits<typename Vector1::word_type>::bits);
        typename Vector1::word_type lastdot = *i & *j;
        lastdot <<= zeroing;
        lastdot >>= zeroing;
        
        t ^= lastdot;
        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class Vector1, class Vector2>
inline bool &DotProductDomain<GF2>::dotSpecializedDSP
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	typename Vector2::const_iterator i;

	res = 0;

	for (i = v2.begin (); i != v2.end (); ++i)
		res ^= v1[*i];

	return res;
}

template <class Vector1, class Vector2>
inline bool &DotProductDomain<GF2>::dotSpecializedDH
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	typename Vector1::word_type t = 0;
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::first_type::const_iterator j_idx = v2.first.begin ();
	typename Vector2::second_type::const_word_iterator j_elt = v2.second.wordBegin ();

	for (; j_idx != v2.first.end (); ++j_idx, ++j_elt)
		t ^= *(i + *j_idx) & *j_elt;
        
        return res = WordTraits<typename Vector1::word_type>::ParallelParity (t);
}

template <class Iterator, class Endianness, class Vector1, class Vector2>
inline BitVectorReference<Iterator, Endianness> DotProductDomain<GF2>::dotSpecializedDD
	(BitVectorReference<Iterator, Endianness> res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	bool tmp;

	return res = dotSpecializedDD (tmp, v1, v2);
}

template <class Iterator, class Endianness, class Vector1, class Vector2>
inline BitVectorReference<Iterator, Endianness> DotProductDomain<GF2>::dotSpecializedDSP
	(BitVectorReference<Iterator, Endianness> res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	typename Vector2::const_iterator i;

	res = 0;

	for (i = v2.begin (); i != v2.end (); ++i)
		res ^= v1[*i];

	return res;
}

template <class Iterator, class Endianness, class Vector1, class Vector2>
inline BitVectorReference<Iterator, Endianness> DotProductDomain<GF2>::dotSpecializedDH
	(BitVectorReference<Iterator, Endianness> res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	bool tmp;

	return res = dotSpecializedDH (tmp, v1, v2);
}

template <class Vector>
std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
						   VectorCategories::DenseZeroOneVectorTag) const
{
	os << "[ ";

	for (typename Vector::const_iterator i = x.begin (); i != x.end (); ++i)
		os << *i << ' ';

	os << ']';

#if 0
	os << "( ";

	for (typename Vector::const_word_iterator i = x.wordBegin (); i != x.wordEnd (); ++i)
		os << *i << ' ';

	os << ')';
#endif

	return os;
}

template <class Vector>
std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
						   VectorCategories::SparseZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type index_type;

	typename Vector::const_iterator i;
	index_type idx = 0;

	os << "[ ";

	for (i = x.begin (); i != x.end (); ++i) {
		while (++idx <= *i)
			os << 0 << ' ';

		os << 1 << ' ';
	}
	os << ']';

	return os;
}

template <class Vector, class Endianness>
std::ostream &VectorDomain<GF2>::writeHybridSpecialized (std::ostream &os, const Vector &x, Endianness) const
{
	typedef typename std::iterator_traits<typename Vector::first_type::const_iterator>::value_type index_type;
	typedef typename std::iterator_traits<typename Vector::second_type::const_word_iterator>::value_type word_type;

	typename Vector::first_type::const_iterator i_idx;
	typename Vector::second_type::const_word_iterator i_elt;
	index_type idx = 0;
	word_type mask;

	os << "[ ";

	for (i_idx = x.first.begin (), i_elt = x.second.wordBegin (); i_idx != x.first.end (); ++i_idx, ++i_elt) {
		while (++idx <= *i_idx << WordTraits<word_type>::logof_size)
			os << "0 ";

		for (mask = Endianness::e_0; mask != 0; mask = Endianness::shift_right (mask, 1)) {
			if (*i_elt & mask)
				os << "1 ";
			else
				os << "0 ";
			++idx;
		}
	}
	os << ']';

	return os;
}

template <class Vector>
std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
						  VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector::iterator i;
	char c;

	do { is >> c ; } while (!std::isdigit (c));

	is.unget ();

	for (i = x.begin (); i != x.end (); ++i)
		is >> *i;

	return is;
}

template <class Vector>
std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
						  VectorCategories::SparseZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type index_type;

	char c;
	index_type idx;

	do { is >> c ; } while (!std::isdigit (c));

	is.unget ();
	x.clear ();

	while (1) {
		is >> c;

		if (!std::isdigit (c) && c != ' ') break;
		is.unget ();
		is >> idx;
		x.push_back (idx);
	}

	return is;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector2::const_iterator>::value_type index_type;

	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();
	index_type idx = 0;

	for (; j != v2.end (); ++j, ++i, ++idx) {
		while (idx < *j) {
			if (*i) return false;
			++idx;
			++i;
		}

		if (!*i) return false;
	}

	for (; i != v1.end (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_word_iterator j = v2.wordBegin ();

	for (; j != v2.wordEnd (); ++j, ++i)
		if (*i != *j) return false;

	return true;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::SparseZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector1::const_iterator i_1;
	typename Vector2::const_iterator i_2;

	if (v1.size () != v2.size ())
		return false;

	for (i_1 = v1.begin (), i_2 = v2.begin (); i_1 != v1.end (); ++i_1, ++i_2)
		if (*i_1 != *i_2)
			return false;

	return true;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::HybridZeroOneVectorTag,
					     VectorCategories::HybridZeroOneVectorTag) const
{
	typename Vector1::first_type::const_iterator i_1_idx;
	typename Vector2::first_type::const_iterator i_2_idx;
	typename Vector1::second_type::const_iterator i_1_elt;
	typename Vector2::second_type::const_iterator i_2_elt;

	if (v1.first.size () != v2.first.size ())
		return false;

	for (i_1_idx = v1.first.begin (), i_2_idx = v2.first.begin (); i_1_idx != v1.first.end (); ++i_1_idx, ++i_2_idx)
		if (*i_1_idx != *i_2_idx)
			return false;

	for (i_1_elt = v1.second.begin (), i_2_elt = v2.second.begin (); i_1_elt != v1.second.end (); ++i_1_elt, ++i_2_elt)
		if (*i_1_elt != *i_2_elt)
			return false;

	return true;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::HybridZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector2::first_type::const_iterator>::value_type index_type;

	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::first_type::const_iterator j_idx = v2.first.begin ();
	typename Vector2::second_type::const_word_iterator j_elt = v2.second.wordBegin ();
	index_type idx = 0;

	for (; j_idx != v2.first.end (); ++j_idx, ++j_elt) {
		while (idx < *j_idx) {
			if (*i) return false;
			++idx;
			++i;
		}

		if (*i++ != *j_elt)
			return false;

		++idx;
	}

	for (; i != v1.wordEnd (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::HybridZeroOneSequenceVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector2::const_iterator>::value_type::first_type index_type;

	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_iterator j = v2.begin ();
	index_type idx = 0;
	int t;

	for (; j != v2.end (); ++j) {
		while (idx < j->first) {
			if (*i) return false;
			++idx;
			++i;
		}

		if (*i++ != j->second)
			return false;

		++idx;
	}

	for (; i != v1.wordEnd (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector>
bool VectorDomain<GF2>::isZeroSpecialized (const Vector &v,
					   VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector::const_word_iterator i;

	for (i = v.wordBegin (); i != v.wordEnd (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::SparseZeroOneVectorTag,
					     VectorCategories::DenseZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector1::const_iterator>::value_type index_type;

	typename Vector2::const_iterator i;
	index_type idx = 0;

	res.clear ();

	for (i = v.begin (); i != v.end (); ++i, ++idx)
		if (*i) res.push_back (idx);

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::HybridZeroOneVectorTag,
					     VectorCategories::DenseZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector1::first_type::const_iterator>::value_type index_type;

	typename Vector2::const_word_iterator i;
	index_type idx = 0;

	res.first.clear ();
	res.second.clear ();

	for (i = v.wordBegin (); i != v.wordEnd (); ++i, ++idx) {
		if (*i) {
			res.first.push_back (idx);
			res.second.push_word_back (*i);
		}
	}

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i;

	std::fill (res.wordBegin (), res.wordEnd (), 0);

	for (i = v.begin (); i != v.end (); ++i)
        	res[*i] = true;

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::HybridZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector2::second_type::const_word_iterator>::value_type word_type;

	typename Vector2::first_type::const_iterator i_idx;
	typename Vector2::second_type::const_word_iterator i_elt;

	std::fill (res.wordBegin (), res.wordEnd (), 0);

	for (i_idx = v.first.begin (), i_elt = v.second.wordBegin (); i_idx != v.first.end (); ++i_idx, ++i_elt)
        	*(res.wordBegin () + *i_idx) = *i_elt;

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::HybridZeroOneSequenceVectorTag) const
{
	typename Vector2::const_iterator i;

	std::fill (res.wordBegin (), res.wordEnd (), 0);

	for (i = v.begin (); i != v.end (); ++i)
        	*(res.wordBegin () + i->first) = i->second;

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::SparseZeroOneVectorTag,
					     VectorCategories::HybridZeroOneVectorTag) const
{
	typedef typename Vector2::second_type::Endianness Endianness;

	typedef typename std::iterator_traits<typename Vector1::const_iterator>::value_type index_type;
	typedef typename std::iterator_traits<typename Vector2::second_type::const_word_iterator>::value_type word_type;

	typename Vector2::first_type::const_iterator i_idx;
	typename Vector2::second_type::const_word_iterator i_elt;
	index_type idx = 0;
	word_type t;

	res.clear ();

	for (i_idx = v.first.begin (), i_elt = v.second.wordBegin (); i_idx != v.first.end (); ++i_idx, ++i_elt)
		for (t = Endianness::e_0, idx = *i_idx << WordTraits<word_type>::logof_size; t != 0; t = Endianness::shift_right (t, 1), ++idx)
			if (*i_elt & t) res.push_back (idx);

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::HybridZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
{
	typedef typename Vector1::second_type::Endianness Endianness;

	typedef typename std::iterator_traits<typename Vector2::const_iterator>::value_type index_type;
	typedef typename std::iterator_traits<typename Vector1::second_type::const_word_iterator>::value_type word_type;

	typename Vector2::const_iterator i;

	res.first.clear ();
	res.second.clear ();

	for (i = v.begin (); i != v.end (); ++i) {
		if (res.first.empty () || res.first.back () != *i >> WordTraits<word_type>::logof_size) {
			res.first.push_back (*i >> WordTraits<word_type>::logof_size);
			res.second.push_word_back (0);
		}

		res.second.back_word () |= Endianness::e_j (*i & WordTraits<word_type>::pos_mask);
	}

	return res;
}

template <class Vector, class Iterator>
inline Vector &VectorDomain<GF2>::permuteSpecialized (Vector &v, Iterator P_start, Iterator P_end,
						      VectorCategories::DenseZeroOneVectorTag) const 
{
	Iterator i;
	bool t;

	for (i = P_start; i != P_end; ++i) {
		t = v[i->first];
		v[i->first] = v[i->second];
		v[i->second] = t;
	}

	return v;
}

template <class Vector, class Iterator>
inline Vector &VectorDomain<GF2>::permuteSpecialized (Vector &v, Iterator P_start, Iterator P_end,
						      VectorCategories::SparseZeroOneVectorTag) const 
{
	typename Vector::iterator j;

	for (j = v.begin (); j != v.end (); ++j)
		*j = permutationImage (*j, P_start, P_end);

	std::sort (v.begin (), v.end ());

	return v;
}

template <class Vector, class Iterator>
inline Vector &VectorDomain<GF2>::permuteSpecialized (Vector &v, Iterator P_start, Iterator P_end,
						      VectorCategories::HybridZeroOneVectorTag) const 
{
	::LinBox::Vector<GF2>::Sparse w;

	copy (w, v);
	permute (w, P_start, P_end);
	copy (v, w);

	return v;
}

template <class Vector1, class Vector2>
bool &VectorDomain<GF2>::dotSpecialized (bool &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();
	res = false;

	while (i != v1.end () || j != v2.end ()) {
		while (i != v1.end () && (j == v2.end () || *i < *j)) { res = !res; ++i; }
		while (j != v2.end () && (i == v1.end () || *j < *i)) { res = !res; ++j; }
		if (i != v1.end () && j != v2.end () && *i == *j) { ++i; ++j; }
	}

	return res;
}

template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::DenseZeroOneVectorTag,
					    VectorCategories::DenseZeroOneVectorTag,
					    VectorCategories::DenseZeroOneVectorTag) const
{
	linbox_check (res.size () == y.size ());
	linbox_check (res.size () == x.size ());

	typename Vector1::word_iterator i = res.wordBegin ();
	typename Vector2::const_word_iterator j = y.wordBegin ();
	typename Vector3::const_word_iterator k = x.wordBegin ();

	for (; i != res.wordEnd (); ++i)
		*i = *j++ ^ *k++;

	return res;
}

template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::SparseZeroOneVectorTag,
					    VectorCategories::SparseZeroOneVectorTag,
					    VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i = y.begin ();
	typename Vector3::const_iterator j = x.begin ();

	res.clear ();

	while (i != y.end () || j != x.end ()) {
		while (i != y.end () && (j == x.end () || *i < *j)) { res.push_back (*i); ++i; }
		while (j != x.end () && (i == y.end () || *j < *i)) { res.push_back (*j); ++j; }
		if (i != y.end () && j != x.end () && *i == *j) { ++i; ++j; }
	}

	return res;
}

template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::HybridZeroOneVectorTag,
					    VectorCategories::HybridZeroOneVectorTag,
					    VectorCategories::HybridZeroOneVectorTag) const
{
	typename Vector2::first_type::const_iterator i_idx = y.first.begin ();
	typename Vector3::first_type::const_iterator j_idx = x.first.begin ();

	typename Vector2::second_type::const_word_iterator i_elt = y.second.wordBegin ();
	typename Vector3::second_type::const_word_iterator j_elt = x.second.wordBegin ();

	res.first.clear ();
	res.second.clear ();

	while (i_idx != y.first.end () || j_idx != x.first.end ()) {
		while (i_idx != y.first.end () && (j_idx == x.first.end () || *i_idx < *j_idx)) {
			res.first.push_back (*i_idx);
			res.second.push_word_back (*i_elt);
			++i_idx;
			++i_elt;
		}
		while (j_idx != x.first.end () && (i_idx == y.first.end () || *j_idx < *i_idx)) {
			res.first.push_back (*j_idx);
			res.second.push_word_back (*j_elt);
			++j_idx;
			++j_elt;
		}
		if (i_idx != y.first.end () && j_idx != x.first.end () && *i_idx == *j_idx) {
			if (*i_elt ^ *j_elt) {
				res.first.push_back (*i_idx);
				res.second.push_word_back (*i_elt ^ *j_elt);
			}

			++i_idx;
			++j_idx;
			++i_elt;
			++j_elt;
		}
	}

	res.second.fix_size (x.second.size ());

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag,
					      VectorCategories::DenseZeroOneVectorTag) const
{
	linbox_check (y.size () == x.size ());

	typename Vector1::word_iterator i = y.wordBegin ();
	typename Vector2::const_word_iterator j = x.wordBegin ();

	for (; i != y.wordEnd (); ++i, ++j)
		*i ^= *j;

	return y;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag,
					      VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i;

	for (i = x.begin (); i != x.end (); ++i)
		y[*i] = !y[*i];

	return y;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag,
					      VectorCategories::HybridZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector2::second_type::const_word_iterator>::value_type word_type;

	typename Vector2::first_type::const_iterator i_idx;
	typename Vector2::second_type::const_word_iterator i_elt;
	typename Vector1::word_iterator j = y.wordBegin ();

	for (i_idx = x.first.begin (), i_elt = x.second.wordBegin (); i_idx != x.first.end (); ++i_idx, ++i_elt)
		*(j + *i_idx) ^= *i_elt;

	return y;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag,
					      VectorCategories::HybridZeroOneSequenceVectorTag) const
{
	typename Vector2::const_iterator i;
	typename Vector1::word_iterator j = y.wordBegin ();

	for (i = x.begin (); i != x.end (); ++i)
		*(j + i->first) ^= i->second;

	return y;
}

template <class word, class Endianness>
inline int VectorDomain<GF2>::firstNonzeroEntryInWord (word v) const
{
	// FIXME: This can be made faster...
	word mask = Endianness::e_0;
	int idx = 0;

	for (; mask != 0; mask = Endianness::shift_right (mask, 1), ++idx)
		if (v & mask)
			return idx;

	return -1;
}

template <class Vector>
inline int VectorDomain<GF2>::firstNonzeroEntrySpecialized (bool &a, const Vector &v,
							    VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector::const_iterator i;

	for (i = v.begin (); i != v.end (); ++i)
		if (*i)
			return i - v.begin ();

	return -1;
}

template <class Vector>
inline int VectorDomain<GF2>::firstNonzeroEntrySpecialized (bool &a, const Vector &v,
							    VectorCategories::SparseZeroOneVectorTag) const
{
	if (v.empty ())
		return -1;
	else
		return v.front ();
}

template <class Vector>
inline int VectorDomain<GF2>::firstNonzeroEntrySpecialized (bool &a, const Vector &v,
							    VectorCategories::HybridZeroOneVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector::second_type::const_word_iterator>::value_type word_type;

	if (v.first.empty ())
		return -1;
	else
		return firstNonzeroEntryInWord<word_type, typename Vector::second_type::Endianness> (v.second.front_word ())
			+ (v.first.front () << WordTraits<word_type>::logof_size);
}

template <class Vector>
inline int VectorDomain<GF2>::firstNonzeroEntrySpecialized (bool &a, const Vector &v,
							    VectorCategories::HybridZeroOneSequenceVectorTag) const
{
	typedef typename std::iterator_traits<typename Vector::const_iterator>::value_type::second_type word_type;

	if (v.empty ())
		return -1;
	else
		return firstNonzeroEntryInWord<word_type, typename Vector::second_type::Endianness> (v.begin ()->second)
			+ (v.begin ()->first << WordTraits<word_type>::logof_size);
}

} // namespace LinBox

#endif // __LINBOX_VECTOR_DOMAIN_GF2_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
