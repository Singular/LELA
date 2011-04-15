/* linbox/vector/vector-domain-gf2.h
 * Copyright 2003, 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_VECTOR_DOMAIN_GF2_H
#define __LINBOX_VECTOR_DOMAIN_GF2_H

#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"

namespace LinBox 
{ 

// Specialization of DotProductDomain for GF2

template <>
class DotProductDomain<GF2> : private virtual VectorDomainBase<GF2>
{
    public:

	typedef bool Element;

	DotProductDomain (const GF2 &F)
		: VectorDomainBase<GF2> (F)
	{}

    protected:
	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const;

	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const;

	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDH (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx) const;

	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecializedDD (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									  size_t start_idx, size_t end_idx) const;

	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecializedDSP (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									   size_t start_idx, size_t end_idx) const;

	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecializedDH (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									  size_t start_idx, size_t end_idx) const;
};

// Specialization of vector domain

template <>
class VectorDomain<GF2> : private virtual VectorDomainBase<GF2>, private DotProductDomain<GF2>
{
    public:
	typedef bool Element;

	VectorDomain (const VectorDomain &VD)
		: VectorDomainBase<GF2> (VD._F), DotProductDomain<GF2> (VD._F)
	{}

	VectorDomain &operator = (const VectorDomain &) { return *this; }

	const GF2 &field () const { return _F; }
    
	template <class Vector>
	inline std::ostream &write (std::ostream &os, const Vector &x) const
		{ return writeSpecialized (os, x, typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector>
	inline std::istream &read (std::istream &is, Vector &x) const
		{ return readSpecialized (is, x, typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &copy (Vector1 &res, const Vector2 &v) const
		{ return copySpecialized (res, v,
					  typename GF2VectorTraits<Vector1>::VectorCategory (),
					  typename GF2VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector>
	inline void swap (Vector &v1, Vector &v2) const
		{ swapSpecialised (v1, v2,
				   typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector>
	inline int firstNonzeroEntry (bool &a, const Vector &v) const
		{ return firstNonzeroEntrySpecialized (a, v, typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &copy (Vector1 &res, const Vector2 &v, size_t i, size_t len = 0) const
		{ return copySpecialized (res, v, i, len,
					  typename GF2VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector, class Iterator>
	inline Vector &permute (Vector   &v,
				Iterator  P_start,
				Iterator  P_end) const
		{ return permuteSpecialized (v, P_start, P_end,
					     typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline bool areEqual (const Vector1 &v1, const Vector2 &v2) const
		{ return areEqualSpecialized (v1, v2,
					      typename GF2VectorTraits<Vector1>::VectorCategory (),
					      typename GF2VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector>
	inline bool isZero (const Vector &v) const
		{ return isZeroSpecialized (v, typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx = 0, size_t end_idx = static_cast<size_t> (-1)) const
		{ return dotSpecialized (res, v1, v2, start_idx, end_idx,
					 typename GF2VectorTraits<Vector1>::VectorCategory (),
					 typename GF2VectorTraits<Vector2>::VectorCategory ()); }

	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dot (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
							     size_t start_idx = 0, size_t end_idx = static_cast<size_t> (-1)) const
		{ return dotSpecialized (res, v1, v2, start_idx, end_idx,
					 typename GF2VectorTraits<Vector1>::VectorCategory (),
					 typename GF2VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dot (res, v1, v2); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &add (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 typename GF2VectorTraits<Vector1>::VectorCategory (),
					 typename GF2VectorTraits<Vector2>::VectorCategory (),
					 typename GF2VectorTraits<Vector3>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &addin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   typename GF2VectorTraits<Vector1>::VectorCategory (),
					   typename GF2VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &sub (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 typename GF2VectorTraits<Vector1>::VectorCategory (),
					 typename GF2VectorTraits<Vector2>::VectorCategory (),
					 typename GF2VectorTraits<Vector3>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &subin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   typename GF2VectorTraits<Vector1>::VectorCategory (),
					   typename GF2VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &neg (Vector1 &res, const Vector2 &x) const
		{ copy (res, x); return res; }

	template <class Vector>
	inline Vector &negin (Vector &y) const
		{ return y; }

	template <class Vector1, class Vector2>
	inline Vector1 &mul (Vector1 &res, const Vector2 &x, const Element a) const
		{ return mulSpecialized (res, x, a, typename GF2VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector>
	inline Vector &mulin (Vector &x, const Element a) const
		{ return mulinSpecialized (x, a, typename GF2VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &axpy (Vector1 &res, const Element a, const Vector2 &x, const Vector3 &y) const
		{ if (a) add (res, x, y); else this->copy (res, y); return res; }

	template <class Vector1, class Vector2>
	inline Vector1 &axpyin (Vector1 &y, const Element a, const Vector2 &x) const
		{ if (a) addin (y, x); return y; }

	VectorDomain (const GF2 &F)
		: VectorDomainBase<GF2> (F), DotProductDomain<GF2> (F)
	{}


	// Specialized function implementations
	template <class Vector> 
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector>
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::HybridZeroOneVectorTag) const;

	template <class Vector>
	std::istream &readSpecialized (std::istream &is, const Vector &x,
				       VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	std::istream &readSpecialized (std::istream &is, const Vector &x,
				       VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::DenseZeroOneVectorTag) const;
	
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag) const
		{ return areEqual (v2, v1); }
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::SparseZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::HybridZeroOneVectorTag,
				  VectorCategories::HybridZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::HybridZeroOneVectorTag,
				  VectorCategories::DenseZeroOneVectorTag) const
		{ return areEqual (v2, v1); }
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::HybridZeroOneVectorTag) const;
    

	template <class Vector>
	bool isZeroSpecialized (const Vector &v, VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	inline bool isZeroSpecialized (const Vector &v,
				       VectorCategories::SparseZeroOneVectorTag) const
		{ return v.empty (); }
	template <class Vector>
	inline bool isZeroSpecialized (const Vector &v,
				       VectorCategories::HybridZeroOneVectorTag) const
		{ return v.empty (); }

	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::DenseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag) const
		{ std::copy (v.wordBegin (), v.wordEnd (), res.wordBegin ()); return res; }
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::SparseZeroOneVectorTag,
				  VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::HybridZeroOneVectorTag,
				  VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::HybridZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const
		{ res = v; return res; }
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::HybridZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::SparseZeroOneVectorTag,
				  VectorCategories::HybridZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::HybridZeroOneVectorTag,
					 VectorCategories::HybridZeroOneVectorTag) const
		{ res = v; return res; }

	template <class Vector>
	inline void swapSpecialised (Vector &v1, Vector &v2,
				     VectorCategories::DenseZeroOneVectorTag) const
		{ std::swap_ranges (v1.wordBegin (), v1.wordEnd (), v2.wordBegin ()); }

	template <class Vector>
	inline void swapSpecialised (Vector &v1, Vector &v2,
				     VectorCategories::GenericVectorTag) const
		{ std::swap (v1, v2); }

	template <class Vector, class Iterator>
	inline Vector &permuteSpecialized (Vector &v, Iterator P_start, Iterator P_end,
					   VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector, class Iterator>
	inline Vector &permuteSpecialized (Vector &v, Iterator P_start, Iterator P_end,
					   VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector, class Iterator>
	inline Vector &permuteSpecialized (Vector &v, Iterator P_start, Iterator P_end,
					   VectorCategories::HybridZeroOneVectorTag) const;

	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2, start_idx, end_idx); }
	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::SparseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2, start_idx, end_idx); }
	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx,
					VectorCategories::SparseZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1, start_idx, end_idx); }
	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::HybridZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDH (res, v1, v2, start_idx, end_idx); }
	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx,
					VectorCategories::HybridZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDH (res, v2, v1, start_idx, end_idx); }
	template <class Vector1, class Vector2>
	Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2, size_t start_idx, size_t end_idx,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const;

	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecialized (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									size_t start_idx, size_t end_idx,
									VectorCategories::DenseZeroOneVectorTag,
									VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2, start_idx, end_idx); }
	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecialized (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									size_t start_idx, size_t end_idx,
									VectorCategories::DenseZeroOneVectorTag,
									VectorCategories::SparseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2, start_idx, end_idx); }
	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecialized (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									size_t start_idx, size_t end_idx,
									VectorCategories::SparseZeroOneVectorTag,
									VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1, start_idx, end_idx); }
	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecialized (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									size_t start_idx, size_t end_idx,
									VectorCategories::DenseZeroOneVectorTag,
									VectorCategories::HybridZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDH (res, v1, v2, start_idx, end_idx); }
	template <class Iterator, class Endianness, class Vector1, class Vector2>
	inline BitVectorReference<Iterator, Endianness> dotSpecialized (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
									size_t start_idx, size_t end_idx,
									VectorCategories::HybridZeroOneVectorTag,
									VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDH (res, v2, v1, start_idx, end_idx); }
	template <class Iterator, class Endianness, class Vector1, class Vector2>
	BitVectorReference<Iterator, Endianness> dotSpecialized (BitVectorReference<Iterator, Endianness> res, const Vector1 &v1, const Vector2 &v2,
								 size_t start_idx, size_t end_idx,
								 VectorCategories::SparseZeroOneVectorTag,
								 VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const
		{ copy (res, y); addin (res, x); }
	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::HybridZeroOneVectorTag,
				 VectorCategories::HybridZeroOneVectorTag,
				 VectorCategories::HybridZeroOneVectorTag) const;

	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag,
				   VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag,
				   VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag,
				   VectorCategories::HybridZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::SparseZeroOneVectorTag,
				   VectorCategories::DenseZeroOneVectorTag) const
		{ Vector1 xp, res; copy (xp, x); add (res, y, xp); copy (y, res); return y; }
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::SparseZeroOneVectorTag,
				   VectorCategories::SparseZeroOneVectorTag) const
		{ static Vector1 res; add (res, y, x); this->copy (y, res); return y; }
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::HybridZeroOneVectorTag,
				   VectorCategories::HybridZeroOneVectorTag) const
		{ static Vector1 res; add (res, y, x); this->copy (y, res); return y; }

	template <class Vector1, class Vector2>
	Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
				 VectorCategories::DenseZeroOneVectorTag tag) const
		{ if (a) this->copy (res, x); else std::fill (res.wordBegin (), res.wordEnd (), 0); return res; }
	template <class Vector1, class Vector2>
	Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
				 VectorCategories::SparseZeroOneVectorTag tag) const
		{ if (a) this->copy (res, x); else res.clear (); return res; }

	template <class Vector>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::DenseZeroOneVectorTag) const
		{ if (!a) std::fill (x.wordBegin (), x.wordEnd (), 0); return x; }

	template <class Vector>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::SparseZeroOneVectorTag tag) const
		{ if (!a) x.clear (); return x; }
	template <class Vector>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::HybridZeroOneVectorTag tag) const
		{ if (!a) { x.clear (); } return x; }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag) const
	{
		typename LinBox::Vector<GF2>::Sparse v;
		typename LinBox::Vector<GF2>::Sparse w;
		typename LinBox::Vector<GF2>::Sparse u;

		copy (v, x);
		copy (w, y);
		add (u, w, v);
		copy (res, u);

		return u;
	}

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag) const
	{
		typename LinBox::Vector<GF2>::Sparse v;
		typename LinBox::Vector<GF2>::Sparse w;
		typename LinBox::Vector<GF2>::Sparse u;

		copy (v, x);
		copy (w, y);
		sub (u, w, v);
		copy (res, u);

		return u;
	}

	template <class word, class Endianness>
	inline int firstNonzeroEntryInWord (word v) const;

	template <class Vector>
	inline int firstNonzeroEntrySpecialized (bool &a, const Vector &v,
						 VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	inline int firstNonzeroEntrySpecialized (bool &a, const Vector &v,
						 VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector>
	inline int firstNonzeroEntrySpecialized (bool &a, const Vector &v,
						 VectorCategories::HybridZeroOneVectorTag) const;
};

} // namespace LinBox

#endif // __LINBOX_VECTOR_DOMAIN_GF2_H

#include "linbox/vector/vector-domain-gf2.inl"

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
