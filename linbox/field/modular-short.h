/* Copyright (C) 2010 LinBox
 * Written by <?>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_modular__int16_H
#define __LINBOX_modular__int16_H

#include <cmath>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/debug.h"

#ifndef LINBOX_MAX_INT16
#define LINBOX_MAX_INT16 32767
#endif

// This is replaced by FieldTraits< Modular<int16> >::maxModulus(integer&)
// #ifndef LINBOX_MAX_INT16_MODULUS
// #define LINBOX_MAX_INT16_MODULUS 32767
// #endif

// Namespace in which all LinBox code resides
namespace LinBox 
{

template <class Element>
class Modular;

template <class Element>
class ModularRandIter;

template <class Field>
class DotProductDomain;

template <class Field>
class FieldAXPY;

template <class Field>
class MVProductDomain;

template <class Ring>
struct ClassifyRing; 
	
template <class Element>
struct ClassifyRing<Modular<Element> >;

template <>
struct ClassifyRing<Modular<short> >
{
	typedef RingCategories::ModularTag categoryTag;
};

template <class Element>
struct ModularFieldTraits;

template <>
struct ModularFieldTraits<int16>
{
	typedef int16 Element;
	typedef uint16 UnsignedElement;
	typedef uint32 DoubleSizedElement;

	static const int16 defaultModulus = 251;
	static const int16 maxModulus = 32767;  // 2^15 - 1
	static const int16 maxInt = 32767;
};

/** \brief Specialization of Modular to short element type with efficient dot product.
 * 
 * Efficient element operations for dot product, mul, axpy, by using floating point
 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
 * 
 * Requires: modulus < 2^15. 
 * Intended use: 2^7 < prime modulus < 2^15.
 \ingroup field
*/
template <>
class Modular<int16> : public FieldInterface
{
public:	       

	typedef int16 Element;

protected:

	Element modulus;
	double modulusinv;

	Element _two64;

	void init_two64 (Element modulus)
	{
		_two64 = (Element) ((uint64) (-1) % (uint64) modulus);
		_two64 += 1;

		if (_two64 >= modulus)
			_two64 -= modulus;
	}

public:	       

	friend class FieldAXPY<Modular<Element> >;
	friend class DotProductDomain<Modular<Element> >;
	friend class MVProductDomain<Modular<Element> >;

	typedef ModularRandIter<Element> RandIter;

	Modular ()
		: modulus (ModularFieldTraits<Element>::defaultModulus)
	{
		modulusinv = 1 / (double) modulus;
		init_two64 (modulus);
	}

	Modular (Element value, int exp = 1)
		: modulus (value)
	{
		modulusinv = 1 / ((double) value); 

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");
		if (value <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (value > FieldTraits< Modular<Element> >::maxModulus (max))
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");

		init_two64 (value);
	}

	Modular (integer value, int exp = 1)
		: modulus (value.get_si ())
	{
		modulusinv = 1 / value.get_d ();

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");
		if (value <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (value > FieldTraits< Modular<Element> >::maxModulus (max))
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");

		init_two64 (value.get_si ());
	}

	Modular (const Modular<Element> &mf)
		: modulus (mf.modulus), modulusinv (mf.modulusinv), _two64 (mf._two64)
		{}

	const Modular &operator = (const Modular<Element> &F)
	{
		modulus = F.modulus;
		modulusinv = F.modulusinv;
		_two64 = F._two64;
		return *this;
	}

	inline integer &cardinality (integer &c) const
		{ return c = modulus; }

	inline integer &characteristic (integer &c) const
		{ return c = modulus; }

	inline integer &convert (integer &x, const Element &y) const
		{ return x = y; }

	inline Element &convert (Element &x, const Element &y) const
		{ return x = y; }

	inline double &convert (double &x, const Element &y) const
		{ return x = (double) y; }

	inline float &convert (float &x, const Element &y) const
		{ return x = (float) y; }
		
	inline std::ostream &write (std::ostream &os) const
		{ return os << "int16 mod " << modulus; }
		
	inline std::istream &read (std::istream &is)
	{
		int prime;

		is >> prime;
		modulus = prime;
		modulusinv = 1 / ((double) modulus);

		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (prime > FieldTraits< Modular<Element> >::maxModulus (max))
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");

		init_two64 (modulus);

		return is;
	}
		
	inline std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }

	inline std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;
		init (x, tmp);

		return is;
	}
		
	inline Element &init (Element &x, const integer &y) const
	{
		x = y.get_ui () % modulus;
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &init (Element &x, int y = 0) const
	{
		x = y % modulus;
		if ( x < 0 ) x += modulus;
		return x;
	}

	inline Element &init (Element &x, unsigned int y = 0) const
	{
		x = y % modulus;
		if ( x < 0 ) x += modulus;
		return x;
	}

	inline Element &init (Element &x, long y) const
	{
		x = y % modulus;
		if ( x < 0 ) x += modulus;
		return x;
	}

	inline Element &init (Element &x, const float &y) const
		{ return init (x, (double) y); }

	inline Element &init (Element &x, const double &y) const
	{
		double z = fmod (y, (double) modulus);

		if (z < 0)
			z += (double) modulus;

		//z += 0.5; // C Pernet Sounds nasty and not necessary

		return x = static_cast<long> (z); //rounds towards 0
	}

	inline Element &assign (Element &x, const Element &y) const
		{ return x = y; }

	inline bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	inline  bool isZero (const Element &x) const
		{ return x == 0; }
		
	inline bool isOne (const Element &x) const
		{ return x == 1; }

	inline Element &add (Element &x, const Element &y, const Element &z) const
	{
		x = y + z;

		if ((ModularFieldTraits<Element>::UnsignedElement) x >= (ModularFieldTraits<Element>::UnsignedElement) modulus)
			x = ((ModularFieldTraits<Element>::UnsignedElement) x) - modulus;

		return x;
	}
 
	inline Element &sub (Element &x, const Element &y, const Element &z) const
	{
		x = y - z;
		if (x < 0) x += modulus;
		return x;
	}
		
	inline Element &mul (Element &x, const Element &y, const Element &z) const
	{
		Element q;

		double ab = ((double) y) * ((double) z);
		q = (Element) (ab * modulusinv);  // q could be off by (+/-) 1
		x = (Element) (ab - ((double) q) * ((double) modulus));

		if (x >= modulus)
			x -= modulus;
		else if (x < 0)
			x += modulus;

		return x;
	}
 
	inline Element &div (Element &x, const Element &y, const Element &z) const
	{
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	inline Element &neg (Element &x, const Element &y) const
	{
		if (y == 0)
			return x = 0;
		else
			return x = modulus - y;
	}
 
	inline Element &inv (Element &x, const Element &y) const
	{
		Element d, t;

		XGCD (d, x, t, y, modulus);

		if (d != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "InvMod: Input is not invertible");

		if (x < 0)
			return x += modulus;
		else
			return x;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const
	{
		Element q;

		double ab = ((double) a) * ((double) x) + (double) y;		
		q = (Element) (ab * modulusinv);  // q could be off by (+/-) 1
		r = (Element) (ab - ((double) q) * ((double) modulus));

		if (r >= modulus)
			r -= modulus;
		else if (r < 0)
			r += modulus;

		return r;
	}

	inline Element &addin (Element &x, const Element &y) const
	{
		x += y;

		if (((ModularFieldTraits<Element>::UnsignedElement) x) >= (ModularFieldTraits<Element>::UnsignedElement) modulus)
			x = ((ModularFieldTraits<Element>::UnsignedElement) x) - modulus;

		return x;
	}
 
	inline Element &subin (Element &x, const Element &y) const
	{
		x -= y;
		if (x < 0) x += modulus;
		return x;
	}
 
	inline Element &mulin (Element &x, const Element &y) const
		{ return mul (x, x, y); }
 
	inline Element &divin (Element &x, const Element &y) const
		{ return div (x, x, y); }
 
	inline Element &negin (Element &x) const
	{
		if (x == 0)
			return x; 
		else
			return x = modulus - x; 
	}
 
	inline Element &invin (Element &x) const
		{ return inv (x, x); }

	inline Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		Element q;

		double ab = ((double) a) * ((double) x) + (double) r;
		q = (Element) (ab * modulusinv);  // q could be off by (+/-) 1
		r = (Element) (ab - ((double) q) * ((double) modulus));

		if (r >= modulus)
			r -= modulus;
		else if (r < 0)
			r += modulus;

		return r;
	}

	static inline Element getMaxModulus()
		{ return ModularFieldTraits<Element>::maxModulus; }

	Element zero () const { return 0; }
	Element one () const { return 1; }
	Element minusOne () const { return modulus - 1; }

private:

	static void XGCD (Element &d, Element &s, Element &t, Element a, Element b)
	{
		Element u, v, u0, v0, u1, v1, u2, v2, q, r;
			
		Element aneg = 0, bneg = 0;
			
		if (a < 0) {
			if (a < -ModularFieldTraits<Element>::maxInt)
				throw PreconditionFailed (__FUNCTION__, __LINE__, "XGCD: integer overflow");

			a = -a;
			aneg = 1;
		}

		if (b < 0) {
			if (b < -ModularFieldTraits<Element>::maxInt)
				throw PreconditionFailed (__FUNCTION__, __LINE__, "XGCD: integer overflow");

			b = -b;
			bneg = 1;
		}

		u1 = 1; v1 = 0;
		u2 = 0; v2 = 1;
		u = a; v = b;
			
		while (v != 0) {
			q = u / v;
			r = u % v;
			u = v;
			v = r;
			u0 = u2;
			v0 = v2;
			u2 = u1 - q*u2;
			v2 = v1 - q*v2;
			u1 = u0;
			v1 = v0;
		}
			
		if (aneg)
			u1 = -u1;
			
		if (bneg)
			v1 = -v1;
			
		d = u;
		s = u1;
		t = v1;
	}
};

template <>
class FieldAXPY<Modular<int16> >
{	  
public:
	  
	typedef int16 Element;
	typedef Modular<int16> Field;
	  
	FieldAXPY (const Field &F)
		: _F (F), _y (0)
		{}

	FieldAXPY (const FieldAXPY &faxpy)
		: _F (faxpy._F), _y (0)
		{}
	  
	FieldAXPY<Modular<Element> > &operator = (const FieldAXPY &faxpy)
	{
		_F = faxpy._F;
		_y = faxpy._y;
		return *this; 
	}
	  
	inline uint64 &mulacc (const Element &a, const Element &x)
	{
		uint64 t = ((uint32) a) * ((uint32) x);
		return _y += t;
	}

	inline uint64 &accumulate (const Element &t)
		{ return _y += t; }

	inline Element &get (Element &y)
	{
		y = _y % (uint64) _F.modulus;
		return y;
	}

	inline FieldAXPY &assign (const Element y)
	{
		_y = y; 
		return *this;
	}

	inline void reset()
		{ _y = 0; }

protected:

	Field _F;
	uint64 _y;
};

template <>
class DotProductDomain<Modular<int16> > : private virtual VectorDomainBase<Modular<int16> >
{

public:	  
	typedef int16 Element;

	DotProductDomain (const Modular<Element> &F)
		: VectorDomainBase<Modular<Element> > (F) {}
	  
protected:
	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		  
		uint64 y = 0;
		// uint64 t;
		  
		for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
			y += ((uint32) *i) * ((uint32) *j);
		}

		y %= (uint64) _F.modulus;

		return res = y;
	}
	  
	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDS (Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		typename Vector1::const_iterator i;
		  
		uint64 y = 0;
		  
		for (i = v1.begin (); i != v1.end (); ++i) {
			y += ((uint32) i->second) * ((uint32) v2[i->first]);
		}

		y %= (uint64) _F.modulus;

		return res = y;
	}
};

// Specialization of MVProductDomain for int16 modular field	

template <>
class MVProductDomain<Modular<int16> >
{
public:

	typedef int16 Element;

protected:
	template <class Vector1, class Matrix, class Vector2>
	inline Vector1 &mulColDense
		(const VectorDomain<Modular<Element> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
	{
		return mulColDenseSpecialized
			(VD, w, A, v, typename DefaultVectorTraits<typename Matrix::Column>::VectorCategory ());
	}

private:
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<Element> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<Element> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseVectorTag) const;

	mutable std::vector<uint64> _tmp;
};

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<int16> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int16> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());
		
	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());
		
	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);
		
	for (j = v.begin (); j != v.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((ModularFieldTraits<Element>::DoubleSizedElement) *k) * ((ModularFieldTraits<Element>::DoubleSizedElement) *j);

			*l += t;
				
		}
	}
		
	typename Vector1::iterator w_j;
		
	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l % VD.field ().modulus;
		
	return w;
}
	
template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<int16> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int16> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());
			
	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;
			
	uint64 t;
			
	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());
			
	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);
			
	for (j = v.begin (); j != v.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((ModularFieldTraits<Element>::DoubleSizedElement) k->second) * ((ModularFieldTraits<Element>::DoubleSizedElement) *j);

			_tmp[k->first] += t;
					
		}
	}
			
	typename Vector1::iterator w_j;
			
	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l % VD.field ().modulus;
			
	return w;
}

} // namespace LinBox

#include "linbox/randiter/modular.h"

#endif //__LINBOX_modular__int16_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
