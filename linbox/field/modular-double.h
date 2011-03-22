/* linbox/field/modular-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_modular_double_H
#define __LINBOX_modular_double_H

#include <math.h>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/randiter/modular.h"

// Namespace in which all LinBox code resides
namespace LinBox
{ 

template <class Element>
class Modular;

template <class Element>
class ModularRandIter;

template <class Field, class RandIter>
class NonzeroRandIter;

template <class Ring>
struct ClassifyRing; 

template <class Element>
struct ClassifyRing<Modular<Element> >;

template <>
struct ClassifyRing<Modular<double> >
{
	typedef RingCategories::ModularTag categoryTag;
};

class MultiModDouble;
	
/// \ingroup field
template <>
class Modular<double> : public FieldInterface
{

protected:

	double modulus;
	unsigned long lmodulus;

public:	       
	friend class FieldAXPY<Modular<double> >;
	friend class DotProductDomain<Modular<double> >;
	friend class MultiModDouble;
			       
	typedef double Element;

public:	       
	typedef ModularRandIter<double> RandIter;
	typedef NonzeroRandIter<Modular<double>, ModularRandIter<double> > NonZeroRandIter;

	static ClassifyRing<Modular<double> >::categoryTag getCategory ()
		{ return ClassifyRing<Modular<double> >::categoryTag (); }

	Modular () {}

	Modular (int32 p, int exp = 1)
		: modulus ((double) p),
		  lmodulus (p)
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");

		integer max;

		if (modulus > FieldTraits<Modular<double> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (double p)
		: modulus (p),
		  lmodulus ((unsigned long) p)
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if (modulus > FieldTraits<Modular<double> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (long int p)
		: modulus ((double) p),
		  lmodulus (p)
	{
		if ((double) modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		integer max;

		if ((double) modulus > FieldTraits<Modular<double> >::maxModulus (max).get_d ())
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (const integer &p)
		: modulus (p.get_d ()),
		  lmodulus (p.get_d ())
	{
		if (modulus <= 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");
		if (modulus > 94906265)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");
	}

	Modular (const Modular<double> &mf)
		: modulus (mf.modulus),
		  lmodulus (mf.lmodulus)
		{}

	const Modular &operator = (const Modular<double> &F)
	{
		modulus = F.modulus;
		lmodulus = F.lmodulus;
		return *this;
	}

	integer &cardinality (integer &c) const
		{ return c = integer (modulus); }

	integer &characteristic (integer &c) const
		{ return c = integer (modulus); }

	integer &convert (integer &x, const Element &y) const
		{ return x = integer (y); }

	double &convert (double &x, const Element &y) const
		{ return x=y; }
		
	std::ostream &write (std::ostream &os) const
		{ return os << "double mod " << (int) modulus; }

	std::istream &read (std::istream &is)
	{
		is >> modulus;

		if (modulus <= 1) 
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be > 1");

		if (modulus > 94906265) 
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus is too big");

		return is;
	}

	std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << (int) x; }

	std::istream &read (std::istream &is, Element &x) const
	{
		integer tmp;

		is >> tmp;
		init (x, tmp);

		return is;
	}

	Element &init (Element &x, const integer &y) const
		{ return x = (Element) (y.get_ui () % lmodulus); }

	inline Element &init (Element &x, double y = 0) const
	{
		x = fmod (y, modulus);
		if (x < 0) x += modulus;
		return x;
	}
		
	inline Element &assign(Element &x, const Element &y) const
		{ return x = y; }

	inline bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

	inline bool isZero (const Element &x) const
		{ return x == 0.; }

	inline bool isOne (const Element &x) const
		{ return x == 1.; }

	inline bool isMinusOne (const Element &x) const
		{ return (x == modulus - 1.); }

	inline Element &add (Element &x, const Element &y, const Element &z) const
	{
		x = y + z;
		if (x >= modulus) x -= modulus;
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
		double tmp = y * z;
		x = fmod (tmp, modulus);
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
		// The extended Euclidean algoritm 
		int x_int, y_int, q, tx, ty, temp;
		x_int = int (modulus);
		y_int = int (y);
		tx = 0; 
		ty = 1;
		  
		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			q = x_int / y_int; // integer quotient
			temp = y_int; y_int = x_int - q * y_int;
			x_int = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}
		  
		if (tx < 0)
			tx += (int) modulus;
		  
		// now x_int = gcd (modulus,residue)
		return x = (double) tx;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const
	{
		double tmp = a * x + y;
		return r = fmod (tmp, modulus); 
	}

	inline Element &addin (Element &x, const Element &y) const
	{
		x += y;
		if (x >= modulus) x -= modulus;
		return x;
	}
 
	inline Element &subin (Element &x, const Element &y) const
	{
		x -= y;
		if (x < 0.) x += modulus;
		return x;
	}
 
	inline Element &mulin (Element &x, const Element &y) const
		{ return mul (x, x, y); }
 
	inline Element &divin (Element &x, const Element &y) const
		{ return div (x, x, y); }
 
	inline Element &negin (Element &x) const
	{
		if (x == 0.)
			return x; 
		else 
			return x = modulus - x; 
	}
		
	inline Element &invin (Element &x) const
	{
		return inv (x, x);
	}
		
	inline Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		double tmp = r + a * x;
		return r = fmod (tmp, modulus); 
	}

	static inline double getMaxModulus()
		{ return 67108864.0; } // 2^26 

	double zero () const { return 0.0; }
	double one () const { return 1.0; }
	double minusOne () const { return modulus - 1.0; }
};

template <>
class FieldAXPY<Modular<double> >
{
public:
	  
	typedef double Element;
	typedef Modular<double> Field;
	  
	FieldAXPY (const Field &F)
		: _F (F),
		  _y(0.),
		  _bound ((double) ((1ULL << 53) - (int) (_F.modulus * _F.modulus)))
		{}
	  
	FieldAXPY (const FieldAXPY &faxpy)
		: _F (faxpy._F),
		  _y (faxpy._y),
		  _bound (faxpy._bound)
		{}
	  
	FieldAXPY<Modular<double> > &operator = (const FieldAXPY &faxpy)
	{
		_F = faxpy._F; 
		_y = faxpy._y;
		_bound = faxpy._bound;
		return *this; 
	}
	  
	inline Element &mulacc (const Element &a, const Element &x)
		{ return accumulate (a * x); }

	inline Element &accumulate (const Element &tmp)
	{
                _y += tmp;

                if (_y > _bound)
			return _y = fmod (_y, _F.modulus);
                else
			return _y;
	}

	inline Element &subumulate (const Element &tmp)
	{
                _y -= tmp;

                if (_y < 0)
			return _y += _F.modulus;
                else
			return _y;
	}

	inline Element &get (Element &y)
	{
		_y = fmod (_y, _F.modulus);
		return y=_y ;
	}

	inline FieldAXPY &assign (const Element y)
	{
		_y = y; 
		return *this;
	}

	inline void reset ()
		{ _y = 0.; }
            
	inline Element& set (const Element &tmp)
	{
                _y = tmp;

                if (_y > _bound)
			return _y = fmod (_y, _F.modulus);
                else
			return _y;
	}
	  
private:
	  
	Field _F;
	double _y;
	double _bound;		
};
	
template <>
class DotProductDomain<Modular<double> > : private virtual VectorDomainBase<Modular<double> >
{
private:
	double _bound;
	size_t _nmax;
	  
public:	  
	typedef double Element;	  
	DotProductDomain (const Modular<double> &F)
		: VectorDomainBase<Modular<double> > (F),
		  _bound ((double) ((1ULL << 53) - (int) (F.modulus * F.modulus)))
	{
		_nmax = (size_t) floor ((double (1 << 26) * double (1 << 26) * 2.) / (F.modulus * F.modulus));
		_nmax = (_nmax > 0 ? _nmax : 1);
	}
	  
protected:
	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		double y = 0.;
		double t = 0.;

		if (v1.size() < _nmax) {
			for (size_t i = 0; i < v1.size(); ++i)
				y += v1[i] * v2[i];

			y = fmod(y, _F.modulus);
		} else {			
			size_t i = 0;

			for (; i < v1.size() - _nmax; i = i + _nmax) {
				for (size_t j = i; j < i + _nmax; ++j)
					y += v1[j] * v2[j];

				t += fmod (y, _F.modulus);
				y = 0.;
			}

			for (; i < v1.size(); ++i)
				y += v1[i] * v2[i];

			t += fmod (y, _F.modulus);
			y = fmod (t, _F.modulus);
		}

		return res = y;
	}

	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		double y = 0.;
		double t = 0.;

		if (v1.first.size() < _nmax) {
			for (size_t i = 0; i < v1.first.size(); ++i)
				y += v1.second[i] * v2[v1.first[i]];

			y = fmod (y, _F.modulus);
		}
		else {
			size_t i = 0;

			for (; i < v1.first.size() - _nmax; i = i + _nmax) {
				for (size_t j = i; j < i + _nmax; ++j)
					y += v1.second[j] * v2[v1.first[j]];

				t += fmod (y, this->_F.modulus);
				y = 0.;
			}

			for (;i < v1.first.size(); ++i)
				y += v1.second[i] * v2[v1.first[i]];

			t += fmod (y, _F.modulus);
			y = fmod (t, _F.modulus);
		}

		return res = y;
	}
};

} // namespace LinBox

#endif //__LINBOX_modular_double_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
