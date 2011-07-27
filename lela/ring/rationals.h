/* lela/ring/rationals.h
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Field of rational numbers
 *
 * ----------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_RING_RATIONALS_H
#define __LELA_RING_RATIONALS_H

#include <iostream>
#include <cctype>

#include <gmp.h>

#include "lela/integer.h"
#include "lela/ring/interface.h"
#include "lela/element/rational.h"
#include "lela/lela-config.h"
#include "lela/util/debug.h"

namespace LELA
{

// Forward declarations
class RationalRandIter;

/** Field of rational numbers
 *
 * \ingroup ring
 */

class Rationals : public RingInterface<RationalElement>
{
    private:

	const integer _cardinality;
	const integer _characteristic;

    public:

	typedef RationalElement Element;

	typedef RationalRandIter RandIter;

	Rationals (const Rationals &) 
		: _cardinality (0), _characteristic (0), _zero (0, 1), _one (1, 1), _minus_one (-1, 1)
		{}

	~Rationals () {}

	Rationals &operator = (const Rationals &F)
	{
		_zero = F._zero;
		_one = F._one;
		_minus_one = F._minus_one;
		return *this;
	}
    
	Element &init (Element &x, const integer &y = 0) const
	{
		mpq_set_z (x.rep, y.get_mpz_t ());
		return x;
	}

	Element &init (Element &x, const integer &num, const integer &den) const
	{
		init (x, num);
		Element y;
	        init (y, den);
	        divin (x, y);
	        return x;
	}
  
	Element &copy (Element &x, const Element &y) const
	{
		mpq_set (x.rep, y.rep);
		return x;
	}

	integer &cardinality (integer &c) const 
		{ return c = _cardinality; }

	integer &characteristic (integer &c) const
		{ return c = _characteristic; }

	bool areEqual (const Element &x, const Element &y) const
		{ return mpq_equal (x.rep, y.rep); }

	Element &add (Element &x, const Element &y, const Element &z) const
	{
		mpq_add (x.rep, y.rep, z.rep);
		return x;
	}
    
	Element &sub (Element &x, const Element &y, const Element &z) const
	{
		mpq_sub (x.rep, y.rep, z.rep);
		return x;
	}
    
	Element &mul (Element &x, const Element &y, const Element &z) const
	{
		mpq_mul (x.rep, y.rep, z.rep);
		return x;
	}
    
	bool div (Element &x, const Element &y, const Element &z) const
	{
		if (!isZero (z)) {
			mpq_div (x.rep, y.rep, z.rep);
			return true;
		} else
			return false;
	}

	Element &neg (Element &x, const Element &y) const
	{
		mpq_neg (x.rep, y.rep);
		return x;
	}

	bool inv (Element &x, const Element &y) const
	{
		if (!isZero (y)) {
			mpq_inv (x.rep, y.rep);
			return true;
		} else
			return false;
	}

	bool isZero (const Element &x) const 
		{ return mpq_sgn (x.rep) == 0; }
    
	bool isOne (const Element &x) const 
		{ return mpq_cmp_ui (x.rep, 1L, 1L) == 0; }
    
	Element &addin (Element &x, const Element &y) const
	{
		mpq_add (x.rep, x.rep, y.rep);
		return x;
	}

	Element &subin (Element &x, const Element &y) const
	{
		mpq_sub (x.rep, x.rep, y.rep);
		return x;
	}
 
	Element &mulin (Element &x, const Element &y) const
	{
		mpq_mul (x.rep, x.rep, y.rep);
		return x;
	}

	Element &axpy (Element &r, const Element &a, const Element &x, const Element &y) const
	{
		mpq_mul (r.rep, a.rep, x.rep);
		mpq_add (r.rep, r.rep, y.rep);
		return r;
	}

	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		Element tmp;
		mpq_mul (tmp.rep, a.rep, x.rep);
		mpq_add (r.rep, r.rep, tmp.rep);
		return r;
	}

	bool divin (Element &x, const Element &y) const
	{
		if (!isZero (y)) {
			mpq_div (x.rep, x.rep, y.rep);
			return true;
		} else
			return false;
	}

	Element &negin (Element &x) const
	{
		mpq_neg (x.rep, x.rep);
		return x;
	}

	bool invin (Element &x) const
	{
		if (!isZero (x)) {
			mpq_inv (x.rep, x.rep);
			return true;
		} else
			return false;
	}
    
	std::ostream &write (std::ostream &os) const 
		{ return os << "QQ"; }
    
	std::istream &read (std::istream &is) { return is; }

	std::ostream &write (std::ostream &os, const Element &x) const;

	std::istream &read (std::istream &is, Element &x) const;

	size_t elementWidth () const
		{ return 10; }

	Rationals (int p = 0, int exp = 1)
		: _cardinality (0), _characteristic (0), _zero (0, 1), _one (1, 1), _minus_one (-1, 1)
	{
		if (p != 0)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "modulus must be 0 (no modulus)");

		if (exp != 1)
			throw PreconditionFailed (__FUNCTION__, __LINE__, "exponent must be 1");
	}
    
	integer &get_num (integer &x, const Element &y) const
	{
		mpq_get_num (x.get_mpz_t (), y.rep);
		return x;

	}

	integer &get_den (integer &x, const Element &y) const
	{
		mpq_get_den (x.get_mpz_t (), y.rep);
		return x;
	}

	int sign (const Element& x) const
		{ return mpq_sgn (x.rep); }
	
        integer &bitsize (integer &bs, const Element q) const
	{
                integer y; get_den (y, q);
                integer x; get_num (x, q);
                bs = mpz_sizeinbase (x.get_mpz_t (), 2) + mpz_sizeinbase (y.get_mpz_t (), 2);
                return bs;
        }

	const Element &zero () const { return _zero; }
	const Element &one () const { return _one; }
	const Element &minusOne () const { return _minus_one; }

private:

	Element _zero, _one, _minus_one;

}; // class Rationals

} // namespace LELA

#include "lela/randiter/rationals.h"

#endif // __LELA_RING_RATIONALS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
