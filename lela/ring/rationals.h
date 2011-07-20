/* lela/ring/rationals.h
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ----------------------------------------
 *
 * See COPYING for license information
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

// Namespace in which all LELA library code resides
namespace LELA
{

// Forward declarations
class RationalRandIter;

/**
 * \brief Field of rational numbers using GMP
 \ingroup ring
 *
 * This is a wrapper for the GMP rational number facility, built to the
 * interface of the ring archetype. 
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
  
	integer &convert (integer &x, const Element &y) const
	{
		mpz_t n, d;

		mpz_init (n);
		mpz_init (d);
		mpq_get_num (n, y.rep);
		mpq_get_den (d, y.rep);

		mpz_divexact (x.get_mpz_t (), n, d);

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
	{ 
		os << "QQ"; 
		return os;
	}
    
	std::istream &read (std::istream &is) { return is; }

	std::ostream &write (std::ostream &os, const Element &x) const 
	{
		char *str;

		str = new char[mpz_sizeinbase (mpq_numref (x.rep), 10) + 2];
		mpz_get_str (str, 10, mpq_numref (x.rep));
		os << str;
		delete[] str;

		if (mpz_cmp_ui (mpq_denref (x.rep), 1L) != 0) {
			str = new char[mpz_sizeinbase (mpq_denref (x.rep), 10) + 2];
			mpz_get_str (str, 10, mpq_denref (x.rep));
			os << '/' << str;
			delete[] str;
		}

		return os;
	}

	std::istream &read (std::istream &is, Element &x) const
	{
		char buffer[65535], endc;
		bool found_space = false;
		int i = 0;			

		do {
			is.get (endc);
		} while (is && !isdigit (endc) && endc != '-' && endc != '.' &&  endc !='e' && endc != 'E');		

		buffer[i]=endc;
		
		while ((buffer[i] == '-' || isdigit (buffer[i])) && i < 65535)  {
			i++;
			is.get (buffer[i]);
		}

		endc = buffer[i];
	       	buffer[i] = '\0';

		if (i > 0)
			mpz_set_str (mpq_numref (x.rep), buffer, 10);
		else
			mpq_set_si (x.rep, 0L, 1L);

		if (endc == ' ') {
			found_space = true;
			while (endc == ' ') is >> endc;
		}

		if (endc == '/') {
			i = 0;

			is.get (endc);
			while (isspace (endc)) is.get (endc);
			is.putback (endc);

			do {
				is.get (buffer[i++]);
			} while (isdigit (buffer[i - 1]) && i < 65536);

			is.putback (buffer[i - 1]);
			buffer[i - 1] = '\0';

			mpz_set_str (mpq_denref (x.rep), buffer, 10);
			mpq_canonicalize (x.rep);
			return is;
		} else {
			 mpz_set_si (mpq_denref (x.rep), 1L);
		}

		if (endc == '.' && !found_space) {
			Element decimal_part;

			mpz_set_si (mpq_denref (x.rep), 1L);
			mpq_set_si (decimal_part.rep, 1L, 1L);
			mpz_set_si (mpq_denref (decimal_part.rep), 1L);

			i = 0;

			do {
				is.get (buffer[i++]);
				if (isdigit (buffer[i - 1]))
					mpz_mul_ui (mpq_denref (decimal_part.rep),
						    mpq_denref (decimal_part.rep), 10L);
			} while (isdigit (buffer[i - 1]) && i < 65536);

			is.putback (buffer[i - 1]);
			buffer[i - 1] = '\0';

			mpz_set_str (mpq_numref (decimal_part.rep), buffer, 10);
			mpq_canonicalize (decimal_part.rep);

			mpq_add (x.rep, x.rep, decimal_part.rep);

			do {
				is.get (endc);
			} while (is && endc == ' ') ;
		}
		
		if ((endc == 'e') || (endc == 'E')) {
	                is.get(endc);

	                bool minus = false;

                        if (endc == '-')
                                minus = true;
                        else if (endc == '+')
                                minus = false;
                        else
	                        is.putback(endc);

	                i=0;

	                do {
		                is.get (buffer[i++]);
		        } while (isdigit (buffer[i-1]) && i < 65536);

		        is.putback (buffer[i-1]);
		        buffer[i-1] = '\0';

			integer pow (buffer), powten = 1;

			for (integer it = 0; it < pow; ++it)
				powten *= 10;

                        if (minus)
	                        div (x, x, Element (powten, 1));
	                else
	                        mul (x, x, Element (powten, 1));
                }
		else {
			is.putback (endc);
		//	mpz_set_si (mpq_denref (x.rep), 1L);
		}

		mpq_canonicalize (x.rep);

		return is;
	}

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
