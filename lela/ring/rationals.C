/* lela/ring/rationals.C
 * Copyright 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Support-routines for rational numbers
 *
 * ----------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/ring/rationals.h"

namespace LELA
{

std::ostream &Rationals::write (std::ostream &os, const Element &x) const 
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

std::istream &Rationals::read (std::istream &is, Element &x) const
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

}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
