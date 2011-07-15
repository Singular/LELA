/* lela/util/error.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
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

#ifndef __LELA_util_error_H
#define __LELA_util_error_H

#include <cstring>
#include <iostream>

namespace LELA
{

// ------------------------------- LELAError
/** base class for execption handling in Givaro
\ingroup util
*/
class LELAError {
	static const size_t max_error_string = 256;
public:
	LELAError (const char* msg = '\0') {
		std::strncpy(strg, msg, max_error_string);
		strg[max_error_string-1] = 0;
	};


	// -- virtual print of the error message
	virtual std::ostream &print (std::ostream &o) const
	{ return o << strg<<std::endl ; }
  
	// -- non virtual output operator
	friend std::ostream &operator << (std::ostream &o, const LELAError &E);

	// - useful to setup a break point on it
	static void throw_error (const LELAError &err)
		{ throw err; }

    	virtual ~LELAError() {}        

    protected:
	char strg[max_error_string]; 
};

class LELAMathError : public LELAError {
 public:
	LELAMathError (const char* msg) : LELAError (msg) {};
};

class LELAMathDivZero : public LELAMathError {
 public:
	LELAMathDivZero (const char* msg) : LELAMathError (msg) {};
};

class LELAMathInconsistentSystem : public LELAMathError {
 public:
	LELAMathInconsistentSystem (const char* msg) : LELAMathError (msg) {};
};

// -- Exception thrown in input of data structure 
class LELABadFormat : public LELAError {
 public:
	LELABadFormat (const char* msg) : LELAError (msg) {};
};

class DiagonalEntryNotInvertible 
{
	friend std::ostream &operator << (std::ostream &os, const DiagonalEntryNotInvertible &e)
		{ os << "Entry on the diagonal is not invertible" << std::endl; return os; }
};
 
}

#endif // __LELA_util_error_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax