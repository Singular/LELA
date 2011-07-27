/* lela/util/error.h
 * Copyright 1994-1997 T. Gautier
 *
 * Written by T. Gautier
 *
 * Support-routines for errors
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_UTIL_ERROR_H
#define __LELA_UTIL_ERROR_H

#include <cstring>
#include <iostream>

namespace LELA
{

/** Base class for exceptions
 *
 * \ingroup util
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

/// Exception class for functions which haven't been implemented
///
/// \ingroup util
class NotImplemented : public LELAError
{
public:
	NotImplemented () : LELAError ("Sorry, the requested function is not yet implemented.") {}
};

class DiagonalEntryNotInvertible 
{
	friend std::ostream &operator << (std::ostream &os, const DiagonalEntryNotInvertible &e)
		{ os << "Entry on the diagonal is not invertible" << std::endl; return os; }
};
 
}

#endif // __LELA_UTIL_ERROR_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
