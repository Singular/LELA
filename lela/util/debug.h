/* lela/util/debug.h
 * Copyright 2001,2010 LELA
 * Copyright 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Some routines for checking preconditions
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#ifndef __LELA_UTIL_DEBUG_H
#define __LELA_UTIL_DEBUG_H

#include <iostream>
#include <sstream>
#include "lela/util/error.h"

#ifndef DEBUG
#  define lela_check(check)
#else
#  ifdef __GNUC__
#    define lela_check(check) \
        if (!(check)) \
                 throw LELA::PreconditionFailed (__func__, __FILE__, __LINE__, #check); //BB : should work on non gnu compilers too
#  else
#    define lela_check(check) \
        if (!(check)) \
                 throw LELA::PreconditionFailed (__FILE__, __LINE__, #check);
#  endif
#endif

namespace LELA
{

/** Exception thrown when a precondition-check fails
 *
 * \ingroup util
 */
class PreconditionFailed //: public LELAError BB: otherwise,  error.h:39 segfaults
{
	static std::ostream *_errorStream;

public:
	PreconditionFailed (const char *function, int line, const char *check) { //BB : not so useful, is it ?
		if (_errorStream == (std::ostream *) 0)
			_errorStream = &std::cerr;

		(*_errorStream) << std::endl << std::endl;
		(*_errorStream) << "ERROR (" << function << ":" << line << "): ";
		(*_errorStream) << "Precondition not met:" << check << std::endl;
	}

	PreconditionFailed (const char* function, const char *file, int line, const char *check) {
		if (_errorStream == (std::ostream *) 0)
			_errorStream = &std::cerr;

		(*_errorStream) << std::endl << std::endl;
		(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
		(*_errorStream) << "Precondition not met:" << check << std::endl;
	}

	static void setErrorStream (std::ostream &stream);

	// -- overload the virtual print of LELAError
	std::ostream &print (std::ostream &o) const { 
		if (std::ostringstream * str = dynamic_cast<std::ostringstream*>(_errorStream))
			return o << str->str() ; 
		else
			throw LELAError("LELA ERROR: PreconditionFailed exception is not initialized correctly");
	}
};

} // namespace LELA

#endif // __LELA_UTIL_DEBUG_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
