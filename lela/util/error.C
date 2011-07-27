/* lela/util/error.C
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

#include "lela/lela-config.h"

#include "lela/util/error.h"

namespace LELA
{

std::ostream& operator<< (std::ostream& o, const LELAError& E) 
{
	E.print(o) ; 
	return o ;
}

} // namespace LELA

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax

