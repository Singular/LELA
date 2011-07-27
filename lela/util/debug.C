/* lela/util/debug.C
 * Copyright 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include "lela/lela-config.h"

#include <iostream>

#include "lela/util/debug.h"

namespace LELA
{

void PreconditionFailed::setErrorStream (std::ostream &stream)
	{ _errorStream = &stream; }

std::ostream *PreconditionFailed::_errorStream;

} // namespace LELA

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
