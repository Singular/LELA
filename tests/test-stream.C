
/* tests/test-stream.C
 * Copyright 2011 Florian Fischer
 *
 * Written by Florian Fischer <florian-fischer@gmx.net>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Generic test suite for splicer routines
 */

#include <iostream>
#include "lela/lela-config.h"
#include "lela/util/commentator.h"
#include <lela/blas/context.h>
#include <lela/vector/stream.h>
#include <lela/blas/level1.h>
#include <lela/vector/traits.h>
#include "lela/blas/context.h"
#include "lela/ring/modular.h"
#include "lela/blas/level1.h"

#include "test-common.h"

#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;
using namespace LELA;

template <class Ring, class Modules>
bool testRandomDenseStream (LELA::Context<Ring, Modules> &ctx, size_t n, int iterations, size_t i, int p)
{
	ostringstream str;
	str << "Testing Random Dense Stream " << ends;
	commentator.start (str.str ().c_str ());

	typedef Modular<uint8> Ring1;
	typedef Ring1::Element Element;

	Ring1 F (p);
	Context<Ring1, GenericModule<Ring> > ctxtest (F);

	bool pass = true;
	bool eqmaster = true;

	for(size_t k = 1; k <= i; ++k)
	{
	        RandomDenseStream<Ring1, Vector<Ring1>::Dense> stream1 (ctxtest.F, n, iterations);
	        Vector<Ring1>::Dense v (n);
	        stream1 >> v;

		Vector<Ring1>::Dense w (n);

	 	if ( k != 1)
		{
			if( LELA::BLAS1::is_zero(ctxtest, v) && LELA::BLAS1::is_zero(ctxtest, w))
			{
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR1: Constructed Dense Vector is Zero" << endl;
				pass = false;
			}

			if ( eqmaster)
			{
				if( !LELA::BLAS1::equal(ctxtest, v, w) )
					eqmaster = false;
			}
		}

		BLAS1::copy(ctxtest, v, w);
		stream1.reset ();
	}
	
	if( eqmaster)
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR2: Constructed Dense Vectors are all equal" << endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Ring, class Modules>
bool testRandomSparseStream (LELA::Context<Ring, Modules> &ctx, size_t n, int iterations, size_t i, int p)
{
	ostringstream str;
	str << "Testing Random Sparse Stream " << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;
	bool eqmaster = true;

	typedef Modular<uint8> Ring1;
	typedef Ring1::Element Element;

	Ring1 F (p);
	Context<Ring1, GenericModule<Ring> > ctxtest (F);

	for(size_t k = 1; k <= i; ++k)
	{
	        RandomSparseStream<Ring1, Vector<Ring1>::Sparse> stream1 (ctxtest.F, 1, n);
	        Vector<Ring1>::Sparse v;
	        stream1 >> v;

		Vector<Ring1>::Sparse w;

		if ( k != 1)
		{
			if( LELA::BLAS1::is_zero(ctxtest, v) )
			{
			        if( LELA::BLAS1::is_zero(ctxtest, w))
				{
				        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					        << "ERROR1: Constructed Sparse Vector is Zero" << endl;
					pass = false;
				}
			}

			if ( eqmaster)
			{
				if( !LELA::BLAS1::equal(ctxtest, v, w) )
					eqmaster = false;
			}
		}

		BLAS1::copy(ctxtest, v, w);
		stream1.reset ();
	}
	
	if( eqmaster)
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR2: Constructed Sparse Vectors are all equal" << endl;
		pass = false;
	}

	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

template <class Ring, class Modules>
bool testStandardBasisStream (LELA::Context<Ring, Modules> &ctx, size_t n, int iterations, size_t i, int p)
{
	ostringstream str;
	str << "Testing Standard Basis Stream " << ends;
	commentator.start (str.str ().c_str ());

	bool pass = true;
	bool eqmaster = true;

	typedef Modular<uint8> Ring1;
	typedef Ring1::Element Element;

	Ring1 F (p);
	Context<Ring1, GenericModule<Ring> > ctxtest (F);

	for(size_t k = 1; k <= i; ++k)
	{
	        StandardBasisStream<Ring1, Vector<Ring1>::Dense> stream1 (ctxtest.F, n);
	        Vector<Ring1>::Dense v (n);
	        stream1 >> v;

		for(size_t j; j <= v.size (); ++j)
		{
		       if(v[j] == 1)
		       {
			     if(eqmaster)
			     {
			           eqmaster = false;
			     }
			     else
			     {
			           commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				        << "ERROR1: Constructed Standard Vectors are not Standard" << endl;
				   pass = false;
			     }
		       }
		
		}

		stream1.reset ();
	}
	
	commentator.stop (MSG_STATUS (pass), (const char *) 0, __FUNCTION__);

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static int p = 101;
	static size_t n = 33;
	static size_t i = 50;
	static int iterations = 1;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef Modular<uint8> Ring;
	typedef Ring::Element Element;

	Ring F (p);
	Context<Ring, GenericModule<Ring> > ctx (F);

	commentator.start ("Testing Stream", "Stream");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	pass =  testRandomDenseStream(ctx, n, iterations, i, p) && pass;
	pass =  testRandomSparseStream(ctx, n, iterations, i, p) && pass;
       	pass =  testStandardBasisStream(ctx, n, iterations, i, p) && pass;

	commentator.stop (MSG_STATUS (pass));
	return pass ? 0 : -1;
}
