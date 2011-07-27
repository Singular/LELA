/* lela/tests/test-common.C
 * Copyright 2001, 2002 Bradford Hovinen
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
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "lela/util/commentator.h"
#include "lela/ring/interface.h"

#include "test-common.h"

using namespace LELA;

/* Display a help message on command usage */

void printHelpMessage (const char *program, Argument *args, bool printDefaults) 
{
	int i, l;

	// Skip past libtool prefix in program name
	if (!strncmp (program, "lt-", strlen ("lt-")))
		program += strlen ("lt-");

	std::cout << "Usage: " << program << " [options] [<report file>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Where [options] are the following:" << std::endl;

	for (i = 0; args[i].c != '\0'; i++) {
		if (args[i].example != 0) {
			std::cout << "  " << args[i].example;
			l = 10 - strlen (args[i].example);
			do std::cout << ' '; while (--l > 0);
		}
		else if (args[i].type == TYPE_NONE)
			std::cout << "  -" << args[i].c << " {YN+-} ";
		else 
			std::cout << "  -" << args[i].c << ' ' << args[i].c << "      ";
			
		std::cout << args[i].helpString;
		if (printDefaults) {
			l = 54 - strlen (args[i].helpString);
			do std::cout << ' '; while (--l > 0);
			std::cout << " (default ";
			switch (args[i].type) {
			case TYPE_NONE:
				cout << ((*(bool *)args[i].data)?"ON":"OFF");
				break;
			case TYPE_INT:
				cout << *(int *) args[i].data;
				break;
			case TYPE_INTEGER:
				cout << *(integer *) args[i].data;
				break;
			case TYPE_DOUBLE:
				cout << *(double *) args[i].data;
				break;
			case TYPE_STRING:
				cout << args[i].data;
				break;
			}
			std::cout << ")";		
		}
		std::cout << std::endl;
	}

	std::cout << "  -h or -?  Display this message" << std::endl;
	std::cout << "For boolean switches, the argument may be omitted, meaning the switch should be ON" << std::endl;
	std::cout << std::endl;
	std::cout << "If <report file> is '-' the report is written to std output.  If <report file> is" << std::endl; 
	std::cout << "not given, then no detailed reporting is done. This is suitable if you wish only" << std::endl;
	std::cout << "to determine whether the tests succeeded." << std::endl;
	std::cout << std::endl;
	std::cout << "[1] N.B. This program does not verify the primality of Q, and does not use a" << std::endl;
	std::cout << "    field extension in the event that Q=p^n, n > 1" << std::endl;
	std::cout << std::endl;
}

/* Find an argument in the argument list for a character */

Argument *findArgument (Argument *args, char c) 
{
	int i;

	for (i = 0; args[i].c != '\0' && args[i].c != c; i++) ;

	if (args[i].c != '\0')
		return &(args[i]);
	else
		return (Argument *) 0;
}

/* Parse command line arguments */

void parseArguments (int argc, char **argv, Argument *args, bool printDefaults)
{
	int i;
	Argument *current;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 0) {
				commentator.setBriefReportStream (cout);
				commentator.setReportStream (cout);
				std::cout << "Writing report data to cout (intermingled with brief report)" << std::endl << std::endl;
				std::cout.flush ();
			}
			else if (argv[i][1] == 'h' || argv[i][1] == '?') {
				printHelpMessage (argv[0], args, printDefaults);
				exit (1);
			}
			else if ((current = findArgument (args, argv[i][1])) != (Argument *) 0) {
				switch (current->type) {
				case TYPE_NONE:
					if (argc == i+1 || (argv[i+1][0] == '-' && argv[i+1][1] != '\0')) {
						// if at last argument, or next argument is a switch, set to true
						*(bool *) current->data = true;
						break;
					}
					*(bool *) current->data = 
						(argv[i+1][0] == '+' 
						 || argv[i+1][0] == 'Y' 
						 || argv[i+1][0] == 'y' 
						 || argv[i+1][0] == 'T' 
						 || argv[i+1][0] == 't') ;
					i++;
					break;

				case TYPE_INT:
					*(int *) current->data = atoi (argv[i+1]);
					i++;
					break;

				case TYPE_INTEGER: 
					{
						integer tmp (argv[i+1]);
						*(integer *) current->data = tmp;
					}
					i++;
					break;

				case TYPE_DOUBLE:
					*(double *) current->data = atof (argv[i+1]);
					i++;
					break;

				case TYPE_STRING:
					current->data = argv[i+1];
					i++;
					break;
				}
			} else {
				std::cerr << "ERROR: Bad argument " << argv[i] << std::endl;
				break;
			}
		} else {
			commentator.setBriefReportStream(cout);
			commentator.setDefaultReportFile (argv[i]);
			std::cout << "Writing report data to " << argv[i] << std::endl << std::endl;
			std::cout.flush ();
		}
	}
}

std::ostream& writeCommandString (std::ostream& os, Argument *args, char* programName)
{
	os << programName;

	for (int i = 0; args[i].c != '\0'; i++) {
		cout << " -" << args[i].c;

		switch (args[i].type) {
		case TYPE_NONE:
			if (! (*(bool *) args[i].data)) os << " N";
			break;
		case TYPE_INT:
			os << ' ' << *(int *) args[i].data;
			break;
		case TYPE_INTEGER:
			os << ' ' << *(integer *) args[i].data;
			break;
		case TYPE_DOUBLE:
			os << ' ' << *(double *) args[i].data;
			break;
		case TYPE_STRING:
			os << ' ' << args[i].data;
			break;
		}
	}

	return os << std::endl;
}

bool isPower (integer n, integer m)
{
	return (n == 1) || (((n % m) == 0) && isPower (n/m, m));
}

inline double incompleteGamma (double a, double x, double tol) 
{
	double xa_ex = pow (x, a) * exp (-x);
	double pi = 1.0;
	double xn = 1.0;
	double sigma = 0.0;
	double last_sigma;

	int n = 0;

	do {
		pi *= a + n;
		last_sigma = sigma;
		sigma += xn / pi;
		xn *= x;
		++n;
	} while (abs (sigma - last_sigma) >= tol) ;

	return sigma * xa_ex;
}

double chiSquaredCDF (double chi_sqr, double df)
{
	return incompleteGamma (df / 2.0, chi_sqr / 2.0, 1e-10) / exp (gamma (df / 2.0));
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
