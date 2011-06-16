/* util/support.C
 * Copyright 2010 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Support-routines for utilities: command-line-processing etc.
 */

#include <cstring>
#include <iostream>
#include <cstdlib>

#include "linbox/util/commentator.h"
#include "linbox/integer.h"

#define __SUPPORT_C

#include "support.h"

using namespace LinBox;

#if 0

// Force instantiations of get_matrix_type
template <> MatrixType get_matrix_type<Modular<uint8> > (const char *str);
template <> MatrixType get_matrix_type<Modular<uint16> > (const char *str);
template <> MatrixType get_matrix_type<Modular<uint32> > (const char *str);
template <> MatrixType get_matrix_type<Modular<float> > (const char *str);
template <> MatrixType get_matrix_type<Modular<double> > (const char *str);

#endif

Argument *findArgument (Argument *args, char c) 
{
	int i;

	for (i = 0; args[i].c != '\0' && args[i].c != c; i++) ;

	if (args[i].c != '\0')
		return &(args[i]);
	else
		return (Argument *) 0;
}

void parseArguments (int argc, char **argv, Argument *args, int freeArgs, ...)
{
	int i;
	Argument *current;
	va_list arg_list;

	va_start (arg_list, freeArgs);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 0) {
				commentator.setBriefReportStream (std::cout);
				commentator.setReportStream (std::cout);
				std::cout << "Writing report data to cout (intermingled with brief report)" << std::endl << std::endl;
				std::cout.flush ();
			}
#if 0 // FIXME
			else if (argv[i][1] == 'h' || argv[i][1] == '?') {
				printHelpMessage (argv[0], args, true);
				exit (1);
			}
#endif
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
					*(char **) current->data = argv[i+1];
					i++;
					break;
				}
			} else {
				std::cerr << "ERROR: Bad argument " << argv[i] << std::endl;
				break;
			}
		} else {
			if (freeArgs > 0) {
				*(va_arg (arg_list, char **)) = argv[i];
				--freeArgs;
			} else {
				commentator.setBriefReportStream(std::cout);
				commentator.setDefaultReportFile (argv[i]);
				std::cout << "Writing report data to " << argv[i] << std::endl << std::endl;
				std::cout.flush ();
			}
		}
	}

	va_end (arg_list);
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
