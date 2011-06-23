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
#include <cstdarg>

#include "linbox/util/commentator.h"
#include "linbox/integer.h"

#define __SUPPORT_C

#include "support.h"

using namespace LinBox;

const char *format_names[] = { "detect", "unknown", "Turner", "one-based", "Dumas", "Maple", "Matlab", "Sage", "pretty", "PNG" };

#if 0

// Force instantiations of get_matrix_type
template <> MatrixType get_matrix_type<Modular<uint8> > (const char *str);
template <> MatrixType get_matrix_type<Modular<uint16> > (const char *str);
template <> MatrixType get_matrix_type<Modular<uint32> > (const char *str);
template <> MatrixType get_matrix_type<Modular<float> > (const char *str);
template <> MatrixType get_matrix_type<Modular<double> > (const char *str);

#endif

/* Strip leading text from program-name */

const char *stripProgram (const char *program)
{
	const char *part = strstr (program, "lt-");

	if (part != NULL)
		return part + strlen ("lt-");
	else
		return program;
}

/* Display a help message on command usage */

void printHelpMessage (const char *program, Argument *args, const char *freeArgsText, bool printDefaults) 
{
	int i, l;

	// Skip past libtool prefix in program name
	program = stripProgram (program);

	std::cout << "Usage: " << program << " [options] " << freeArgsText << " [<report file>]" << std::endl;
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
				std::cout << ((*(bool *)args[i].data)?"ON":"OFF");
				break;
			case TYPE_INT:
				std::cout << *(int *) args[i].data;
				break;
			case TYPE_INTEGER:
				std::cout << *(integer *) args[i].data;
				break;
			case TYPE_DOUBLE:
				std::cout << *(double *) args[i].data;
				break;
			case TYPE_STRING:
				std::cout << *(char **) args[i].data;
				break;
			}
			std::cout << ")";		
		}
		std::cout << std::endl;
	}

	std::cout << "  -v        Verbose output" << std::endl;
	std::cout << "  -h or -?  Display this message" << std::endl;
	std::cout << std::endl;
	std::cout << "If <report file> is '-' the report is written to std output.  If <report file> is" << std::endl; 
	std::cout << "not given, then the detailed report is not written." << std::endl;
	std::cout << std::endl;
}

Argument *findArgument (Argument *args, char c) 
{
	int i;

	for (i = 0; args[i].c != '\0' && args[i].c != c; i++) ;

	if (args[i].c != '\0')
		return &(args[i]);
	else
		return (Argument *) 0;
}

void parseArguments (int argc, char **argv, Argument *args, const char *freeArgsText, int freeArgs, ...)
{
	int i;
	Argument *current;
	va_list arg_list;

	va_start (arg_list, freeArgs);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 0) {
				commentator.setReportStream (std::cout);
				std::cout << "Writing report data to cout" << std::endl << std::endl;
				std::cout.flush ();
			}
			if (argv[i][1] == 'v') {
				commentator.setBriefReportStream (std::cout);
				commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, true, true, true);
				commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (4);
				commentator.getMessageClass (BRIEF_REPORT).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
				std::cout.flush ();
			}
			else if (argv[i][1] == 'h' || argv[i][1] == '?') {
				printHelpMessage (argv[0], args, freeArgsText, true);
				exit (1);
			}
			else if ((current = findArgument (args, argv[i][1])) != (Argument *) 0) {
				switch (current->type) {
				case TYPE_NONE:
					// if at last argument, or next argument is a switch, set to true
					*(bool *) current->data = true;
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
