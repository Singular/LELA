# Check for libm4ri
# Bradford Hovinen, 2011-01-25
# Modified by Pascal Giorgi, 2003-12-03
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_M4RI ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for libm4ri and define M4RI_CFLAGS and M4RI_LIBS

AC_DEFUN([LB_CHECK_M4RI],
[

AC_ARG_WITH(m4ri,
[  --with-m4ri=<path>|yes Use libm4ri. This is a library for extremely fast
   			 calculations with dense matrices over GF(2). Enabling
			 it creates a wrapper for this library in LELA so
			 that it can be used. Without it, calculations over
			 GF(2) are much slower.
],		  
	     [if test "$withval" = yes ; then
			M4RI_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			M4RI_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
	     fi],
	     [M4RI_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl min_m4ri_version=ifelse([$1], ,3.2.12,$1)


dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for M4RI)

for M4RI_HOME in ${M4RI_HOME_PATH} 
 do	
if test -r "$M4RI_HOME/include/m4ri/m4ri.h"; then
	if test "x$M4RI_HOME" != "x/usr"; then
		M4RI_CFLAGS="-I${M4RI_HOME}/include"
		M4RI_LIBS="-L${M4RI_HOME}/lib -Wl,-rpath,${M4RI_HOME}/lib -lm4ri"
	else
		M4RI_CFLAGS=
		M4RI_LIBS="-L${M4RI_HOME}/lib -Wl,-rpath,${M4RI_HOME}/lib -lm4ri"
	fi	
	CXXFLAGS="${BACKUP_CXXFLAGS} ${M4RI_CFLAGS} ${GMP_CFLAGS}" 
	LIBS="${BACKUP_LIBS} ${M4RI_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK(
	[#include <m4ri/m4ri.h>],
	[mzd_t *matrix;],
	[
	AC_TRY_RUN(
	[#include <m4ri/m4ri.h>	 
	 int main () { return 0; }
	],[
	m4ri_found="yes"	
	break
	],[	
	m4ri_problem="$problem $M4RI_HOME"	
	unset M4RI_CFLAGS
	unset M4RI_LIBS
	],[
	m4ri_found="yes"
	m4ri_cross="yes"
	
	break
	])	
	],
	[
	m4ri_found="no"
	m4ri_checked="$checked $M4RI_HOME"
	unset M4RI_CFLAGS
	unset M4RI_LIBS
	
	])
else
	m4ri_found="no"
fi	
done

if test "x$m4ri_found" = "xyes" ; then		
	AC_SUBST(M4RI_CFLAGS)
	AC_SUBST(M4RI_LIBS)
	AC_DEFINE(HAVE_M4RI,1,[Define if M4RI is installed])
	HAVE_M4RI=yes
	if test "x$m4ri_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your M4RI version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])
elif test -n "$m4ri_problem"; then
	AC_MSG_RESULT(problem)
	echo "Sorry, your M4RI version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test "x$m4ri_found" = "xno" ; then
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi

dnl Check m4ri version
if test "x$m4ri_found" = "xyes" ; then
   	AC_MSG_CHECKING(whether version of M4RI is at least 20110601)
	AC_LANG_PUSH([C])
	BACKUP_CFLAGS=${CFLAGS}
	CFLAGS="${CFLAGS} -std=c99 -I${M4RI_HOME}/include"
	AC_COMPILE_IFELSE(
		AC_LANG_PROGRAM([#include <m4ri/m4ri.h>],[mzd_t M; wi_t i = M.rowstride;]),
		[m4ri_new_version=yes],
		[m4ri_new_version=no])

	AC_LANG_POP
	CFLAGS=${BACKUP_CFLAGS}
	if test "x$m4ri_new_version" = "xyes"; then
	   	AC_DEFINE_UNQUOTED(HAVE_M4RI_GE_20110601,1,[Version of M4RI is at least 20110601])
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
	fi
fi

AM_CONDITIONAL(LELA_HAVE_M4RI, test "x$HAVE_M4RI" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
