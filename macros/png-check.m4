# Check for PNG

dnl LB_CHECK_PNG
dnl
dnl Test for the libpng and define PNG_CFLAGS, PNG_LIBS, and HAVE_LIBPNG

AC_DEFUN([LB_CHECK_PNG],
[
PNG_LIBS=
PNG_CFLAGS=

AC_ARG_WITH(png,[
   --with-png= <path>|yes Enable use of PNG library.
],[
if test "$withval" = yes ; then
    AC_CHECK_LIB([png], [png_create_write_struct_2], [
        AC_DEFINE(HAVE_LIBPNG,1,Enable use of libpng)
    	PNG_LIBS="-lpng"
    ])
elif test "$withval" != no ; then
    BACKUP_LIBS=$LIBS
    LIBS=$withval/lib
    AC_CHECK_LIB([png], [png_create_write_struct_2],[
            AC_DEFINE(HAVE_LIBPNG,1,Enable use of libpng)
	    PNG_CFLAGS="-I$withval/include -lpng"
	    PNG_LIBS="-L$withval/lib -lpng"
    ],[
	    AC_MSG_WARN([Could not find libpng in path $withval, disabling])
    ])
    LIBS=$BACKUP_LIBS
fi
],[
    AC_CHECK_LIB([png], [png_create_write_struct_2], [
    	AC_DEFINE(HAVE_LIBPNG,1,Enable use of libpng)
    	PNG_LIBS="-lpng"
    ],[AC_MSG_RESULT(not found)])
])

AC_SUBST(PNG_CFLAGS)
AC_SUBST(PNG_LIBS)
])
