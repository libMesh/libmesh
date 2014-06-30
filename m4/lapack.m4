# ----------------------------------------------------------------------------
# @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# This macro looks for a library that implements the LAPACK
# linear-algebra interface (see http://www.netlib.org/lapack/).
# On success, it sets the LAPACK_LIBS output variable to
# hold the requisite library linkages.
#
# To link with LAPACK, you should link with:
#
#     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
#
# in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
# macro, called automatically.  FLIBS is the output variable of the
# AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
# and is sometimes necessary in order to link with F77 libraries.
# Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
# manual), for the same reason.
#
# The user may also use --with-lapack=<lib> in order to use some
# specific LAPACK library <lib>.  In order to link successfully,
# however, be aware that you will probably need to use the same
# Fortran compiler (which can be set via the F77 env. var.) as
# was used to compile the LAPACK and BLAS libraries.
#
# ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
# library is found, and ACTION-IF-NOT-FOUND is a list of commands
# to run it if it is not found.  If ACTION-IF-FOUND is not specified,
# the default action will define HAVE_LAPACK.
#
# @version acsite.m4,v 1.3 2002/08/02 09:28:12 steve Exp
# @author Steven G. Johnson <stevenj@alum.mit.edu>

AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
            AS_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>]))
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])# ACX_LAPACK
