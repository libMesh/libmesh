dnl -------------------------------------------------------------
dnl Read/Write Compressed Streams with gzstream
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GZ], 
[

dnl Initialize variables
GZSTREAM_INCLUDE=""
GZSTREAM_LIB=""

if (test $enablegz = yes); then
  dnl First check for C++-compatible system headers and libraries,
  dnl then to make sure the user really has the contrib directory
  AC_LANG_PUSH([C++])
  OLD_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="-Icontrib/gzstream $CPPFLAGS"
  AC_CHECK_HEADERS(zlib.h,
    [AC_CHECK_LIB(z, gzopen,
      [AC_CHECK_HEADER(gzstream.h,
                       [enablegz=yes], [enablegz=no])])])
  CPPFLAGS=$OLD_CPPFLAGS
  AC_LANG_POP([C++])
fi


if (test $enablegz = yes); then

  dnl If both tests succeded, continue the configuration process.
  if (test "$have_zlib_h" = yes -a "$have_libz" = yes) ; then
    GZSTREAM_INCLUDE="-I$PWD/contrib/gzstream"
    GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(libext) -lz"
    AC_DEFINE(HAVE_GZSTREAM, 1, [Flag indicating whether or not gzstreams are available])
    AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)

  dnl Otherwise do not enable gzstreams
  else
    enablegz=no;
  fi
fi

AC_SUBST(GZSTREAM_INCLUDE)
AC_SUBST(GZSTREAM_LIB)	
AC_SUBST(enablegz)

])
