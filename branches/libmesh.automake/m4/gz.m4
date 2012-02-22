dnl -------------------------------------------------------------
dnl Read/Write Compressed Streams with gzstream
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GZ], 
[
  AC_ARG_ENABLE(gzstreams,
                AC_HELP_STRING([--enable-gzstreams],
                               [build with gzstreams compressed I/O suppport]),
		[case "${enableval}" in
		  yes)  enablegz=yes ;;
		   no)  enablegz=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-gz) ;;
		 esac],
		 [enablegz=yes])
		 
	  

if (test $enablegz = yes); then
  dnl First check for the required system headers and libraries
  AC_CHECK_HEADERS(zlib.h, have_zlib_h=yes)
  AC_CHECK_LIB(z, gzopen, have_libz=yes)

  dnl If both tests succeded, continue the configuration process.
  if (test "$have_zlib_h" = yes -a "$have_libz" = yes) ; then
    GZSTREAM_INCLUDE="-I\$(top_srcdir)/contrib/gzstream"
    #GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(libext) -lz"
    AC_DEFINE(HAVE_GZSTREAM, 1, [Flag indicating whether or not gzstreams are available])
    AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)

  dnl Otherwise do not enable gzstreams
  else
    enablegz=no;
  fi
fi

AC_SUBST(GZSTREAM_INCLUDE)
#AC_SUBST(GZSTREAM_LIB)	
AC_SUBST(enablegz)

AM_CONDITIONAL(ENABLE_GZSTREAMS, test x$enablegz = xyes)
AC_CONFIG_FILES([contrib/gzstream/Makefile])
])
