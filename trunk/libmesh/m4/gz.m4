# -------------------------------------------------------------
# Read/Write Compressed Streams with gzstream
# -------------------------------------------------------------
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
		 [enablegz=$enableoptional])
		 
	  

if (test $enablegz = yes); then
  # First check for the required system headers and libraries
  AC_CHECK_HEADERS(zlib.h, have_zlib_h=yes)
  AC_CHECK_LIB(z, gzopen, have_libz=yes)

  # If both tests succeded, continue the configuration process.
  if (test "$have_zlib_h" = yes -a "$have_libz" = yes) ; then
    GZSTREAM_INCLUDE="-I\$(top_srcdir)/contrib/gzstream"
    #GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(libext) -lz"
    AC_DEFINE(HAVE_GZSTREAM, 1, [Flag indicating whether or not gzstreams are available])
    AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)     
    libmesh_contrib_INCLUDES="$GZSTREAM_INCLUDE $libmesh_contrib_INCLUDES"
    libmesh_optional_LIBS="-lz $libmesh_optional_LIBS"
  # Otherwise do not enable gzstreams
  else
    enablegz=no;
  fi
fi

AC_SUBST(GZSTREAM_INCLUDE)
#AC_SUBST(GZSTREAM_LIB)	
AC_SUBST(enablegz)

AM_CONDITIONAL(LIBMESH_ENABLE_GZSTREAMS, test x$enablegz = xyes)
])
