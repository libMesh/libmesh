# -------------------------------------------------------------
# Read/Write Compressed Streams with gzstream
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GZ],
[
  AC_ARG_ENABLE(gzstreams,
                AS_HELP_STRING([--disable-gzstreams],
                               [build without gzstreams compressed I/O support]),
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
fi

if (test "$have_zlib_h" != yes -o "$have_libz" != yes) ; then
  enablegz=no
fi

# If both tests succeeded, continue the configuration process.
if (test "$enablegz" = yes) ; then
  GZSTREAM_INCLUDE="-I\$(top_srcdir)/contrib/gzstream"
  GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(libext) -lz"
  AC_DEFINE(HAVE_GZSTREAM, 1, [Flag indicating whether or not gzstreams are available])
  AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)
# Otherwise do not enable gzstreams
else
  GZSTREAM_INCLUDE=""
  GZSTREAM_LIB=""
fi

AC_SUBST(GZSTREAM_INCLUDE)
AC_SUBST(GZSTREAM_LIB)
])
