# -------------------------------------------------------------
# Read/Write Compressed Streams with gzstream
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GZ],
[
  AC_ARG_ENABLE(gzstreams,
                AS_HELP_STRING([--disable-gzstreams],
                               [build without gzstreams compressed I/O support]),
                [AS_CASE("${enableval}",
                         [yes], [enablegz=yes],
                         [no],  [enablegz=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-gz)])],
                [enablegz=$enableoptional])

  dnl First check for the required system headers and libraries
  AS_IF([test "x$enablegz" = "xyes"],
        [
          AC_CHECK_HEADERS(zlib.h, have_zlib_h=yes)
          AC_CHECK_LIB(z, gzopen, have_libz=yes)
        ])

  AS_IF([test "$have_zlib_h" != yes || test "$have_libz" != yes], [enablegz=no])

  dnl If both tests succeeded, continue the configuration process.
  AS_IF([test "x$enablegz" = "xyes"],
        [
          GZSTREAM_INCLUDE="-I\$(top_srcdir)/contrib/gzstream"
          GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(libext) -lz"
          AC_DEFINE(HAVE_GZSTREAM, 1, [Flag indicating whether or not gzstreams are available])
          AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)
        ],
        [
          dnl Otherwise do not enable gzstreams
          GZSTREAM_INCLUDE=""
          GZSTREAM_LIB=""
        ])

  AC_SUBST(GZSTREAM_INCLUDE)
  AC_SUBST(GZSTREAM_LIB)
])
