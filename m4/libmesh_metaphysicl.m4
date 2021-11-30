# -------------------------------------------------------------
# MetaPhysicL - a C++ header-only numerics utility library
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METAPHYSICL],
[
  AC_ARG_ENABLE(metaphysicl,
                AS_HELP_STRING([--disable-metaphysicl],
                               [build without MetaPhysicL suppport]),
                [AS_CASE("${enableval}",
                         [yes], [enablemetaphysicl=yes],
                         [no],  [enablemetaphysicl=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-metaphysicl)])],
                [enablemetaphysicl=$enablenested])

  dnl Setting --enable-metaphysicl-required causes an error to be
  dnl emitted during configure if the MetaPhysicL library is not
  dnl detected successfully.  This is useful for app codes which require
  dnl MetaPhysicL (like MOOSE-based apps), since it prevents situations
  dnl where libmesh is accidentally built without MetaPhysicL support
  dnl (which may take a very long time), and then the app fails to
  dnl compile, requiring you to redo everything.
  AC_ARG_ENABLE(metaphysicl-required,
                AC_HELP_STRING([--enable-metaphysicl-required],
                               [Error if MetaPhysicL is not detected by configure]),
                [AS_CASE("${enableval}",
                         [yes], [metaphysiclrequired=yes],
                         [no],  [metaphysiclrequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-metaphysicl-required)])],
                     [metaphysiclrequired=no])

  dnl Let --enable-metaphysicl-required override the value of $enablenested. Also, if the user
  dnl provides the nonsensical combination "--disable-metaphysicl --enable-metaphysicl-required"
  dnl we'll set $enablemetaphysicl to yes instead.
  AS_IF([test "x$enablemetaphysicl" = "xno" && test "x$metaphysiclrequired" = "xyes"],
        [enablemetaphysicl=yes])

  dnl Check for existence of a file that should always be in metaphysicl
  AS_IF([test "x$enablemetaphysicl" = "xyes" && test -r $top_srcdir/contrib/metaphysicl/README],
        [enablemetaphysicl=yes],
        [
          AC_MSG_RESULT([>>> Configuring metaphysicl failed, you may need to run 'git submodule update --init' first <<<])
          enablemetaphysicl=no
        ])

  dnl If metaphysicl was required but isn't available, throw an error.
  dnl We return a non-unity error code here, since 0 means success and 1 is
  dnl indistinguishable from other errors.  Ideally, all of the
  dnl AC_MSG_ERROR calls in our m4 files would return a different
  dnl error code, but currently this is not implemented.
  AS_IF([test "x$enablemetaphysicl" = "xno" && test "x$metaphysiclrequired" = "xyes"],
        [AC_MSG_ERROR([*** MetaPhysicL was not found, but --enable-metaphysicl-required was specified.], 5)])

  dnl The MetaPhysicL API is distributed with libmesh, so we don't have
  dnl to guess where it might be installed.  This needs to be replaced
  dnl someday with an option to include an external version instead.
  AS_IF([test "x$enablemetaphysicl" = "xyes"],
        [
          METAPHYSICL_INCLUDE="-I\$(top_srcdir)/contrib/metaphysicl/src/numerics/include -I\$(top_srcdir)/contrib/metaphysicl/src/core/include -I\$(top_srcdir)/contrib/metaphysicl/src/utilities/include"
          AC_DEFINE(HAVE_METAPHYSICL, 1, [Flag indicating whether the library will be compiled with MetaPhysicL support])
          AC_MSG_RESULT(<<< Configuring library with MetaPhysicL support >>>)

          dnl MetaPhysicL needs to know which TIMPI library to link to
          dnl for "make check"
          my_method=dbg
          AS_IF([test "x$build_devel" = xyes],
                [my_method=devel])
          AS_IF([test "x$build_prof" = xyes],
                [my_method=prof])
          AS_IF([test "x$build_oprof" = xyes],
                [my_method=oprof])
          AS_IF([test "x$build_opt" = xyes],
                [my_method=opt])

          dnl Autoconf doesn't define $abs_top_srcdir at this point;
          dnl here's a trick from GraphicsMagick:
          my_top_srcdir="$(cd $srcdir && pwd)"

          dnl FIXME: setting TIMPI_DIR, even to something invalid, is
          dnl currently the best way to get MetaPhysicL to recognize
          dnl that we really want TIMPI.

          AX_SUBDIRS_CONFIGURE([contrib/metaphysicl],[[--with-cxx-std=20$acsm_cxx_version],[CPPFLAGS=-I$my_top_srcdir/contrib/timpi/src/algorithms/include/ -I$my_top_srcdir/contrib/timpi/src/parallel/include/ -I$my_top_srcdir/contrib/timpi/src/utilities/include/ -I$ac_abs_top_builddir/contrib/timpi/src/utilities/include/ $CPPFLAGS],[LDFLAGS=-L$ac_abs_top_builddir/contrib/timpi/src/ $LDFLAGS],[TIMPI_DIR=$my_top_srcdir/contrib/timpi/],[--with-timpi-method=$my_method]])
        ],
        [
          METAPHYSICL_INCLUDE=""
          enablemetaphysicl=no
        ])

  AC_SUBST(METAPHYSICL_INCLUDE)
])
