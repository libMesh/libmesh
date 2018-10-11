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
          AC_CONFIG_SUBDIRS([contrib/metaphysicl])
        ],
        [
          METAPHYSICL_INCLUDE=""
          enablemetaphysicl=no
        ])

  AC_SUBST(METAPHYSICL_INCLUDE)
])
