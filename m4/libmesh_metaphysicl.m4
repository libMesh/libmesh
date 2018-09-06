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

  dnl Check for existence of a file that should always be in metaphysicl
  AS_IF([test "x$enablemetaphysicl" = "xyes" && test -r $top_srcdir/contrib/metaphysicl/README],
        [enablemetaphysicl=yes],
        [
          AC_MSG_RESULT([>>> Configuring metaphysicl failed, you may need to run 'git submodule update --init' first <<<])
          enablemetaphysicl=no
        ])

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
