dnl -------------------------------------------------------------
dnl nemesis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NEMESIS],
[
  AC_ARG_ENABLE(nemesis,
                AS_HELP_STRING([--disable-nemesis],
                               [build without NemesisII API support]),
                [AS_CASE("${enableval}",
                         [yes|new|v522], [enablenemesis=yes ; nemesisversion="v5.22"],
                         [old|v309],     [enablenemesis=yes ; nemesisversion="v3.09"],
                         [no],           [enablenemesis=no  ; nemesisversion=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-nemesis)])],
                [
                  # if unspecified, depend on exodus
                  enablenemesis=$enableexodus ;
                  AS_IF([test "x$exodusversion" = "xv5.22"],
                        [nemesisversion="v5.22"],
                        [nemesisversion="v3.09"])
                ])

  # Trump --enable-nemesis with --disable-mpi
  AS_IF([test "x$enablempi" = xno],
        [nemesisversion=no])

  AS_CASE("${nemesisversion}",
          ["v3.09"],
          [
            dnl The NEMESIS API is distributed with libmesh, so we don't have to guess
            dnl where it might be installed...
            NEMESIS_INCLUDE="-I\$(top_srcdir)/contrib/nemesis/$nemesisversion"
            AC_DEFINE(HAVE_NEMESIS_API, 1, [Flag indicating whether the library will be compiled with Nemesis support])
            AC_MSG_RESULT(<<< Configuring library with Nemesis version $nemesisversion support >>>)
          ],
          ["v5.22"],
          [
            NEMESIS_INCLUDE="-I\$(top_srcdir)/contrib/nemesis/$nemesisversion/nemesis"
            AC_DEFINE(HAVE_NEMESIS_API, 1, [Flag indicating whether the library will be compiled with Nemesis support])
            AC_MSG_RESULT(<<< Configuring library with Nemesis version $nemesisversion support >>>)
          ],
          [
            NEMESIS_INCLUDE=""
            enablenemesis=no
          ])

  AC_CONFIG_FILES([contrib/nemesis/v3.09/Makefile])
  AC_CONFIG_FILES([contrib/nemesis/v5.22/nemesis/Makefile])
  AC_SUBST(NEMESIS_INCLUDE)
])
