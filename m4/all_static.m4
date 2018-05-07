# Configure option to pass -all-static to libtool's link mode.  So
# far, this was only required on the BG/Q machine Mira to get static
# linking to work.  Off by default.  Other than LDFLAGS, AM_LDFLAGS
# seems to be the only flag that makes it into the Makefile.in files
# after the --mode=link arguments to libtool.
AC_DEFUN([AX_ALL_STATIC],
[
  AC_ARG_ENABLE(all-static,
                [AS_HELP_STRING([--enable-all-static],[Pass -all-static to libtool's link mode])],
                enableallstatic=$enableval,
                enableallstatic=no)

  AS_IF([test "x$enableallstatic" = "xyes"],
        [
          AS_IF([test "x$libmesh_LDFLAGS" != "x"],
                [
                  dnl Append to whatever the user has for libmesh_LDFLAGS in their environment
                  libmesh_LDFLAGS="$libmesh_LDFLAGS -all-static"
                ],
                [
                  dnl Set libmesh_LDFLAGS to -all-static
                  libmesh_LDFLAGS=-all-static
                ])

          dnl Substitute into Makefiles where relevant.
          AC_SUBST(libmesh_LDFLAGS)
        ])
])
