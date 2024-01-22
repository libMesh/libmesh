dnl -------------------------------------------------------------
dnl NetGen mesh generation library
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NETGEN],
[
  AC_ARG_ENABLE(netgen,
                AS_HELP_STRING([--enable-netgen],
                               [build with NetGen mesh generation library support]),
                [AS_CASE("${enableval}",
                         [yes], [enablenetgen=yes],
                         [no],  [enablenetgen=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-netgen)])],
                [enablenetgen=no])

  dnl The Netgen API is distributed with libmesh, so we don't
  dnl currently have to guess where it might be installed, though we
  dnl ought to add support later for linking to a pre-installed system
  dnl Netgen.
  NETGEN_INCLUDE=""
  NETGEN_LIBS=""
  NETGEN_BUILD_LDFLAGS=""
  AS_IF([test "x$enablenetgen" = "xyes"],
        [
          dnl Autoconf doesn't define $abs_top_srcdir at this point;
          dnl here's a trick from GraphicsMagick:
          my_top_srcdir="$(cd $srcdir && pwd)"

          dnl I thought abs_top_builddir was defined here but
          dnl apparently not?
          my_top_builddir="$(pwd)"

          dnl Wipe the old build directory, because it can corrupt a
          dnl "fresh" cmake, because cmake sucks
          rm -rf contrib/netgen/build
          mkdir -p contrib/netgen/build

          dnl Make a fake install directory, because we need that
          dnl prefix to be created and writable or we get a grossly
          dnl misleading error message, because cmake sucks
          mkdir -p contrib/netgen/install

          dnl Refer to the install directory by an absolute rather
          dnl than a relative path, because cmake emits errors if
          dnl given the latter.  Guess why?
          AS_IF([(cd contrib/netgen/build && \
                   cmake -DCMAKE_INSTALL_PREFIX:PATH=$my_top_builddir/contrib/netgen/install \
                   $my_top_srcdir/contrib/netgen/netgen)],
                [
                  NETGEN_INCLUDE="-I\$(top_srcdir)/contrib/netgen/ -I\$(top_builddir)/contrib/netgen/build/"
                  NETGEN_LIBS="-lnglib -lngcore"
                  NETGEN_BUILD_LDFLAGS="-L\$(abs_top_builddir)/contrib/netgen/build/netgen/ -L\$(abs_top_builddir)/contrib/netgen/build/netgen/libsrc/core/"
                  AS_IF([test "x$RPATHFLAG" != "x"],
                        [NETGEN_BUILD_LDFLAGS="$NETGEN_BUILD_LDFLAGS ${RPATHFLAG}\$(abs_top_builddir)/contrib/netgen/build/netgen/ ${RPATHFLAG}\$(abs_top_builddir)/contrib/netgen/build/netgen/libsrc/core/"]);
                  AC_DEFINE(HAVE_NETGEN, 1, [Flag indicating whether the library will be compiled with Netgen support])
                  AC_MSG_RESULT(<<< Configuring library with Netgen support >>>)
                ],
                [
                  AC_MSG_RESULT(<<< Failed cmake configuration of Netgen submodule >>>)
                  enablenetgen=no
                ])
        ],
        [
          enablenetgen=no
        ])

  AC_SUBST(NETGEN_INCLUDE)
  AC_SUBST(NETGEN_LIBS)
  AC_SUBST(NETGEN_BUILD_LDFLAGS)
])
