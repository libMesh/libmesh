# -------------------------------------------------------------
# qhull
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_QHULL],
[
  AC_ARG_ENABLE(qhull,
                AS_HELP_STRING([--disable-qhull],
                               [build without Qhull API support]),
                [AS_CASE("${enableval}",
                         [yes], [enableqhull=yes],
                         [no],  [enableqhull=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-qhull)])],
                [enableqhull=$enableoptional]) # if unspecified, depend on enableoptional

  dnl Setting --enable-qhull-required causes an error to be emitted during
  dnl configure if Qhull is not successfully configured.
  AC_ARG_ENABLE(qhull-required,
                AS_HELP_STRING([--enable-qhull-required],
                               [Error if Qhull is not detected by configure]),
                [AS_CASE("${enableval}",
                         [yes], [enableqhullrequired=yes
                                 enableqhull=yes],
                         [no],  [enableqhullrequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-qhull-required)])],
                [enableqhullrequired=no])

  dnl The Qhull library is distributed with libmesh as a git submodule.
  QHULL_INCLUDE=""
  QHULL_LIBS=""
  QHULL_BUILD_LDFLAGS=""
  AS_IF([test "x$enableqhull" = "xyes"],
        [
          dnl Autoconf doesn't define $abs_top_srcdir at this point;
          dnl here's a trick from GraphicsMagick:
          my_top_srcdir="$(cd $srcdir && pwd)"

          dnl I thought abs_top_builddir was defined here but
          dnl apparently not?
          my_top_builddir="$(pwd)"

          dnl Wipe the old build directory (along with our .buildstamp
          dnl that lets us know a build is complete there), because an
          dnl old build can corrupt a "fresh" cmake, because cmake
          dnl sucks
          rm -f contrib/qhull/.buildstamp
          rm -rf contrib/qhull/build
          mkdir -p contrib/qhull/build

          dnl Make a fake install directory, because we need that
          dnl prefix to be created and writable or we get a grossly
          dnl misleading error message, because cmake sucks
          mkdir -p contrib/qhull/install

          dnl Build static libraries only to avoid RPATH issues in
          dnl the build tree.  Disable tests and applications to
          dnl keep the build lean.  Force -fPIC on the static
          dnl archives so they can be linked into libmesh's shared
          dnl library; qhull's CMakeLists only sets
          dnl POSITION_INDEPENDENT_CODE on libqhullcpp, not on
          dnl libqhullstatic_r.
          AS_IF([(cd contrib/qhull/build && \
                   cmake -DBUILD_SHARED_LIBS=OFF \
                         -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
                         -DQHULL_ENABLE_TESTING=OFF \
                         -DBUILD_APPLICATIONS=OFF \
                         -DCMAKE_INSTALL_PREFIX:PATH=$my_top_builddir/contrib/qhull/install \
                   $my_top_srcdir/contrib/qhull/qhull)],
                [
                  dnl Headers live entirely in the source tree; no
                  dnl cmake-generated headers to worry about.
                  QHULL_INCLUDE="-I\$(top_srcdir)/contrib/qhull/qhull/src"
                  QHULL_LIBS="-lqhullcpp -lqhullstatic_r"
                  QHULL_BUILD_LDFLAGS="-L\$(abs_top_builddir)/contrib/qhull/build/"
                  AC_DEFINE(HAVE_QHULL_API, 1, [Flag indicating whether the library will be compiled with Qhull support])
                  AC_MSG_RESULT(<<< Configuring library with Qhull support >>>)
                ],
                [
                  AS_IF([test "x$enableqhullrequired" = "xyes"],
                        [AC_MSG_ERROR(<<< Failed cmake configuration of user-required Qhull submodule >>>)],
                        [AC_MSG_RESULT(<<< Failed cmake configuration of Qhull submodule >>>)])

                  enableqhull=no
                ])
        ],
        [
          enableqhull=no
        ])

  AC_SUBST(QHULL_INCLUDE)
  AC_SUBST(QHULL_LIBS)
  AC_SUBST(QHULL_BUILD_LDFLAGS)
])
