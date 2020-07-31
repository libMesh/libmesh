# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_CONFIGURE_OPTIONAL_PACKAGES],
[
dnl We need to make sure that we've done AC_ARG_ENABLE(optional), AC_ARG_ENABLE(mpi), AC_ARG_WITH(mpi),
dnl and AC_ARG_ENABLE(petsc) which all occur in ACSM_COMPILER_CONTROL_ARGS
AC_REQUIRE([ACSM_COMPILER_CONTROL_ARGS])

dnl We also need to ensure that we've set our compilers, which is where query for a valid
dnl PETSc configuration
AC_REQUIRE([LIBMESH_SET_COMPILERS])

# initialize these empty - append below
# note that
# libmesh_optional_INCLUDES and
# libmesh_optional_LIBS should point to third party packages
# outside the libMesh source and installation tree, and will
# be exported to the installation environment.
#
# By contrast, libmesh_contrib_INCLUDES point inside the
# source tree for building contributed packages that do not
# need to be exported as part of the installation environment.
#
# libmesh_subpackage_arguments is a list of configure arguments
# that will be passed down to any subpackages that we are nesting.
#
# libmesh_pkgconfig_requires is a list of pkgconfig requirements
# we will add
#
# libmesh_installed_LIBS are libraries that we install and need to
# link with, usually not needed.  presently only for Tecplot's binary
# library blob
libmesh_optional_INCLUDES=""
libmesh_optional_LIBS=""
libmesh_contrib_INCLUDES=""
libmesh_subpackage_arguments=""
libmesh_pkgconfig_requires=""
libmesh_installed_LIBS=""

# --------------------------------------------------------------
# TIMPI is required - so we should be able to see our TIMPI submodule.
# If our version of TIMPI doesn't have AM MAINTAINER MODE set, then we
# need its autoconf submodule initialized too.
# --------------------------------------------------------------
AS_IF([test -r $top_srcdir/contrib/timpi/README &&
       (test -r $top_srcdir/contrib/timpi/m4/autoconf-submodule/acsm_mpi.m4 ||
        grep "^AM""_MAINTAINER_MODE" $top_srcdir/contrib/timpi/configure.ac >/dev/null)],
[
  libmesh_contrib_INCLUDES="-I\$(top_srcdir)/contrib/timpi/src/algorithms/include $libmesh_contrib_INCLUDES"
  libmesh_contrib_INCLUDES="-I\$(top_srcdir)/contrib/timpi/src/parallel/include $libmesh_contrib_INCLUDES"
  libmesh_contrib_INCLUDES="-I\$(top_srcdir)/contrib/timpi/src/utilities/include $libmesh_contrib_INCLUDES"
  # Including timpi_config.h
  libmesh_contrib_INCLUDES="-I\$(top_builddir)/contrib/timpi/src/utilities/include $libmesh_contrib_INCLUDES"
],
[
  AC_MSG_ERROR([You must run "git submodule update --init --recursive" before configuring libmesh])
])


# Note that even when optional packages are disabled we need to
# run their m4 macros to get proper AM_CONDITIONALs.  Just be
# quiet about it...
AS_IF([test "$enableoptional" != no],
      [
        AC_MSG_RESULT(---------------------------------------------)
        AC_MSG_RESULT(----- Configuring for optional packages -----)
        AC_MSG_RESULT(---------------------------------------------)
      ])


# --------------------------------------------------------------
# Allow user to specify --disable-strict-lgpl
# By default libmesh is built only with LGPL-compatible contrib
# libraries, but the user can pass --disable-strict-lgpl to configure
# to turn on Laspack, Triangle, and Space-filling curves library
# support.
# --------------------------------------------------------------
AC_ARG_ENABLE(strict-lgpl,
              AS_HELP_STRING([--disable-strict-lgpl],
                             [Compile libmesh with even non-LGPL-compatible contrib libraries]),
              [AS_CASE("${enableval}",
                       [yes], [enablestrictlgpl=yes],
                       [no],  [enablestrictlgpl=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-strict-lgpl)])],
              [enablestrictlgpl=yes])


# --------------------------------------------------------------
# Allow for disable-nested
# --------------------------------------------------------------
AC_ARG_ENABLE(nested,
              AS_HELP_STRING([--disable-nested],
                             [Do not use nested autoconf subpackages]),
              [AS_CASE("${enableval}",
                       [yes], [enablenested=yes],
                       [no],  [enablenested=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-nested)])],
              [enablenested=$enableoptional])


# --------------------------------------------------------------
# XDR binary IO support - enabled by default
# This used to be tested in libmesh_core_features.m4 since your
# system either had it or it didn't. Now it's possible for the
# XDR headers to be in different places, so it's more convenient
# to test for it here.
# --------------------------------------------------------------
AC_ARG_ENABLE(xdr,
              AS_HELP_STRING([--disable-xdr],
                             [build without XDR platform-independent binary I/O]),
              enablexdr=$enableval,
              enablexdr=yes)

AS_IF([test "$enablexdr" != no],
      [
      dnl Check whether the system/compiler has xdr.h in /usr/include and glibc.
      AC_MSG_CHECKING([for built-in XDR support])
      CONFIGURE_XDR

      dnl Check for headers in /usr/include/tirpc. (Fedora 28 does this.)
      AS_IF([test "x$enablexdr" = "xno"],
            [
              AC_MSG_CHECKING([for XDR support in /usr/include/tirpc])
              old_CPPFLAGS="$CPPFLAGS"
              old_LIBS="$LIBS"
              CPPFLAGS="$CPPFLAGS -I/usr/include/tirpc"
              LIBS="$LIBS -ltirpc"

              CONFIGURE_XDR

              dnl If that worked, append the required include paths and libraries as necessary.
              AS_IF([test "x$enablexdr" = "xyes"],
                    [
                      libmesh_optional_INCLUDES="$libmesh_optional_INCLUDES -I/usr/include/tirpc"
                      libmesh_optional_LIBS="$libmesh_optional_LIBS -ltirpc"
                    ])

              dnl Reset flags after testing
              CPPFLAGS="$old_CPPFLAGS"
              LIBS="$old_LIBS"
           ])
      ])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Boost -- enabled by default
# -------------------------------------------------------------
CONFIGURE_BOOST
LIBMESH_TEST_BOOST_MOVELIB_UNIQUE_PTR
AC_CONFIG_FILES([contrib/boost/include/Makefile])
# --------------------------------------------------------------

# -------------------------------------------------------------
# MPI -- enabled by default
# -------------------------------------------------------------
AS_IF([test "x$enablempi" = xyes],
      [
        ACSM_MPI
        AS_IF([test "x$enablempi" = xyes],
              [
                AS_IF([test x"$MPI_INCLUDES" = x],,[libmesh_optional_INCLUDES="$MPI_INCLUDES $libmesh_optional_INCLUDES"])
                AS_IF([test x"$MPI_LIBS" != x], [libmesh_optional_LIBS="$MPI_LIBS $libmesh_optional_LIBS"])
                AS_IF([test x"$MPI_LDFLAGS" != x], [libmesh_optional_LIBS="$MPI_LDFLAGS $libmesh_optional_LIBS"])
              ])
      ])


# -------------------------------------------------------------------
# Petsc -- We arelady called CONFIGURE_PETSC in LIBMESH_SET_COMPILERS
# -------------------------------------------------------------------
AS_IF([test $enablepetsc != no],
      [
        CONFIGURE_PETSC
        libmesh_optional_INCLUDES="$PETSCINCLUDEDIRS $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$PETSCLINKLIBS $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_PETSC, test x$enablepetsc = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# Check for inconsistencies between PETSc and libmesh's scalar
# and index datatypes.
# -------------------------------------------------------------
AS_IF([test $enablepetsc != no],
      [
        petsc_use_64bit_indices=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_USE_64BIT_INDICES`

        dnl If PETSc is using 64-bit indices, make sure that
        dnl $dof_bytes==8, or *make* $dof_bytes=8 if the user didn't
        dnl specify otherwise, or else give an informative error.
        AS_IF([test $petsc_use_64bit_indices -gt 0 && test "$dof_bytes_setting" != "explicit"],
              [AS_IF([test $dof_bytes != "8"],
                     [AC_MSG_RESULT([>>> adopting PETSc dof_id size: 8])])
               dof_bytes=8
               dof_bytes_setting="explicit"
               AC_DEFINE(DOF_ID_BYTES, 8, [size of dof_id])
              ])
        AS_IF([test $petsc_use_64bit_indices -gt 0 && test "$dof_bytes" != "8"],
              [AC_MSG_ERROR([<<< PETSc is using 64-bit indices, you must configure libmesh with --with-dof-id-bytes=8. >>>])])

        dnl If PETSc is using 32-bit indices, make sure that
        dnl libmesh's $dof_bytes<=4, or *make* $dof_bytes=4 if the
        dnl user didn't specify otherwise, or else give an informative
        dnl error.
        AS_IF([test $petsc_use_64bit_indices = "0" && test "$dof_bytes_setting" != "explicit"],
              [AS_IF([test $dof_bytes != "4"],
                     [AC_MSG_RESULT([>>> adopting PETSc dof_id size: 4])])
               dof_bytes=4
               dof_bytes_setting="explicit"
              ])
        AS_IF([test "$petsc_use_64bit_indices" = "0" && test $dof_bytes -gt 4],
              [AC_MSG_ERROR([<<< PETSc is using 32-bit indices, you must configure libmesh with --with-dof-id-bytes=<1|2|4>. >>>])])

        dnl Libmesh must use {complex,real} scalars when PETSc uses {complex,real} scalars.
        petsc_use_complex=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_USE_COMPLEX`

        AS_IF([test $petsc_use_complex -gt 0 && test "$enablecomplex" = "no"],
              [AC_MSG_ERROR([<<< PETSc was built with complex scalars, you must configure libmesh with --enable-complex. >>>])])

        AS_IF([test "$petsc_use_complex" = "0" && test "$enablecomplex" = "yes"],
              [AC_MSG_ERROR([<<< PETSc was built with real scalars, you must configure libmesh with --disable-complex. >>>])])
      ])



# -------------------------------------------------------------
# SLEPc -- enabled by default
# -------------------------------------------------------------
CONFIGURE_SLEPC
AS_IF([test $enableslepc = yes],
      [
        libmesh_optional_INCLUDES="$SLEPC_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$SLEPC_LIBS $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_SLEPC, test x$enableslepc = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# "Trilinos" -- enabled by default unless we're building with
#               complex numbers.
# -------------------------------------------------------------
CONFIGURE_TRILINOS
AS_IF([test "$enabletrilinos" = yes],
      [
        libmesh_optional_INCLUDES="$TRILINOS_INCLUDES $AZTECOO_INCLUDES $NOX_INCLUDES $ML_INCLUDES $TPETRA_INCLUDES $DTK_INCLUDES $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$TRILINOS_LIBS $AZTECOO_LIBS $NOX_LIBS $ML_LIBS $TPETRA_INCLUDES $DTK_INCLUDES $libmesh_optional_LIBS"
      ])
# -------------------------------------------------------------


# -------------------------------------------------------------
# Choose between TBB, OpenMP, and pthreads thread models.
# The user can control this by configuring with
#
# --with-thread-model={tbb,pthread,auto,none}
#
# where "auto" will try to automatically detect the best possible
# version (see threads.m4).
# -------------------------------------------------------------
ACX_BEST_THREAD



# -------------------------------------------------------------
# LASPACK iterative solvers -- enabled unless
# --enable-strict-lgpl is specified
# -------------------------------------------------------------
AS_IF([test $enablestrictlgpl = yes],
      [
        AC_MSG_RESULT([<<< Laspack support is disabled, configure with --disable-strict-lgpl to enable it >>>])
        enablelaspack=no;
      ],
      [
        CONFIGURE_LASPACK
        AS_IF([test $enablelaspack = yes],
              [libmesh_contrib_INCLUDES="$LASPACK_INCLUDE $libmesh_contrib_INCLUDES"])
      ])

AM_CONDITIONAL(LIBMESH_ENABLE_LASPACK, test x$enablelaspack = xyes)
AC_CONFIG_FILES([contrib/laspack/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Space filling curves -- enabled unless
# --enable-strict-lgpl is specified
# -------------------------------------------------------------
AS_IF([test $enablestrictlgpl = yes],
      [
        AC_MSG_RESULT([<<< The space filling curves partitioner is disabled, configure with --disable-strict-lgpl to enable it >>>])
        enablesfc=no;
      ],
      [
        CONFIGURE_SFC
        AS_IF([test $enablesfc = yes],
              [libmesh_contrib_INCLUDES="$SFC_INCLUDE $libmesh_contrib_INCLUDES"])
      ])

AM_CONDITIONAL(LIBMESH_ENABLE_SFC, test x$enablesfc = xyes)
AC_CONFIG_FILES([contrib/sfcurves/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Compressed Streams with gzstream -- enabled by default
# -------------------------------------------------------------
CONFIGURE_GZ
AS_IF([test "$enablegz" = yes],
      [
        libmesh_contrib_INCLUDES="$GZSTREAM_INCLUDE $libmesh_contrib_INCLUDES"
        libmesh_optional_LIBS="-lz $libmesh_optional_LIBS"
      ])

AM_CONDITIONAL(LIBMESH_ENABLE_GZSTREAMS, test x$enablegz = xyes)
AC_CONFIG_FILES([contrib/gzstream/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Compressed Files with bzip2
# -------------------------------------------------------------
AC_ARG_ENABLE(bzip2,
              AS_HELP_STRING([--disable-bzip2],
                             [build without bzip2 compressed I/O support]),
              enablebz2=$enableval,
              enablebz2=$enableoptional)

AS_IF([test "$enablebz2" != no],
      [
        AC_CHECK_PROG(BZIP2,bzip2,bzip2,none,$PATH)
        AS_IF([test "$BZIP2" = bzip2],
              [
                AC_CHECK_PROG(BUNZIP2,bunzip2,bunzip2,none,$PATH)
                AS_IF([test "$BUNZIP2" = bunzip2],
                      [
                        AC_MSG_RESULT(<<< Using bzip2/bunzip2 for writing/reading compressed .bz2 files >>>)
                        AC_DEFINE(HAVE_BZIP, 1, [Flag indicating bzip2/bunzip2 are available for handling compressed .bz2 files])
                      ])
              ])
      ])
# -------------------------------------------------------------


# -------------------------------------------------------------
# Compressed Files with xz
# -------------------------------------------------------------
AC_ARG_ENABLE(xz,
              AS_HELP_STRING([--disable-xz],
                             [build without xz compressed I/O support]),
              enablexz=$enableval,
              enablexz=$enableoptional)

AS_IF([test "$enablexz" != no],
      [
        AC_CHECK_PROG(XZ,xz,xz,none,$PATH)
        AS_IF([test "$XZ" = xz],
              [
                AC_MSG_RESULT(<<< Using xz for writing/reading compressed .xz files >>>)
                AC_DEFINE(HAVE_XZ, 1, [Flag indicating xz is available for handling compressed .xz files])
              ])
      ])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Tecplot, from source -- enabled by default
# -------------------------------------------------------------
CONFIGURE_TECIO
AS_IF([test $enabletecio = yes],
      [libmesh_contrib_INCLUDES="$TECIO_INCLUDE $libmesh_contrib_INCLUDES"])
AM_CONDITIONAL(LIBMESH_ENABLE_TECIO, test x$enabletecio = xyes)
AC_CONFIG_FILES([contrib/tecplot/tecio/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Tecplot, vendor provided libraries -- disabled by default
# -------------------------------------------------------------
CONFIGURE_TECPLOT
AS_IF([test $enabletecplot = yes],
      [
        libmesh_contrib_INCLUDES="$TECPLOT_INCLUDE $libmesh_contrib_INCLUDES"
        libmesh_installed_LIBS="$libmesh_installed_LIBS -ltecio_vendor"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_TECPLOT, test x$enabletecplot = xyes)
AC_CONFIG_FILES([contrib/tecplot/binary/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Metis Partitioning -- enabled by default
# -------------------------------------------------------------
CONFIGURE_METIS
AS_IF([test $enablemetis = yes],
      [libmesh_contrib_INCLUDES="$METIS_INCLUDE $libmesh_contrib_INCLUDES"
       libmesh_optional_LIBS="$METIS_LIB $libmesh_optional_LIBS"])
AM_CONDITIONAL(LIBMESH_ENABLE_METIS, test x$enablemetis = xyes)
AC_CONFIG_FILES([contrib/metis/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Parmetis Partitioning -- enabled by default
# -------------------------------------------------------------
AS_IF([test $enablemetis = yes],
      [CONFIGURE_PARMETIS],
      [enableparmetis=no])
AS_IF([test $enableparmetis = yes],
      [libmesh_contrib_INCLUDES="$PARMETIS_INCLUDE $libmesh_contrib_INCLUDES"
       libmesh_optional_LIBS="$PARMETIS_LIB $libmesh_optional_LIBS"])
AM_CONDITIONAL(LIBMESH_ENABLE_PARMETIS, test x$enableparmetis = xyes)
AC_CONFIG_FILES([contrib/parmetis/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Doxygen - look for doxygen (a documentation tool)
# -------------------------------------------------------------
AC_PATH_PROG(DOXYGEN, doxygen)
AC_SUBST(DOXYGEN)
AS_IF([test "x$DOXYGEN" != x],
      [
        dnl -----------------------------------------------------------
        dnl Dot -- lets doxygen generate pretty class diagrams
        dnl -----------------------------------------------------------
        AC_PATH_PROG(DOT, dot)
        HAVE_DOT=NO
        AS_IF([test "x$DOT" != x],
              [
                HAVE_DOT=YES
                DOTPATH=$PWD/doc
                AC_SUBST(DOT)
                AC_SUBST(DOTPATH)
              ])
        AC_SUBST(HAVE_DOT)
      ])
# -------------------------------------------------------------



# -------------------------------------------------------------
# TetGen -- enabled unless --enable-strict-lgpl is specified
# -------------------------------------------------------------
AS_IF([test $enablestrictlgpl = yes],
      [
        AC_MSG_RESULT([<<< Tetgen support is disabled, configure with --disable-strict-lgpl to enable it >>>])
        enabletetgen=no;
      ],
      [
        CONFIGURE_TETGEN
        AS_IF([test $enabletetgen = yes],
              [libmesh_contrib_INCLUDES="$TETGEN_INCLUDE $libmesh_contrib_INCLUDES"])
      ])

AM_CONDITIONAL(LIBMESH_ENABLE_TETGEN, test x$enabletetgen = xyes)
AC_CONFIG_FILES([contrib/tetgen/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Triangle -- enabled unless --enable-strict-lgpl is specified
# -------------------------------------------------------------
AS_IF([test $enablestrictlgpl = yes],
      [
        AC_MSG_RESULT([<<< Triangle meshing support is disabled, configure with --disable-strict-lgpl to enable it >>>])
        enabletriangle=no;
      ],
      [
        CONFIGURE_TRIANGLE
        AS_IF([test $enabletriangle = yes],
              [libmesh_contrib_INCLUDES="$TRIANGLE_INCLUDE $libmesh_contrib_INCLUDES"])
      ])

AM_CONDITIONAL(LIBMESH_ENABLE_TRIANGLE, test x$enabletriangle = xyes)
AC_CONFIG_FILES([contrib/triangle/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Qhull -- enabled by default
# -------------------------------------------------------------
CONFIGURE_QHULL
AS_IF([test $enableqhull = yes],
      [libmesh_contrib_INCLUDES="$QHULL_INCLUDE $libmesh_contrib_INCLUDES"])
AM_CONDITIONAL(LIBMESH_ENABLE_QHULL, test x$enableqhull = xyes)
AC_CONFIG_FILES([contrib/qhull/qhull/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# GMV -- file I/O API is enabled by default (it is distributed in contrib)
# -------------------------------------------------------------
CONFIGURE_GMV
AS_IF([test x$enablegmv = xyes],
      [libmesh_contrib_INCLUDES="$GMV_INCLUDE $libmesh_contrib_INCLUDES"])
AM_CONDITIONAL(LIBMESH_ENABLE_GMV, test x$enablegmv = xyes)
AC_CONFIG_FILES([contrib/gmv/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# VTK -- Mesh I/O API is enabled by default
# -------------------------------------------------------------
CONFIGURE_VTK
AS_IF([test x$enablevtk = xyes],
      [
        libmesh_optional_INCLUDES="$VTK_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$VTK_LIBRARY $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_VTK, test x$enablevtk = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# Eigen -- Optimized linear algebra routines, enabled by default
# -------------------------------------------------------------

# we require Eigen/Sparse support if we're going to enable Eigen
enableeigensparse=yes
# we test with Eigen 3.1.2, so if the user has their own Eigen it
# should be at least that new.
CONFIGURE_EIGEN(3.1.2,no)
AS_IF([test "x$enableeigen" = xyes],
      [
        dnl if we are installing our own Eigen, add it to the contrib search path
        dnl which is not exported during install
        AS_IF([test "x$install_internal_eigen" = xyes],
              [libmesh_contrib_INCLUDES="$EIGEN_INCLUDE $libmesh_contrib_INCLUDES"],
              dnl if we are depending on an external Eigen, add it to the optional search
              dnl path, which gets exported during install
              [libmesh_optional_INCLUDES="$EIGEN_INCLUDE $libmesh_optional_INCLUDES"])
      ])
AC_CONFIG_FILES([contrib/eigen/eigen/Makefile])
AM_CONDITIONAL(LIBMESH_ENABLE_EIGEN, test x$enableeigen = xyes)
AM_CONDITIONAL(LIBMESH_INSTALL_INTERNAL_EIGEN, test x$install_internal_eigen = xyes)
#--------------------------------------------------------------



# -------------------------------------------------------------
# GLPK -- Needed by the SCM routine of rbOOmit.
# Enabled by default.
# -------------------------------------------------------------
CONFIGURE_GLPK
AS_IF([test x$enableglpk = xyes],
      [
        libmesh_optional_INCLUDES="$GLPK_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$GLPK_LIBRARY $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_GLPK, test x$enableglpk = xyes)
# -------------------------------------------------------------


# -------------------------------------------------------------
# NLOPT -- A library of nonlinear optimization routines.
# -------------------------------------------------------------
CONFIGURE_NLOPT
AS_IF([test x$enablenlopt = xyes],
      [
        libmesh_optional_INCLUDES="$NLOPT_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$NLOPT_LIBRARY $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_NLOPT, test x$enablenlopt = xyes)
# -------------------------------------------------------------


# -------------------------------------------------------------
# CAPNPROTO -- Serialization library.
# -------------------------------------------------------------
CONFIGURE_CAPNPROTO
AS_IF([test x$enablecapnproto = xyes],
      [
        libmesh_optional_INCLUDES="$CAPNPROTO_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$CAPNPROTO_LIBRARY $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_CAPNPROTO, test x$enablecapnproto = xyes)
AC_CONFIG_FILES([contrib/capnproto/Makefile])
# -------------------------------------------------------------

# -------------------------------------------------------------
# libcurl -- enabled by default
# Note: I tried to use the m4 files ax_lib_curl.m4 and
# ax_path_generic.m4 from the autoconf-archive for this, but they
# would not work (bootstrap failed!) on either Linux or OSX.
# -------------------------------------------------------------
CONFIGURE_CURL
AS_IF([test x$enablecurl = xyes],
      [
        libmesh_optional_INCLUDES="$CURL_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$CURL_LIBRARY $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_CURL, test x$enablecurl = xyes)
# -------------------------------------------------------------



# --------------------------------------------------------------
# HDF5 -- enabled by default
# --------------------------------------------------------------
CONFIGURE_HDF5
AS_IF([test $enablehdf5 = yes],
      [
        libmesh_optional_INCLUDES="$HDF5_CPPFLAGS $libmesh_optional_INCLUDES"

        dnl If the HDF5 C++ interface was found, add the C++ library to the link line.
        AS_IF([test "$hdf5_has_cxx" = yes],
              [libmesh_optional_LIBS="$HDF5_CXXLIBS $libmesh_optional_LIBS"])

        dnl And add the HDF5 C library to the link line.
        libmesh_optional_LIBS="$HDF5_LIBS $libmesh_optional_LIBS"
      ])
AM_CONDITIONAL(LIBMESH_ENABLE_HDF5, test x$enablehdf5 = xyes)


# --------------------------------------------------------------
# netCDF -- enabled by default (it is distributed in contrib)
# --------------------------------------------------------------
CONFIGURE_NETCDF
AS_IF([test $enablenetcdf = yes],
      [libmesh_contrib_INCLUDES="$NETCDF_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_NETCDF,    test x$enablenetcdf  = xyes)
AM_CONDITIONAL(LIBMESH_ENABLE_NETCDF_V4, test x$netcdfversion = x4)

# -------------------------------------------------------------
# ExodusII -- enabled by default (it is distributed in contrib)
# (note that ExodusII requires netCDF
# -------------------------------------------------------------
CONFIGURE_EXODUS
AS_IF([test $enableexodus = yes],
      [libmesh_contrib_INCLUDES="$EXODUS_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS,      test x$enableexodus  = xyes)
AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS_V509, test x$exodusversion = xv5.09)
AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS_V522, test x$exodusversion = xv5.22)

# -------------------------------------------------------------
# Nemesis -- enabled by default (it is distributed in contrib)
# (note that Nemesis requires netCDF and exodus)
# -------------------------------------------------------------
CONFIGURE_NEMESIS
AS_IF([test $enablenemesis = yes],
      [libmesh_contrib_INCLUDES="$NEMESIS_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS,      test x$enablenemesis  = xyes)
AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS_V309, test x$nemesisversion = xv3.09)
AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS_V522, test x$nemesisversion = xv5.22)



# -------------------------------------------------------------
# libHilbert -- distributed in ./contrib,
#               enabled by default
# -------------------------------------------------------------
CONFIGURE_LIBHILBERT
AS_IF([test $enablelibhilbert = yes],
      [libmesh_contrib_INCLUDES="$LIBHILBERT_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_LIBHILBERT, test x$enablelibhilbert = xyes)
AC_CONFIG_FILES([contrib/libHilbert/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# fparser -- distributed in ./contrib,
#            enabled by default
# -------------------------------------------------------------
CONFIGURE_FPARSER
AS_IF([test $enablefparser = yes],
      [libmesh_contrib_INCLUDES="$FPARSER_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_FPARSER, test x$enablefparser = xyes)
AC_CONFIG_FILES([contrib/fparser/Makefile])
AC_CONFIG_FILES([contrib/fparser/extrasrc/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# cppunit C++ unit testing -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(cppunit,
              AS_HELP_STRING([--disable-cppunit],
                             [Build without cppunit C++ unit testing support]),
              [AS_CASE("${enableval}",
                       [yes], [enablecppunit=yes],
                       [no],  [enablecppunit=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-cppunit)])],
              [enablecppunit=yes])
AS_IF([test "$enablecppunit" = yes],
      [AM_PATH_CPPUNIT])

AM_CONDITIONAL(LIBMESH_ENABLE_CPPUNIT, test x$enablecppunit = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# nanoflann -- enabled by default
# -------------------------------------------------------------
CONFIGURE_NANOFLANN
AS_IF([test $enablenanoflann = yes],
      [libmesh_contrib_INCLUDES="$NANOFLANN_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_NANOFLANN, test x$enablenanoflann = xyes)
AC_CONFIG_FILES([contrib/nanoflann/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# MetaPhysicL -- enabled by default
# -------------------------------------------------------------
CONFIGURE_METAPHYSICL
AS_IF([test $enablemetaphysicl = yes],
      [libmesh_contrib_INCLUDES="$METAPHYSICL_INCLUDE $libmesh_contrib_INCLUDES"])

AM_CONDITIONAL(LIBMESH_ENABLE_METAPHYSICL, test x$enablemetaphysicl = xyes)
# -------------------------------------------------------------




AS_IF([test "$enableoptional" != no],
      [
        AC_MSG_RESULT(----------------------------------------------)
        AC_MSG_RESULT(--- Done configuring for optional packages ---)
        AC_MSG_RESULT(----------------------------------------------)
      ])

# clean up values, if we have perl.  This step is purely cosmetic, but
# helps create readable (and easier to debug) compile and link lines
# by stripping out repeated entries.  This can happen for example when
# several optional packages all want to include and link agains the
# same MPI.
AS_IF([test -x $PERL],
      [
        AS_IF([test -f $srcdir/contrib/bin/strip_dup_incl_paths.pl],
              [
                AC_MSG_RESULT(removing duplicate include paths...)
                libmesh_optional_INCLUDES=`$PERL $srcdir/contrib/bin/strip_dup_incl_paths.pl $libmesh_optional_INCLUDES`
              ])

        AS_IF([test -f $srcdir/contrib/bin/strip_dup_libs.pl],
              [
                AC_MSG_RESULT(removing duplicate libraries...)
                libmesh_optional_LIBS=`$PERL $srcdir/contrib/bin/strip_dup_libs.pl $libmesh_optional_LIBS`
              ])
      ])

# substitute values
AC_SUBST(libmesh_optional_INCLUDES)
AC_SUBST(libmesh_optional_LIBS)
AC_SUBST(libmesh_contrib_INCLUDES)
AC_SUBST(libmesh_pkgconfig_requires)
AC_SUBST(libmesh_installed_LIBS)
])
