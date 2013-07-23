# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_CONFIGURE_OPTIONAL_PACKAGES],
[


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
# Allow for disable-optional
# --------------------------------------------------------------
AC_ARG_ENABLE(optional,
              AC_HELP_STRING([--enable-optional],
                             [en/disable optional external libraries]),
	      [case "${enableval}" in
	      	  yes) enableoptional=yes ;;
		   no) enableoptional=no ;;
 		    *) AC_MSG_ERROR(bad value ${enableval} for --enable-optional) ;;
	       esac],
              [enableoptional=yes])

# Note that even when optional packages are disabled we need to
# run their m4 macros to get proper AM_CONDITIONALs.  Just be
# quiet about it...
if test "$enableoptional" != no ; then
   AC_MSG_RESULT(---------------------------------------------)
   AC_MSG_RESULT(----- Configuring for optional packages -----)
   AC_MSG_RESULT(---------------------------------------------)
fi


# --------------------------------------------------------------
# Allow for disable-nested
# --------------------------------------------------------------
AC_ARG_ENABLE(nested,
              AC_HELP_STRING([--enable-nested],
                             [en/disable nested autoconf subpackages]),
	      [case "${enableval}" in
	      	  yes) enablenested=yes ;;
		   no) enablenested=no ;;
 		    *) AC_MSG_ERROR(bad value ${enableval} for --enable-nested) ;;
	       esac],
              [enablenested=$enableoptional])

# -------------------------------------------------------------
# Boost -- enabled by default
# -------------------------------------------------------------
CONFIGURE_BOOST
AC_CONFIG_FILES([contrib/boost/include/Makefile])
# --------------------------------------------------------------



# -------------------------------------------------------------
# Petsc -- enabled by default
# -------------------------------------------------------------
CONFIGURE_PETSC
if (test $enablempi != no) ; then
  libmesh_optional_INCLUDES="$MPI_INCLUDES_PATHS $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$MPI_LIBS_PATHS $MPI_LIBS $libmesh_optional_LIBS"
fi
if (test $enablepetsc != no) ; then
  libmesh_optional_INCLUDES="$PETSCINCLUDEDIRS $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$PETSCLINKLIBS $libmesh_optional_LIBS"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_PETSC, test x$enablepetsc = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# SLEPc -- enabled by default
# -------------------------------------------------------------
CONFIGURE_SLEPC
if (test $enableslepc = yes ) ; then
  libmesh_optional_INCLUDES="$SLEPC_INCLUDE $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$SLEPC_LIBS $libmesh_optional_LIBS"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_SLEPC, test x$enableslepc = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# "Trilinos" -- enabled by default unless we're building with
#               complex numbers.
# -------------------------------------------------------------
CONFIGURE_TRILINOS
if ( test "$enabletrilinos" = yes ); then
  libmesh_optional_INCLUDES="$TRILINOS_INCLUDES $AZTECOO_INCLUDES $NOX_INCLUDES $ML_INCLUDES $TPETRA_INCLUDES $DTK_INCLUDES $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$TRILINOS_LIBS $AZTECOO_LIBS $NOX_LIBS $ML_LIBS $TPETRA_INCLUDES $DTK_INCLUDES $libmesh_optional_LIBS"
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Intel's Threading Building Blocks -- enabled by default
# -------------------------------------------------------------
CONFIGURE_TBB
if (test $enabletbb = yes); then
  libmesh_optional_INCLUDES="$TBB_INCLUDE $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$TBB_LIBRARY $libmesh_optional_LIBS"
fi
# -------------------------------------------------------------

# -------------------------------------------------------------
# Pthread support -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(pthreads,
              AC_HELP_STRING([--enable-pthreads],
                             [build with threading support via POSIX threads (pthreads)]),
		[case "${enableval}" in
		  yes)  enablepthreads=yes ;;
		   no)  enablepthreads=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-pthreads) ;;
		 esac],
		 [enablepthreads=$enableoptional])

if (test "$enablepthreads" != no) ; then
  AX_PTHREAD
fi

if (test x$ax_pthread_ok = xyes); then
  AC_MSG_RESULT(<<< Configuring library with pthread support >>>)
  libmesh_optional_INCLUDES="$PTHREAD_CFLAGS $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$PTHREAD_LIBS $libmesh_optional_LIBS"
else
  enablepthreads=no
fi
# -------------------------------------------------------------


# -------------------------------------------------------------
# C++ Thread Support  -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(cppthreads,
             AC_HELP_STRING([--enable-cppthreads],
                            [Build with C++ std::thread support]),
             enablecppthreads=$enableval,
             enablecppthreads=yes)
if (test "$enablecppthreads" != no) ; then
  ACX_BEST_THREAD
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# LASPACK iterative solvers -- enabled by default
# -------------------------------------------------------------
CONFIGURE_LASPACK
if (test $enablelaspack = yes); then
  libmesh_contrib_INCLUDES="$LASPACK_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_LASPACK, test x$enablelaspack = xyes)
AC_CONFIG_FILES([contrib/laspack/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Space filling curves -- enabled by default
# -------------------------------------------------------------
CONFIGURE_SFC
if (test $enablesfc = yes); then
  libmesh_contrib_INCLUDES="$SFC_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_SFC, test x$enablesfc = xyes)
AC_CONFIG_FILES([contrib/sfcurves/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Compressed Streams with gzstream -- enabled by default
# -------------------------------------------------------------
CONFIGURE_GZ
if (test "$enablegz" = yes) ; then
  libmesh_contrib_INCLUDES="$GZSTREAM_INCLUDE $libmesh_contrib_INCLUDES"
  libmesh_optional_LIBS="-lz $libmesh_optional_LIBS"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_GZSTREAMS, test x$enablegz = xyes)
AC_CONFIG_FILES([contrib/gzstream/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Compressed Files with bzip2
# -------------------------------------------------------------
AC_ARG_ENABLE(bzip2,
              AC_HELP_STRING([--enable-bzip2],
                             [build with bzip2 compressed I/O suppport]),
              enablebz2=$enableval,
              enablebz2=$enableoptional)

if (test "$enablebz2" != no) ; then
   #           Var   | look for | name if found | name if not | where
   AC_CHECK_PROG(BZIP2,  bzip2,      bzip2,           none,      $PATH)
   if test "$BZIP2" = bzip2; then
       AC_CHECK_PROG(BUNZIP2,  bunzip2,     bunzip2,        none,      $PATH)
       if test "$BUNZIP2" = bunzip2; then
         AC_MSG_RESULT(<<< Using bzip2/bunzip2 for writing/reading compressed .bz2 files >>>)
         AC_DEFINE(HAVE_BZIP, 1,
                   [Flag indicating bzip2/bunzip2 are available for handling compressed .bz2 files])
       fi
   fi
fi
# -------------------------------------------------------------


# -------------------------------------------------------------
# Compressed Files with xz
# -------------------------------------------------------------
AC_ARG_ENABLE(xz,
              AC_HELP_STRING([--enable-xz],
                             [build with xz compressed I/O suppport]),
              enablexz=$enableval,
              enablexz=$enableoptional)

if (test "$enablexz" != no) ; then
   #           Var   | look for | name if found | name if not | where
   AC_CHECK_PROG(XZ,  xz,      xz,           none,      $PATH)
   if test "$XZ" = xz; then
      AC_MSG_RESULT(<<< Using xz for writing/reading compressed .xz files >>>)
      AC_DEFINE(HAVE_XZ, 1,
                [Flag indicating xz is available for handling compressed .xz files])
   fi
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Tecplot, from source -- enabled by default
# -------------------------------------------------------------
CONFIGURE_TECIO
if (test $enabletecio = yes); then
  libmesh_contrib_INCLUDES="$TECIO_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_TECIO, test x$enabletecio = xyes)
AC_CONFIG_FILES([contrib/tecplot/tecio/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Tecplot, vendor provided libraries -- disabled by default
# -------------------------------------------------------------
CONFIGURE_TECPLOT
if (test $enabletecplot = yes); then
  libmesh_contrib_INCLUDES="$TECPLOT_INCLUDE $libmesh_contrib_INCLUDES"
  libmesh_installed_LIBS="$libmesh_installed_LIBS -ltecio_vendor"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_TECPLOT, test x$enabletecplot = xyes)
AC_CONFIG_FILES([contrib/tecplot/binary/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Metis Partitioning -- enabled by default
# -------------------------------------------------------------
CONFIGURE_METIS
if (test $enablemetis = yes); then
  libmesh_contrib_INCLUDES="$METIS_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_METIS, test x$enablemetis = xyes)
AC_CONFIG_FILES([contrib/metis/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Parmetis Partitioning -- enabled by default
# -------------------------------------------------------------
CONFIGURE_PARMETIS
if (test $enableparmetis = yes); then
  libmesh_contrib_INCLUDES="$PARMETIS_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_PARMETIS, test x$enableparmetis = xyes)
AC_CONFIG_FILES([contrib/parmetis/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Doxygen - look for doxygen (a documentation tool)
# -------------------------------------------------------------
AC_PATH_PROG(DOXYGEN, doxygen)
AC_SUBST(DOXYGEN)
if test "x$DOXYGEN" != x ; then
  # -----------------------------------------------------------
  # Dot -- lets doxygen generate pretty class diagrams
  # -----------------------------------------------------------
  AC_PATH_PROG(DOT, dot)
  HAVE_DOT=NO
  if test "x$DOT" != x ; then
    HAVE_DOT=YES
    DOTPATH=$PWD/doc
    AC_SUBST(DOT)
    AC_SUBST(DOTPATH)
  fi
  AC_SUBST(HAVE_DOT)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# TetGen -- enabled by default
# -------------------------------------------------------------
CONFIGURE_TETGEN
if (test $enabletetgen = yes); then
  libmesh_contrib_INCLUDES="$TETGEN_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_TETGEN, test x$enabletetgen = xyes)
AC_CONFIG_FILES([contrib/tetgen/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# Triangle -- enabled by default (it is distributed in contrib)
# -------------------------------------------------------------
CONFIGURE_TRIANGLE
if (test $enabletriangle = yes); then
  libmesh_contrib_INCLUDES="$TRIANGLE_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_TRIANGLE, test x$enabletriangle = xyes)
AC_CONFIG_FILES([contrib/triangle/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# GMV -- file I/O API is enabled by default (it is distributed in contrib)
# -------------------------------------------------------------
CONFIGURE_GMV
if (test x$enablegmv = xyes); then
  libmesh_contrib_INCLUDES="$GMV_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_GMV, test x$enablegmv = xyes)
AC_CONFIG_FILES([contrib/gmv/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# VTK -- Mesh I/O API is enabled by default
# -------------------------------------------------------------
CONFIGURE_VTK
if (test x$enablevtk = xyes); then
  libmesh_optional_INCLUDES="$VTK_INCLUDE $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$VTK_LIBRARY $libmesh_optional_LIBS"
fi
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
if (test x$enableeigen = xyes); then
  # if we are installing our own Eigen, add it to the contrib search path
  # which is not exported during install
  if (test x$install_internal_eigen = xyes); then
    libmesh_contrib_INCLUDES="$EIGEN_INCLUDE $libmesh_contrib_INCLUDES"
  # if we are depending on an external Eigen, add it to the optional search
  # path, which gets exported during install
  else
    libmesh_optional_INCLUDES="$EIGEN_INCLUDE $libmesh_optional_INCLUDES"
  fi
fi
AC_CONFIG_FILES([contrib/eigen/eigen/Makefile])
AM_CONDITIONAL(LIBMESH_ENABLE_EIGEN, test x$enableeigen = xyes)
AM_CONDITIONAL(LIBMESH_INSTALL_INTERNAL_EIGEN, test x$install_internal_eigen = xyes)
#--------------------------------------------------------------



# -------------------------------------------------------------
# GLPK -- Needed by the SCM routine of rbOOmit.
# Enabled by default.
# -------------------------------------------------------------
CONFIGURE_GLPK
if (test x$enableglpk = xyes); then
  libmesh_optional_INCLUDES="$GLPK_INCLUDE $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$GLPK_LIBRARY $libmesh_optional_LIBS"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_GLPK, test x$enableglpk = xyes)
# -------------------------------------------------------------



# --------------------------------------------------------------
# HDF5 -- enabled by default
# --------------------------------------------------------------
CONFIGURE_HDF5
if (test $enablehdf5 = yes); then
  libmesh_optional_INCLUDES="$HDF5_CPPFLAGS $libmesh_optional_INCLUDES"
  libmesh_optional_LIBS="$HDF5_LIBS $libmesh_optional_LIBS"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_HDF5, test x$enablehdf5 = xyes)


# --------------------------------------------------------------
# netCDF -- enabled by default (it is distributed in contrib)
# --------------------------------------------------------------
CONFIGURE_NETCDF
if (test $enablenetcdf = yes); then
  libmesh_contrib_INCLUDES="$NETCDF_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_NETCDF,    test x$enablenetcdf  = xyes)
AM_CONDITIONAL(LIBMESH_ENABLE_NETCDF_V3, test x$netcdfversion = x3)
AM_CONDITIONAL(LIBMESH_ENABLE_NETCDF_V4, test x$netcdfversion = x4)

   # -------------------------------------------------------------
   # ExodusII -- enabled by default (it is distributed in contrib)
   # (note that ExodusII requires netCDF
   # -------------------------------------------------------------
   CONFIGURE_EXODUS
   if (test $enableexodus = yes); then
     libmesh_contrib_INCLUDES="$EXODUS_INCLUDE $libmesh_contrib_INCLUDES"
   fi
   AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS,      test x$enableexodus  = xyes)
   AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS_V509, test x$exodusversion = xv5.09)
   AM_CONDITIONAL(LIBMESH_ENABLE_EXODUS_V522, test x$exodusversion = xv5.22)

      # -------------------------------------------------------------
      # Nemesis -- enabled by default (it is distributed in contrib)
      # (note that Nemesis requires netCDF and exodus)
      # -------------------------------------------------------------
      CONFIGURE_NEMESIS
      if (test $enablenemesis = yes); then
         libmesh_contrib_INCLUDES="$NEMESIS_INCLUDE $libmesh_contrib_INCLUDES"
      fi
      AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS,      test x$enablenemesis  = xyes)
      AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS_V309, test x$nemesisversion = xv3.09)
      AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS_V522, test x$nemesisversion = xv5.22)
      # -------------------------------------------------------------
   # -------------------------------------------------------------
# -------------------------------------------------------------



# -------------------------------------------------------------
# libHilbert -- distributed in ./contrib,
#               enabled by default
# -------------------------------------------------------------
CONFIGURE_LIBHILBERT
if (test $enablelibhilbert = yes); then
  libmesh_contrib_INCLUDES="$LIBHILBERT_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_LIBHILBERT, test x$enablelibhilbert = xyes)
AC_CONFIG_FILES([contrib/libHilbert/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# fparser -- distributed in ./contrib,
#            enabled by default
# -------------------------------------------------------------
CONFIGURE_FPARSER
if (test $enablefparser = yes); then
  libmesh_contrib_INCLUDES="$FPARSER_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_FPARSER, test x$enablefparser = xyes)
AC_CONFIG_FILES([contrib/fparser/Makefile])
AC_CONFIG_FILES([contrib/fparser/extrasrc/Makefile])
# -------------------------------------------------------------



# -------------------------------------------------------------
# cppunit C++ unit testing -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(cppunit,
             AC_HELP_STRING([--enable-cppunit],
                            [Build with cppunit C++ unit testing support]),
		[case "${enableval}" in
		  yes)  enablecppunit=yes ;;
		   no)  enablecppunit=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cppunit) ;;
		 esac],
		 [enablecppunit=$enableoptional])
if (test "$enablecppunit" = yes) ; then
  AM_PATH_CPPUNIT([1.10.0],[enablecppunit=yes],[enablecppunit=no])
fi
AM_CONDITIONAL(LIBMESH_ENABLE_CPPUNIT, test x$enablecppunit = xyes)
# -------------------------------------------------------------



# -------------------------------------------------------------
# nanoflann -- enabled by default
# -------------------------------------------------------------
CONFIGURE_NANOFLANN
if (test $enablenanoflann = yes); then
  libmesh_contrib_INCLUDES="$NANOFLANN_INCLUDE $libmesh_contrib_INCLUDES"
fi
AM_CONDITIONAL(LIBMESH_ENABLE_NANOFLANN, test x$enablenanoflann = xyes)
AC_CONFIG_FILES([contrib/nanoflann/Makefile])
# -------------------------------------------------------------



if test "$enableoptional" != no ; then
   AC_MSG_RESULT(----------------------------------------------)
   AC_MSG_RESULT(--- Done configuring for optional packages ---)
   AC_MSG_RESULT(----------------------------------------------)
fi

# clean up values, if we have perl.  This step is purely cosmetic, but
# helps create readable (and easier to debug) compile and link lines
# by stripping out repeated entries.  This can happen for example when
# several optional packages all want to include and link agains the
# same MPI.
if (test -x $PERL); then
  if (test -f $srcdir/contrib/bin/strip_dup_incl_paths.pl); then
     AC_MSG_RESULT(removing duplicate include paths...)
     libmesh_optional_INCLUDES=`$PERL $srcdir/contrib/bin/strip_dup_incl_paths.pl $libmesh_optional_INCLUDES`
  fi
  if (test -f $srcdir/contrib/bin/strip_dup_libs.pl); then
     AC_MSG_RESULT(removing duplicate libraries...)
     libmesh_optional_LIBS=`$PERL $srcdir/contrib/bin/strip_dup_libs.pl $libmesh_optional_LIBS`
  fi
fi
# substitute values
AC_SUBST(libmesh_optional_INCLUDES)
AC_SUBST(libmesh_optional_LIBS)
AC_SUBST(libmesh_contrib_INCLUDES)
AC_SUBST(libmesh_pkgconfig_requires)
AC_SUBST(libmesh_installed_LIBS)
])
