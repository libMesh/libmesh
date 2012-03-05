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
libmesh_optional_INCLUDES=""
libmesh_optional_LIBS=""
libmesh_contrib_INCLUDES=""



# --------------------------------------------------------------
# Allow for disable-optional
# --------------------------------------------------------------
AC_ARG_ENABLE(optional,
              AC_HELP_STRING([--enable-optional],
                             [en/disable optional external libraries]),
              enableoptional=$enableval,
              enableoptional=yes)
AC_SUBST(enableoptional)	

# Note that even when optional packages are disabled we need to
# run their m4 macros to get proper AM_CONDITIONALs.  Just be
# quiet about it...
if test "$enableoptional" != no ; then
   AC_MSG_RESULT(---------------------------------------------)
   AC_MSG_RESULT(----- Configuring for optional packages -----)
   AC_MSG_RESULT(---------------------------------------------)
fi


# -------------------------------------------------------------
# Petsc -- enabled by default
# -------------------------------------------------------------
CONFIGURE_PETSC
# -------------------------------------------------------------



# -------------------------------------------------------------
# SLEPc -- enabled by default
# -------------------------------------------------------------
CONFIGURE_SLEPC
# -------------------------------------------------------------



# -------------------------------------------------------------
# "Trilinos" -- enabled by default unless we're building with
#               complex numbers.
# -------------------------------------------------------------
CONFIGURE_TRILINOS
# -------------------------------------------------------------



# -------------------------------------------------------------
# Intel's Threading Building Blocks -- enabled by default
# -------------------------------------------------------------
CONFIGURE_TBB
# -------------------------------------------------------------


# -------------------------------------------------------------
# LASPACK iterative solvers -- enabled by default
# -------------------------------------------------------------
CONFIGURE_LASPACK
# -------------------------------------------------------------




# -------------------------------------------------------------
# Space filling curves -- enabled by default
# -------------------------------------------------------------
CONFIGURE_SFC
# -------------------------------------------------------------



# -------------------------------------------------------------
# Compressed Streams with gzstream -- enabled by default
# -------------------------------------------------------------
CONFIGURE_GZ
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
# Tecplot -- enabled by default
# -------------------------------------------------------------
CONFIGURE_TECPLOT
# -------------------------------------------------------------



# -------------------------------------------------------------
# Metis Partitioning -- enabled by default
# -------------------------------------------------------------
CONFIGURE_METIS
# -------------------------------------------------------------



# -------------------------------------------------------------
# Parmetis Partitioning -- enabled by default
# -------------------------------------------------------------
CONFIGURE_PARMETIS
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
# -------------------------------------------------------------



# -------------------------------------------------------------
# Triangle -- enabled by default (it is distributed in contrib)
# -------------------------------------------------------------
CONFIGURE_TRIANGLE
# -------------------------------------------------------------



# -------------------------------------------------------------
# GMV -- file I/O API is enabled by default (it is distributed in contrib)
# -------------------------------------------------------------
CONFIGURE_GMV
# -------------------------------------------------------------



# -------------------------------------------------------------
# VTK -- Mesh I/O API is enabled by default
# -------------------------------------------------------------
CONFIGURE_VTK
# -------------------------------------------------------------



# -------------------------------------------------------------
# Eigen -- Optimized linear algebra routines, enabled by default
# -------------------------------------------------------------
CONFIGURE_EIGEN
# -------------------------------------------------------------



# -------------------------------------------------------------
# GLPK -- Needed by the SCM routine of rbOOmit.
# Enabled by default.
# -------------------------------------------------------------
CONFIGURE_GLPK
# -------------------------------------------------------------



# --------------------------------------------------------------
# netCDF -- enabled by default (it is distributed in contrib)
# --------------------------------------------------------------
CONFIGURE_NETCDF

   # -------------------------------------------------------------
   # ExodusII -- enabled by default (it is distributed in contrib)
   # (note that ExodusII requires netCDF
   # -------------------------------------------------------------
   CONFIGURE_EXODUS

      # -------------------------------------------------------------
      # Nemesis -- enabled by default (it is distributed in contrib)
      # (note that Nemesis requires netCDF and exodus)
      # -------------------------------------------------------------
      CONFIGURE_NEMESIS
      # -------------------------------------------------------------
   # -------------------------------------------------------------
# -------------------------------------------------------------


   
# -------------------------------------------------------------
# libHilbert -- distributed in ./contrib,
#               enabled by default
# -------------------------------------------------------------
CONFIGURE_LIBHILBERT
# -------------------------------------------------------------



# -------------------------------------------------------------
# fparser -- distributed in ./contrib,
#            enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(fparser,
              AC_HELP_STRING([--enable-fparser],
                             [build with fparser by Juha Nieminen, Joel Yliluoma]),
              enablefparser=$enableval,
              enablefparser=$enableoptional)

CONFIGURE_FPARSER
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

])