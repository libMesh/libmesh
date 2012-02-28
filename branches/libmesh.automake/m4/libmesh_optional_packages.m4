# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_CONFIGURE_OPTIONAL_PACKAGES],
[

libmesh_optional_INCLUDES=""

# --------------------------------------------------------------
# Allow for disable-optional
# --------------------------------------------------------------
AC_ARG_ENABLE(optional,
              AC_HELP_STRING([--enable-optional],
                             [en/disable optional external libraries]),
              enableoptional=$enableval,
              enableoptional=yes)
AC_SUBST(enableoptional)	

if test "$enableoptional" != no ; then
   AC_MSG_RESULT(---------------------------------------------)
   AC_MSG_RESULT(----- Configuring for optional packages -----)
   AC_MSG_RESULT(---------------------------------------------)
fi


# -------------------------------------------------------------
# Petsc -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(petsc,
              AC_HELP_STRING([--enable-petsc],
                             [build with PETSc iterative solver suppport]),
              enablepetsc=$enableval,
              enablepetsc=$enableoptional)


# -------------------------------------------------------------
# configure MPI separately only if PETSc fails
if (test "$enablepetsc" !=  no) ; then
  CONFIGURE_PETSC
else
  if (test "$enableoptional" != no) ; then
    ACX_MPI
  fi	
fi
# -------------------------------------------------------------




# -------------------------------------------------------------
# SLEPc  -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(slepc,
              AC_HELP_STRING([--enable-slepc],
                             [build with SLEPc eigen solver support]),
     enableslepc=$enableval,
     enableslepc=$enableoptional)

if test "$enableslepc" != no ; then
   CONFIGURE_SLEPC
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# "Trilinos" -- enabled by default unless we're building with
#               complex numbers. This looks for the AztecOO
#               linear solvers and requisite Epetra parallel 
#               data structures
# -------------------------------------------------------------
AC_ARG_ENABLE(trilinos,
              AC_HELP_STRING([--enable-trilinos],
                             [build with Trilinos support]),
              enabletrilinos=$enableval,
              enabletrilinos=$enableoptional)
if test "$enablecomplex" = no ; then
   if test "$enabletrilinos" != no ; then
      # -- try Trilinos 10 first
      CONFIGURE_TRILINOS_10
      # -- then Trlinos 9
      if test "$enabletrilinos10" = no ; then
        CONFIGURE_TRILINOS_9
      fi
   fi
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Intel's Threading Building Blocks  -- enabled by default
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

# substitute values
AC_SUBST(libmesh_optional_INCLUDES)

])