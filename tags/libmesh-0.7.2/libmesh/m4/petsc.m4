dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl PETSc
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PETSC], 
[
  dnl AC_REQUIRE:
  dnl If the M4 macro AC_PROG_F77 has not already been called, call
  dnl it (without any arguments). Make sure to quote AC_PROG_F77 with
  dnl square brackets. AC_PROG_F77 must have been defined using
  dnl AC_DEFUN or else contain a call to AC_PROVIDE to indicate
  dnl that it has been called.
  AC_REQUIRE([AC_PROG_F77])

  dnl If the user doesn't have any PETSC directory specified, let's check to
  dnl see if it's installed via Ubuntu module
  if test "x$PETSC_DIR" = x ; then
    PETSC_DIR=/usr/lib/petsc
    PETSC_ARCH=linux-gnu-c-opt
  fi

  AC_CHECK_FILE($PETSC_DIR/include/petsc.h,
                PETSC_H_PATH=$PETSC_DIR/include/petsc.h)

  dnl Grab PETSc version and substitute into Makefile.
  dnl If version 2.x, also check that PETSC_ARCH is set
  if (test -r $PETSC_DIR/include/petsc.h) ; then
    dnl Some tricks to discover the version of petsc.
    dnl You have to have grep and sed for this to work.
    petscmajor=`grep "define PETSC_VERSION_MAJOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
    petscminor=`grep "define PETSC_VERSION_MINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MINOR[ ]*//g"`
    petscsubminor=`grep "define PETSC_VERSION_SUBMINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_SUBMINOR[ ]*//g"`
    petscversion=$petscmajor.$petscminor.$petscsubminor
    petscmajorminor=$petscmajor.$petscminor.x

    AC_SUBST(petscversion)
    AC_SUBST(petscmajor)
    AC_SUBST(petscmajorminor)

    AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_MAJOR, [$petscmajor],
      [PETSc's major version number, as detected by LibMesh])
	      
    AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_MINOR, [$petscminor],
      [PETSc's minor version number, as detected by LibMesh])
	      
    AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_SUBMINOR, [$petscsubminor],
      [PETSc's subminor version number, as detected by LibMesh])

    if test $petscmajor = 2; then
      if test "x$PETSC_ARCH" = x ; then
        enablepetsc=no
        AC_MSG_RESULT([<<< PETSc 2.x detected and "\$PETSC_ARCH" not set.  PETSc disabled. >>>])
        dnl PETSc config failed.  We will try MPI at the end of this function.
        dnl ACX_MPI
      fi
    fi

  else dnl petsc.h was not readable
      enablepetsc=no
  fi




  dnl If we haven't been disabled yet, carry on!
  if (test $enablepetsc != no) ; then
dnl      dnl look for a decent F77 compiler or honor --with-77=...
dnl      FC_TRY_LIST="g77 gfortran ifort f77 xlf frt pgf77 fort77 fl32 af77 f90 xlf90 pgf90 epcf90 f95 fort xlf95 ifc efc pgf95 lf95"
dnl      AC_ARG_WITH([f77],
dnl      	    AC_HELP_STRING([--with-f77=F77],
dnl                                 [Fortran compiler to use]),
dnl      	    [F77="$withval"],
dnl      	    [])
dnl      
dnl      dnl --------------------------------------------------------------
dnl      dnl Determine a F77 compiler to use.
dnl      dnl --------------------------------------------------------------
dnl      AC_PROG_F77([$F77_TRY_LIST])
dnl      AC_F77_LIBRARY_LDFLAGS

      AC_SUBST(PETSC_ARCH) dnl Note: may be empty...
      AC_SUBST(PETSC_DIR)
      AC_DEFINE(HAVE_PETSC, 1,
  	      [Flag indicating whether or not PETSc is available])

      dnl Check for snoopable MPI
      if (test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf) ; then           dnl 2.3.x	
      	 PETSC_MPI=`grep MPIEXEC $PETSC_DIR/bmake/$PETSC_ARCH/petscconf | grep -v mpiexec.uni` 
      elif (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then dnl 3.0.x
      	 PETSC_MPI=`grep MPIEXEC $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | grep -v mpiexec.uni`
      elif (test -r $PETSC_DIR/conf/petscvariables) ; then dnl 3.0.x
      	 PETSC_MPI=`grep MPIEXEC $PETSC_DIR/conf/petscvariables | grep -v mpiexec.uni`
      fi		 
      if test "x$PETSC_MPI" != x ; then
        AC_DEFINE(HAVE_MPI, 1,
  	        [Flag indicating whether or not MPI is available])
        MPI_IMPL="petsc_snooped"      
	AC_MSG_RESULT(<<< Configuring library with MPI from PETSC config >>>)
      else
	AC_MSG_RESULT(<<< Warning: configuring in serial - no MPI in PETSC config >>>)
      fi

      dnl Print informative message about the version of PETSc we detected
      AC_MSG_RESULT([<<< Configuring library with PETSc version $petscversion support >>>])

dnl      PETSCLINKLIBS=`cd $PETSC_DIR ; make getlinklibs`
dnl      PETSCINCLUDEDIRS=`cd $PETSC_DIR ; make getincludedirs`
dnl
dnl      AC_SUBST(PETSCLINKLIBS)
dnl      AC_SUBST(PETSCINCLUDEDIRS)

      AC_SUBST(MPI_IMPL)

      dnl Check for Hypre
      if (test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf) ; then           dnl 2.3.x	
      	 HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/bmake/$PETSC_ARCH/petscconf` 
      elif (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then dnl 3.0.x
      	 HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables`
      elif (test -r $PETSC_DIR/conf/petscvariables) ; then dnl 3.0.x 
         HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/conf/petscvariables`
      fi
		 
      if test "x$HYPRE_LIB" != x ; then
        AC_DEFINE(HAVE_PETSC_HYPRE, 1, [Flag indicating whether or not PETSc was compiled with Hypre support])
	AC_MSG_RESULT(<<< Configuring library with Hypre support >>>)
      fi
  
  else 
      dnl PETSc config failed.  Try MPI.
      AC_MSG_RESULT(<<< PETSc disabled.  Will try configuring MPI now... >>>)
      ACX_MPI
  
  fi 

  
  AC_SUBST(enablepetsc)
])



dnl ----------------------------------------------------------------------------
dnl check for the required PETSc library
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_PETSc], [
AC_REQUIRE([ACX_LAPACK])
BLAS_LIBS="$BLAS_LIBS $FLIBS"
LAPACK_LIBS="$LAPACK_LIBS $BLAS_LIBS"
AC_PATH_XTRA
X_LIBS="$X_PRE_LIBS $X_LIBS -lX11 $X_EXTRA_LIBS"

# Set variables...
AC_ARG_WITH([PETSc],
	    AC_ARG_HELP([--with-PETSc=PATH],
                        [Prefix where PETSc is installed (PETSC_DIR)]),
	    [PETSc="$withval"],
	    [
              if test $PETSC_DIR; then
		PETSc="$PETSC_DIR"
		echo "note: assuming PETSc library is in $PETSc (/lib,/include) as specified by environment variable PETSC_DIR"
	      else
		PETSc="/usr/local"
		echo "note: assuming PETSc library is in /usr/local (/lib,/include)"
	      fi
            ])

AC_ARG_WITH([BOPT],
	    AC_ARG_HELP([--with-BOPT=VAL],[BOPT setting for PETSc (BOPT)]),
 	    [BOPT="$withval"],
	    [
              echo "note: assuming BOPT to O"
	      BOPT="O"
            ])

AC_ARG_WITH([PETSc_ARCH],
	    AC_ARG_HELP([--with-PETSc_ARCH=VAL],[PETSc hardware architecture (PETSC_ARCH)]),
	    [PETSc_ARCH="$withval"],
	    [
              if test $PETSC_ARCH; then
		PETSc_ARCH="$PETSC_ARCH"
		echo "note: assuming PETSc hardware architecture to be $PETSc_ARCH as specified by environment variable PETSC_ARCH"
	      else
		PETSc_ARCH=`uname -p`
		echo "note: assuming PETSc hardware architecture to be $PETSc_ARCH"
	      fi
            ])

PETSc_LIBS_PATH="$PETSc/lib/lib$BOPT/$PETSc_ARCH"
PETSc_INCLUDES_PATH="$PETSc/include"

# Check that the compiler uses the library we specified...
if test -e $PETSc_LIBS_PATH/libpetsc.a || test -e $PETSc_LIBS_PATH/libpetsc.so; then
	echo "note: using $PETSc_LIBS_PATH/libpetsc (.a/.so)"
else
	AC_MSG_ERROR( [Could not physically find PETSc library... exiting] )
fi 
if test -e $PETSc_INCLUDES_PATH/petsc.h; then
	echo "note: using $PETSc_INCLUDES_PATH/petsc.h"
else
	AC_MSG_ERROR( [Could not physically find PETSc header file... exiting] )
fi 

# Ensure the comiler finds the library...
tmpLIBS=$LIBS
tmpCPPFLAGS=$CPPFLAGS
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_CHECK_LIB(
	[dl],
	[dlopen],
	[DL_LIBS="-ldl"],
	[DL_LIBS=""; echo "libdl not found, assuming not needed for this architecture"] )
LIBS="-L$PETSc_LIBS_PATH $MPI_LIBS_PATHS $MPI_LIBS $LAPACK_LIBS $X_LIBS $LIBS -lm $DL_LIBS"
CPPFLAGS="$MPI_INCLUDES_PATHS -I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH $CPPFLAGS"
echo "cppflags=$CPPFLAGS"

AC_CHECK_LIB(
	[petsc],
	[PetscError],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc library... exiting] )] )
AC_CHECK_LIB(
	[petscvec],
	[ISCreateGeneral],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscvec library... exiting] )] )
AC_CHECK_LIB(
	[petscmat],
	[MAT_Copy],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscmat library... exiting] )] )
AC_CHECK_LIB(
	[petscdm],
	[DMInitializePackage],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscdm library... exiting] )] )
AC_CHECK_LIB(
	[petscsles],
	[SLESCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscsles library... exiting] )] )
AC_CHECK_LIB(
	[petscsnes],
	[SNESCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscsnes library... exiting] )] )
AC_CHECK_LIB(
	[petscts],
	[TSCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscts library... exiting] )] )
AC_CHECK_LIB(
	[petscmesh],
	[MESH_CreateFullCSR],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscmesh library... exiting] )] )
AC_CHECK_LIB(
	[petscgrid],
	[GridCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscgrid library... exiting] )] )
AC_CHECK_LIB(
	[petscgsolver],
	[GSolverInitializePackage],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscgsolver library... exiting] )] )
AC_CHECK_LIB(
	[petscfortran],
	[meshcreate_],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc library... exiting] )] )
	AC_CHECK_LIB(
	[petsccontrib],
	[SDACreate1d],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petsccontrib library... exiting] )] )
AC_CHECK_HEADER(
	[petsc.h],
	[AC_DEFINE( 
		[HAVE_PETSC],,
		[Define to 1 if you have the <petsc.h> header file.])],
	[AC_MSG_ERROR( [Could not compile in the PETSc headers... exiting] )] )
PETSc_LIBS="-lpetsc -lpetscvec -lpetscmat -lpetscdm -lpetscsles -lpetscsnes \
	-lpetscts -lpetscmesh -lpetscgrid -lpetscgsolver -lpetscfortran -lpetsccontrib \
	$PETSc_ARCH_LIBS"
PETSc_LIBS_PATHS="-L$PETSc_LIBS_PATH"
PETSc_INCLUDES_PATHS="-I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH"

# Save variables...
AC_LANG_RESTORE
LIBS=$tmpLIBS
CPPFLAGS=$tmpCPPFLAGS
AC_SUBST( PETSc )
AC_SUBST( PETSc_IMPL )
AC_SUBST( PETSc_LIBS )
AC_SUBST( PETSc_LIBS_PATH )
AC_SUBST( PETSc_LIBS_PATHS )
AC_SUBST( PETSc_INCLUDES_PATH )
AC_SUBST( PETSc_INCLUDES_PATHS )
])dnl ACX_PETSc ------------------------------------------------------------
