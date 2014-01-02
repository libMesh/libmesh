# -------------------------------------------------------------
# PETSc
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PETSC],
[
  AC_ARG_ENABLE(petsc,
                AC_HELP_STRING([--enable-petsc],
                               [build with PETSc iterative solver suppport]),
		[case "${enableval}" in
		  yes)  enablepetsc=yes ;;
		   no)  enablepetsc=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-petsc) ;;
		 esac],
		 [enablepetsc=$enableoptional])


  AC_ARG_VAR([PETSC_DIR],  [path to PETSc installation])
  AC_ARG_VAR([PETSC_ARCH], [PETSc build architecture])

  if (test "$enablepetsc" !=  no) ; then
    # AC_REQUIRE:
    # If the M4 macro AC_PROG_F77 has not already been called, call
    # it (without any arguments). Make sure to quote AC_PROG_F77 with
    # square brackets. AC_PROG_F77 must have been defined using
    # AC_DEFUN or else contain a call to AC_PROVIDE to indicate
    # that it has been called.
    AC_REQUIRE([AC_PROG_F77])

    # If the user doesn't have any PETSC directory specified, let's check to
    # see if it's installed via Ubuntu module
    if (test "x$PETSC_DIR" = x); then
      AC_PATH_PROG(PETSCARCH, petscarch)
      if (test "x$PETSCARCH" != x); then
        export PETSC_DIR=/usr/lib/petsc
        export PETSC_ARCH=`$PETSCARCH`
	if (test -d $PETSC_DIR); then
  	  AC_MSG_RESULT([using system-provided PETSC_DIR $PETSC_DIR])
	  AC_MSG_RESULT([using system-provided PETSC_ARCH $PETSC_ARCH])
        fi
      fi
    fi

    # Let's use a C compiler for the AC_CHECK_HEADER test, although this is
    # not strictly necessary...
    AC_LANG_PUSH(C)
    AC_CHECK_HEADER($PETSC_DIR/include/petscversion.h,
                    [enablepetsc=yes],
                    [enablepetsc=no])
    AC_LANG_POP

    # Grab PETSc version and substitute into Makefile.
    # If version 2.x, also check that PETSC_ARCH is set
    # This if-test used to be: if (test -r $PETSC_DIR/include/petsc.h) ; then
    if (test "$enablepetsc" !=  no) ; then
      # Some tricks to discover the version of petsc.
      # You have to have grep and sed for this to work.
      petscmajor=`grep "define PETSC_VERSION_MAJOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
      petscminor=`grep "define PETSC_VERSION_MINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MINOR[ ]*//g"`
      petscsubminor=`grep "define PETSC_VERSION_SUBMINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_SUBMINOR[ ]*//g"`
      petscrelease=`grep "define PETSC_VERSION_RELEASE" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_RELEASE[ ]*//g"`
      petscversion=$petscmajor.$petscminor.$petscsubminor
      petscmajorminor=$petscmajor.$petscminor.x

      AC_SUBST(petscversion)
      AC_SUBST(petscmajor)
      AC_SUBST(petscmajorminor)

      AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_MAJOR, [$petscmajor],
        [PETSc's major version number, as detected by petsc.m4])

      AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_MINOR, [$petscminor],
        [PETSc's minor version number, as detected by petsc.m4])

      AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_SUBMINOR, [$petscsubminor],
        [PETSc's subminor version number, as detected by petsc.m4])

      AC_DEFINE_UNQUOTED(DETECTED_PETSC_VERSION_RELEASE, [$petscrelease],
        [PETSc release (1) or petsc-dev (0), as detected by petsc.m4])

      if test $petscmajor = 2; then
        if test "x$PETSC_ARCH" = x ; then
          enablepetsc=no
          AC_MSG_RESULT([<<< PETSc 2.x detected and "\$PETSC_ARCH" not set.  PETSc disabled. >>>])
          # PETSc config failed.  We will try MPI at the end of this function.
        fi
      fi

    else # petscversion.h was not readable
        enablepetsc=no
    fi




    # If we haven't been disabled yet, carry on!
    if (test $enablepetsc != no) ; then

        AC_SUBST(PETSC_ARCH) # Note: may be empty...
        AC_SUBST(PETSC_DIR)
        AC_DEFINE(HAVE_PETSC, 1,
    	      [Flag indicating whether or not PETSc is available])

        # Check for snoopable MPI
        if (test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf) ; then           # 2.3.x
        	 PETSC_MPI=`grep MPIEXEC $PETSC_DIR/bmake/$PETSC_ARCH/petscconf | grep -v mpiexec.uni`
        elif (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then # 3.0.x
        	 PETSC_MPI=`grep MPIEXEC $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | grep -v mpiexec.uni`
        elif (test -r $PETSC_DIR/conf/petscvariables) ; then # 3.0.x
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

        # Print informative message about the version of PETSc we detected
        AC_MSG_RESULT([<<< Configuring library with PETSc version $petscversion support >>>])


        # If we have a full petsc distro with a makefile query it for
        # what we can
	if (test -r $PETSC_DIR/makefile); then
          PETSCLINKLIBS=`make -s -C $PETSC_DIR getlinklibs`
          PETSCINCLUDEDIRS=`make -s -C $PETSC_DIR getincludedirs`

	# create a simple makefile to provide other targets we want,
	# then query it.
 	  cat <<EOF >Makefile_config_petsc
include $PETSC_DIR/conf/variables

getPETSC_CC_INCLUDES:
	echo \$(PETSC_CC_INCLUDES)

getPETSC_FC_INCLUDES:
	echo \$(PETSC_FC_INCLUDES)
EOF
          PETSC_CC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_CC_INCLUDES`
          PETSC_FC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_FC_INCLUDES`
	  rm -f Makefile_config_petsc

  	elif (test -r $PETSC_DIR/conf/variables); then
 	  cat <<EOF >Makefile_config_petsc
include $PETSC_DIR/conf/variables
getincludedirs:
	echo -I\$(PETSC_DIR)/include -I\$(PETSC_DIR)/\$(PETSC_ARCH)/include \$(BLOCKSOLVE_INCLUDE) \$(HYPRE_INCLUDE) \$(PACKAGES_INCLUDES)

getPETSC_CC_INCLUDES:
	echo \$(PETSC_CC_INCLUDES)

getPETSC_FC_INCLUDES:
	echo \$(PETSC_FC_INCLUDES)

getlinklibs:
	echo \$(PETSC_SNES_LIB)
EOF
          PETSCLINKLIBS=`make -s -f Makefile_config_petsc getlinklibs`
          PETSCINCLUDEDIRS=`make -s -f Makefile_config_petsc getincludedirs`
          PETSC_CC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_CC_INCLUDES`
          PETSC_FC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_FC_INCLUDES`
	  rm -f Makefile_config_petsc
	fi
        #echo ""
        #echo "PETSCLINKLIBS=$PETSCLINKLIBS"
        #echo "PETSCINCLUDEDIRS=$PETSCINCLUDEDIRS"
        #echo ""

        # We sometimes need the full CC_INCLUDES to access a
        # PETSc-snooped MPI
        PETSCINCLUDEDIRS="$PETSCINCLUDEDIRS $PETSC_CC_INCLUDES"

        AC_SUBST(PETSCLINKLIBS)
        AC_SUBST(PETSCINCLUDEDIRS)
        AC_SUBST(PETSC_CC_INCLUDES)
        AC_SUBST(PETSC_FC_INCLUDES)

        AC_SUBST(MPI_IMPL)

        # Check for Hypre
        if (test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf) ; then           # 2.3.x
        	 HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/bmake/$PETSC_ARCH/petscconf`
        elif (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then # 3.0.x
        	 HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables`
        elif (test -r $PETSC_DIR/conf/petscvariables) ; then # 3.0.x
           HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/conf/petscvariables`
        fi

        if test "x$HYPRE_LIB" != x ; then
          AC_DEFINE(HAVE_PETSC_HYPRE, 1, [Flag indicating whether or not PETSc was compiled with Hypre support])
  	AC_MSG_RESULT(<<< Configuring library with Hypre support >>>)
        fi

    else
        # PETSc config failed.  Try MPI, unless directed otherwise
	if (test "$enablempi" != no); then
          AC_MSG_RESULT(<<< PETSc disabled.  Will try configuring MPI now... >>>)
          ACX_MPI
	fi
    fi


  else # --disable-petsc
    if (test "$enablempi" != no) ; then
      ACX_MPI
    fi
  fi

  AC_SUBST(enablepetsc)
])

