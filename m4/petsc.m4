# -------------------------------------------------------------
# PETSc
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PETSC],
[
  AC_ARG_ENABLE(petsc,
                AS_HELP_STRING([--disable-petsc],
                               [build without PETSc iterative solver support]),
                [case "${enableval}" in
                  yes)  enablepetsc=yes ;;
                   no)  enablepetsc=no ;;
                    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-petsc) ;;
                 esac],
                [enablepetsc=$enableoptional])

  # Setting --enable-petsc-required causes an error to be emitted
  # during configure if PETSc is not detected successfully during
  # configure.  This is useful for app codes which require PETSc (like
  # MOOSE-based apps), since it prevents situations where libmesh is
  # accidentally built without PETSc support (which may take a very
  # long time), and then the app fails to compile, requiring you to
  # redo everything.
  AC_ARG_ENABLE(petsc-required,
                AC_HELP_STRING([--enable-petsc-required],
                               [Error if PETSc is not detected by configure]),
                               [case "${enableval}" in
                     yes) petscrequired=yes ;;
                     no)  petscrequired=no ;;
                     *)   AC_MSG_ERROR(bad value ${enableval} for --enable-petsc-required) ;;
                     esac],
                     [petscrequired=no])

  # Trump --enable-petsc with --disable-mpi
  if (test "x$enablempi" = xno); then
    enablepetsc=no
  fi

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

      if test $petscmajor = 2; then
        if test "x$PETSC_ARCH" = x ; then
          enablepetsc=no
          AC_MSG_RESULT([<<< PETSc 2.x detected and "\$PETSC_ARCH" not set.  PETSc disabled. >>>])
          # PETSc config failed.  We will try MPI at the end of this function.
        fi
      fi

      # We look for petscconf.h in both $PETSC_DIR/include and
      # $PETSC_DIR/$PETSC_ARCH/include, since it can appear in either.
      petsc_use_debug=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_USE_DEBUG`
      petsc_have_superlu_dist=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_SUPERLU_DIST`
      petsc_have_mumps=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_MUMPS`
      petsc_have_metis=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_METIS`

      # We have had a slightly different way of checking for Hypre for
      # quite some time, see below.
      # petsc_have_hypre=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_HYPRE`

    else # petscversion.h was not readable
        enablepetsc=no
    fi

    # If we haven't been disabled yet, carry on!
    if (test $enablepetsc != no) ; then

        # Check for snoopable MPI
        if (test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf) ; then           # 2.3.x
          PETSC_MPI=`grep MPIEXEC $PETSC_DIR/bmake/$PETSC_ARCH/petscconf | grep -v mpiexec.uni`

        elif (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then # 3.0.x
          PETSC_MPI=`grep MPIEXEC $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | grep -v mpiexec.uni`
        elif (test -r $PETSC_DIR/conf/petscvariables) ; then # 3.0.x
          PETSC_MPI=`grep MPIEXEC $PETSC_DIR/conf/petscvariables | grep -v mpiexec.uni`

        elif (test -r $PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables) ; then # 3.6.x
          PETSC_MPI=`grep MPIEXEC $PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables | grep -v mpiexec.uni`
        elif (test -r $PETSC_DIR/lib/petsc/conf/petscvariables) ; then # 3.6.x
          PETSC_MPI=`grep MPIEXEC $PETSC_DIR/lib/petsc/conf/petscvariables | grep -v mpiexec.uni`
        fi

        # If we couldn't snoop MPI from PETSc, fall back on ACX_MPI.
        if test "x$PETSC_MPI" != x ; then
          MPI_IMPL="petsc_snooped"
          AC_MSG_RESULT(<<< Attempting to configure library with MPI from PETSC config... >>>)
        else
          AC_MSG_RESULT(<<< PETSc did not define MPIEXEC.  Will try configuring MPI now... >>>)
          ACX_MPI
        fi

        # Print informative message about the version of PETSc we detected
        AC_MSG_RESULT([<<< Found PETSc $petscversion installation in $PETSC_DIR ... >>>])

        # Figure out whether this PETSC_DIR is a PETSc source tree or an installed PETSc.
        if (test -r ${PETSC_DIR}/makefile -a -r ${PETSC_DIR}/${PETSC_ARCH}/conf/variables); then # pre-3.6.0 non-installed PETSc
          PREFIX_INSTALLED_PETSC=no
          PETSC_VARS_FILE=${PETSC_DIR}/${PETSC_ARCH}/conf/variables
        elif (test -r ${PETSC_DIR}/makefile -a -r ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/variables); then # 3.6.0+ non-installed PETSc
          PREFIX_INSTALLED_PETSC=no
          PETSC_VARS_FILE=${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/variables
        elif (test -r ${PETSC_DIR}/conf/variables); then # pre 3.6.0 prefix-installed PETSc
          PREFIX_INSTALLED_PETSC=yes
          PETSC_VARS_FILE=${PETSC_DIR}/conf/variables
        elif (test -r ${PETSC_DIR}/lib/petsc/conf/variables); then # 3.6.0 prefix-installed PETSc
          PREFIX_INSTALLED_PETSC=yes
          PETSC_VARS_FILE=${PETSC_DIR}/lib/petsc/conf/variables
        # Support having a non-prefix-installed PETSc with an
        # *incorrectly* set PETSC_ARCH environment variable.  This is
        # a less desirable configuration, but we need to support it
        # for backwards compatibility.
        elif (test -r ${PETSC_DIR}/makefile -a -r ${PETSC_DIR}/conf/variables); then # pre-3.6.0 non-installed PETSc with invalid $PETSC_ARCH
          PREFIX_INSTALLED_PETSC=no
          PETSC_VARS_FILE=${PETSC_DIR}/conf/variables
        elif (test -r ${PETSC_DIR}/makefile -a -r ${PETSC_DIR}/lib/petsc/conf/variables); then # 3.6.0+ non-installed PETSc with invalid $PETSC_ARCH
          PREFIX_INSTALLED_PETSC=no
          PETSC_VARS_FILE=${PETSC_DIR}/lib/petsc/conf/variables
        else
          AC_MSG_RESULT([<<< Could not find a viable PETSc Makefile to determine PETSC_CC_INCLUDES, etc. >>>])
          enablepetsc=no
        fi

        # We can skip the rest of the tests because they aren't going to pass
        if (test $enablepetsc != no) ; then
          if (test $PREFIX_INSTALLED_PETSC = no) ; then
            PETSCLINKLIBS=`make -s -C $PETSC_DIR getlinklibs`
            PETSCINCLUDEDIRS=`make -s -C $PETSC_DIR getincludedirs`

            printf '%s\n' "include $PETSC_VARS_FILE" > Makefile_config_petsc

            printf '%s\n' "getPETSC_CC_INCLUDES:" >> Makefile_config_petsc
            printf '\t%s\n' "echo \$(PETSC_CC_INCLUDES)" >> Makefile_config_petsc

            printf '%s\n' "getPETSC_FC_INCLUDES:" >> Makefile_config_petsc
            printf '\t%s\n' "echo \$(PETSC_FC_INCLUDES)" >> Makefile_config_petsc

            PETSC_CC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_CC_INCLUDES`
            PETSC_FC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_FC_INCLUDES`
            rm -f Makefile_config_petsc
          else
            printf '%s\n' "include $PETSC_VARS_FILE" > Makefile_config_petsc

            printf '%s\n' "getincludedirs:" >> Makefile_config_petsc
            printf '\t%s\n' "echo -I\$(PETSC_DIR)/include -I\$(PETSC_DIR)/\$(PETSC_ARCH)/include \$(BLOCKSOLVE_INCLUDE) \$(HYPRE_INCLUDE) \$(PACKAGES_INCLUDES)" >> Makefile_config_petsc

            printf '%s\n'  "getPETSC_CC_INCLUDES:" >> Makefile_config_petsc
            printf '\t%s\n' "echo \$(PETSC_CC_INCLUDES)" >> Makefile_config_petsc

            printf '%s\n' "getPETSC_FC_INCLUDES:" >> Makefile_config_petsc
            printf '\t%s\n' "echo \$(PETSC_FC_INCLUDES)" >> Makefile_config_petsc

            printf '%s\n' "getlinklibs:" >> Makefile_config_petsc
            printf '\t%s\n' "echo \$(PETSC_SNES_LIB)" >> Makefile_config_petsc

            PETSCLINKLIBS=`make -s -f Makefile_config_petsc getlinklibs`
            PETSCINCLUDEDIRS=`make -s -f Makefile_config_petsc getincludedirs`
            PETSC_CC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_CC_INCLUDES`
            PETSC_FC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_FC_INCLUDES`
            rm -f Makefile_config_petsc
          fi
        fi

        # We can skip the rest of the tests because they aren't going to pass
        if (test $enablepetsc != no) ; then

          # Debugging: see what actually got set for PETSCINCLUDEDIRS
          # echo ""
          # echo "PETSCLINKLIBS=$PETSCLINKLIBS"
          # echo "PETSCINCLUDEDIRS=$PETSCINCLUDEDIRS"
          # echo ""

          # We sometimes need the full CC_INCLUDES to access a
          # PETSc-snooped MPI
          PETSCINCLUDEDIRS="$PETSCINCLUDEDIRS $PETSC_CC_INCLUDES"

          # Check for Hypre
          if (test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf) ; then           # 2.3.x
            HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/bmake/$PETSC_ARCH/petscconf`
          elif (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then # 3.0.x
            HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables`
          elif (test -r $PETSC_DIR/conf/petscvariables) ; then # 3.0.x
            HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/conf/petscvariables`
          elif (test -r $PETSC_DIR/lib/petsc/conf/petscvariables) ; then # 3.6.x
            HYPRE_LIB=`grep "HYPRE_LIB" $PETSC_DIR/lib/petsc/conf/petscvariables`
          fi

          if test "x$HYPRE_LIB" != x ; then
            AC_MSG_RESULT(<<< Configuring library with Hypre support >>>)
          fi

          # Try to compile a trivial PETSc program to check our
          # configuration... this should handle cases where we slipped
          # by the tests above with an invalid PETSCINCLUDEDIRS
          # variable, which happened when PETSc 3.6 came out.
          AC_MSG_CHECKING(whether we can compile a trivial PETSc program)
          AC_LANG_PUSH([C++])

          # Save the original CXXFLAGS contents
          saveCXXFLAGS="$CXXFLAGS"

          # Append PETSc include paths to the CXXFLAGS variables
          CXXFLAGS="$saveCXXFLAGS $PETSCINCLUDEDIRS"

          AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
          @%:@include <petsc.h>
          static char help[]="";

          int main(int argc, char **argv)
          {
            PetscInitialize(&argc, &argv, (char*)0,help);
            PetscFinalize();
            return 0;
          }
          ]])],[
            AC_MSG_RESULT(yes)
          ],[
            AC_MSG_RESULT(no)
            enablepetsc=no
          ])

          # Return CXXFLAGS to their original state.
          CXXFLAGS="$saveCXXFLAGS"

          AC_LANG_POP([C++])


          # PETSc >= 3.5.0 should have TAO built-in, we don't currently support any other type of TAO installation.
          petsc_have_tao=no
          AC_MSG_CHECKING(for TAO support via PETSc)
          AC_LANG_PUSH([C])

          # Save the original CFLAGS contents
          saveCFLAGS="$CFLAGS"

          # Append PETSc include paths to the CFLAGS variables
          CFLAGS="$saveCFLAGS $PETSCINCLUDEDIRS"

          AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
          @%:@include <petsctao.h>
          static char help[]="";

          int main(int argc, char **argv)
          {
            Tao tao;
            PetscInitialize(&argc, &argv, (char*)0,help);
            TaoCreate(PETSC_COMM_WORLD,&tao);
            TaoDestroy(&tao);
            PetscFinalize();
            return 0;
          }
          ]])],[
            AC_MSG_RESULT(yes)
            petsc_have_tao=yes
          ],[
            AC_MSG_RESULT(no)
          ])

          # Return C flags to their original state.
          CFLAGS="$saveCFLAGS"

          AC_LANG_POP([C])
        fi # if (test $enablepetsc != no)

        if (test $enablepetsc = no); then
          # PETSc config failed.  Try MPI, unless directed otherwise
          if (test "$enablempi" != no); then
            AC_MSG_RESULT(<<< PETSc disabled.  Will try configuring MPI now... >>>)
            ACX_MPI
          fi
        fi # if (test $enablepetsc = no)
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

  # If PETSc is still enabled at this point, do the required AC_SUBST
  # and AC_DEFINE commands.  This prevents libmesh_config.h from having
  # confusing information in it if the test compilation steps fail.
  if (test $enablepetsc != no) ; then
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

    # Set a #define if PETSc was built with debugging enabled.  Note
    # that this token will appear as LIBMESH_PETSC_USE_DEBUG in our
    # header, so it won't collide with PETSc's.  It is safe to test
    # the value of $petsc_use_debug, it is guaranteed to be either 0
    # or nonzero.
    if (test $petsc_use_debug -gt 0) ; then
      AC_DEFINE(PETSC_USE_DEBUG, 1, [Flag indicating whether or not PETSc was configured with debugging enabled])
    fi

    # Set a #define if PETSc was built with SuperLU_dist support
    if (test $petsc_have_superlu_dist -gt 0) ; then
      AC_DEFINE(PETSC_HAVE_SUPERLU_DIST, 1, [Flag indicating whether or not PETSc was configured with SuperLU_dist support])
    fi

    # Set a #define if PETSc was built with MUMPS support
    if (test $petsc_have_mumps -gt 0) ; then
      AC_DEFINE(PETSC_HAVE_MUMPS, 1, [Flag indicating whether or not PETSc was configured with MUMPS support])
    fi

    # Set a #define if PETSc was built with Hypre support.  We have
    # been using the slightly different "LIBMESH_HAVE_PETSC_HYPRE" for
    # quite a while, I may update this in a future branch.
    # if (test $petsc_have_hypre -gt 0) ; then
    #   AC_DEFINE(PETSC_HAVE_HYPRE, 1, [Flag indicating whether or not PETSc was configured with Hypre support])
    # fi

    AC_SUBST(PETSC_ARCH) # Note: may be empty...
    AC_SUBST(PETSC_DIR)

    AC_SUBST(PETSCLINKLIBS)
    AC_SUBST(PETSCINCLUDEDIRS)
    AC_SUBST(PETSC_CC_INCLUDES)
    AC_SUBST(PETSC_FC_INCLUDES)
    AC_SUBST(MPI_IMPL)

    if test "x$PETSC_MPI" != x ; then
      AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
    fi

    if test "x$HYPRE_LIB" != x ; then
      AC_DEFINE(HAVE_PETSC_HYPRE, 1, [Flag indicating whether or not PETSc was compiled with Hypre support])
    fi

    AC_DEFINE(HAVE_PETSC, 1, [Flag indicating whether or not PETSc is available])

    if (test $petsc_have_tao != no) ; then
      AC_DEFINE(HAVE_PETSC_TAO, 1, [Flag indicating whether or not the Toolkit for Advanced Optimization (TAO) is available via PETSc])
    fi
  fi

  # If PETSc is not enabled, but it *was* required, error out now
  # instead of compiling libmesh in an invalid configuration.
  if (test $enablepetsc = no -a $petscrequired = yes) ; then
    dnl We return error code 3 here, since 0 means success and 1 is
    dnl indistinguishable from other errors.  Ideally, all of the
    dnl AC_MSG_ERROR calls in our m4 files would return a different
    dnl error code, but currently this is not implemented.
    AC_MSG_ERROR([*** PETSc was not found, but --enable-petsc-required was specified.], 3)
  fi
])
