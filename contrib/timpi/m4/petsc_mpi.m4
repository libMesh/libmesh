# -------------------------------------------------------------
# PETSc MPI detection
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PETSC_MPI],
[
  AC_ARG_ENABLE(petsc,
                AS_HELP_STRING([--disable-petsc],
                               [build without trying to use PETSc MPI]),
                [AS_CASE("${enableval}",
                         [yes], [enablepetsc=yes],
                         [no],  [enablepetsc=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-petsc)])],
                [enablepetsc=$enableoptional])

  AC_ARG_VAR([PETSC_DIR],  [path to PETSc installation])
  AC_ARG_VAR([PETSC_ARCH], [PETSc build architecture])

  AS_IF([test "$enablepetsc" !=  no],
        [
    AS_ECHO(["Trying to find MPI via PETSc configuration"])

    # AC_REQUIRE:
    # If the M4 macro AC_PROG_F77 has not already been called, call
    # it (without any arguments). Make sure to quote AC_PROG_F77 with
    # square brackets. AC_PROG_F77 must have been defined using
    # AC_DEFUN or else contain a call to AC_PROVIDE to indicate
    # that it has been called.
    AC_REQUIRE([AC_PROG_F77])

    # If the user doesn't have any PETSC directory specified, let's check to
    # see if it's installed via Ubuntu module
    AS_IF([test "x$PETSC_DIR" = x],
          [
            AC_PATH_PROG(PETSCARCH, petscarch)
            AS_IF([test "x$PETSCARCH" != x],
                  [
                    export PETSC_DIR=/usr/lib/petsc
                    export PETSC_ARCH=`$PETSCARCH`
                    AS_IF([test -d $PETSC_DIR],
                          [
                            AC_MSG_RESULT([using system-provided PETSC_DIR $PETSC_DIR])
                            AC_MSG_RESULT([using system-provided PETSC_ARCH $PETSC_ARCH])
                          ])
                  ])
          ])

    # Let's use a C compiler for the AC_CHECK_HEADER test, although this is
    # not strictly necessary...
    AC_LANG_PUSH(C)
    AC_CHECK_HEADER($PETSC_DIR/include/petscversion.h,
                    [enablepetsc=yes],
                    [enablepetsc=no])
    AC_LANG_POP

    dnl We now have a -gt check for this that occurs at the end of the file, so make sure
    dnl it is initialized to some sensible value to avoid syntax errors.
    petsc_have_hypre=0

    # Grab PETSc version and substitute into Makefile.
    # If version 2.x, also check that PETSC_ARCH is set
    AS_IF([test "$enablepetsc" !=  no],
          [
            dnl Some tricks to discover the version of petsc.
            dnl You have to have grep and sed for this to work.
            petscmajor=`grep "define PETSC_VERSION_MAJOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
            petscminor=`grep "define PETSC_VERSION_MINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MINOR[ ]*//g"`
            petscsubminor=`grep "define PETSC_VERSION_SUBMINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_SUBMINOR[ ]*//g"`
            petscrelease=`grep "define PETSC_VERSION_RELEASE" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_RELEASE[ ]*//g"`
            petscversion=$petscmajor.$petscminor.$petscsubminor
            petscmajorminor=$petscmajor.$petscminor.x

            AS_IF([test "$petscmajor" = "2" && test "x$PETSC_ARCH" = "x"],
                  [
                    dnl PETSc config failed.
                    enablepetsc=no
                    AC_MSG_RESULT([<<< PETSc 2.x detected and "\$PETSC_ARCH" not set.  PETSc disabled. >>>])
                  ])
          ],
          [enablepetsc=no])

    # If we haven't been disabled yet, carry on!
    AS_IF([test $enablepetsc != no],
          [
        dnl Check for snoopable MPI
        AS_IF([test -r $PETSC_DIR/bmake/$PETSC_ARCH/petscconf], dnl 2.3.x
              [PETSC_MPI=`grep MPIEXEC $PETSC_DIR/bmake/$PETSC_ARCH/petscconf | grep -v mpiexec.uni`],
              [test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables], dnl 3.0.x
              [PETSC_MPI=`grep MPIEXEC $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | grep -v mpiexec.uni`],
              [test -r $PETSC_DIR/conf/petscvariables], dnl 3.0.x
              [PETSC_MPI=`grep MPIEXEC $PETSC_DIR/conf/petscvariables | grep -v mpiexec.uni`],
              [test -r $PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables], dnl 3.6.x
              [PETSC_MPI=`grep MPIEXEC $PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables | grep -v mpiexec.uni`],
              [test -r $PETSC_DIR/lib/petsc/conf/petscvariables], dnl 3.6.x
              [PETSC_MPI=`grep MPIEXEC $PETSC_DIR/lib/petsc/conf/petscvariables | grep -v mpiexec.uni`])

        AS_IF([test "x$PETSC_MPI" != x],
              [
                MPI_IMPL="petsc_snooped"
                AC_MSG_RESULT(<<< Attempting to configure library with MPI from PETSC config... >>>)
              ])

        # Print informative message about the version of PETSc we detected
        AC_MSG_RESULT([<<< Found PETSc $petscversion installation in $PETSC_DIR ... >>>])

        # Figure out whether this PETSC_DIR is a PETSc source tree or an installed PETSc.
        AS_IF(dnl pre-3.6.0 non-installed PETSc
              [test -r ${PETSC_DIR}/makefile && test -r ${PETSC_DIR}/${PETSC_ARCH}/conf/variables],
              [PREFIX_INSTALLED_PETSC=no
               PETSC_VARS_FILE=${PETSC_DIR}/${PETSC_ARCH}/conf/variables],

              dnl 3.6.0+ non-installed PETSc
              [test -r ${PETSC_DIR}/makefile && test -r ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/variables],
              [PREFIX_INSTALLED_PETSC=no
               PETSC_VARS_FILE=${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/variables],

              dnl pre 3.6.0 prefix-installed PETSc
              [test -r ${PETSC_DIR}/conf/variables],
              [PREFIX_INSTALLED_PETSC=yes
               PETSC_VARS_FILE=${PETSC_DIR}/conf/variables],

              dnl 3.6.0 prefix-installed PETSc
              [test -r ${PETSC_DIR}/lib/petsc/conf/variables],
              [PREFIX_INSTALLED_PETSC=yes
               PETSC_VARS_FILE=${PETSC_DIR}/lib/petsc/conf/variables],

              dnl Support having a non-prefix-installed PETSc with an
              dnl *incorrectly* set PETSC_ARCH environment variable.  This is
              dnl a less desirable configuration, but we need to support it
              dnl for backwards compatibility.

              dnl pre-3.6.0 non-installed PETSc with invalid $PETSC_ARCH
              [test -r ${PETSC_DIR}/makefile && test -r ${PETSC_DIR}/conf/variables],
              [PREFIX_INSTALLED_PETSC=no
               PETSC_VARS_FILE=${PETSC_DIR}/conf/variables],

              dnl 3.6.0+ non-installed PETSc with invalid $PETSC_ARCH
              [test -r ${PETSC_DIR}/makefile && test -r ${PETSC_DIR}/lib/petsc/conf/variables],
              [PREFIX_INSTALLED_PETSC=no
               PETSC_VARS_FILE=${PETSC_DIR}/lib/petsc/conf/variables],

              dnl If nothing else matched
              [AC_MSG_RESULT([<<< Could not find a viable PETSc Makefile to determine PETSC_CC_INCLUDES, etc. >>>])
               enablepetsc=no])

        dnl Set some include and link variables by building and running temporary Makefiles.
        AS_IF([test "$enablepetsc" != "no" && test "$PREFIX_INSTALLED_PETSC" = "no"],
              [
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
              ],
              [
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
              ])

        # We can skip the rest of the tests because they aren't going to pass
        AS_IF([test $enablepetsc != no],
              [
          dnl Debugging: see what actually got set for PETSCINCLUDEDIRS
          dnl echo ""
          dnl echo "PETSCLINKLIBS=$PETSCLINKLIBS"
          dnl echo "PETSCINCLUDEDIRS=$PETSCINCLUDEDIRS"
          dnl echo ""

          # We sometimes need the full CC_INCLUDES to access a
          # PETSc-snooped MPI
          PETSCINCLUDEDIRS="$PETSCINCLUDEDIRS $PETSC_CC_INCLUDES"

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
        ])
    ])
  ])

  AC_SUBST(enablepetsc)

  # If PETSc is still enabled at this point, do the required AC_SUBST
  # and AC_DEFINE commands.  This prevents libmesh_config.h from having
  # confusing information in it if the test compilation steps fail.
  AS_IF([test $enablepetsc != no],
        [
          AC_SUBST(PETSCLINKLIBS)
          AC_SUBST(PETSCINCLUDEDIRS)
          AC_SUBST(PETSC_CC_INCLUDES)
          AC_SUBST(MPI_IMPL)

          AS_IF([test "x$PETSC_MPI" != x],
                [AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                 enablempi=yes
                ])
        ])
])
