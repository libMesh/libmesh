# -------------------------------------------------------------
# PETSc
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PETSC],
[
  dnl Require having scraped PETSc CXX, includes, and libs
  AC_REQUIRE([ACSM_SCRAPE_PETSC_CONFIGURE])
  dnl Require having set a compiler so we can check headers,
  dnl try compiles, try links, etc.
  AC_REQUIRE([LIBMESH_SET_COMPILERS])

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
                [AS_CASE("${enableval}",
                         [yes], [petscrequired=yes],
                         [no],  [petscrequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-petsc-required)])],
                     [petscrequired=no])

  # Setting --enable-petsc-hypre-required causes an error to be
  # emitted during configure if PETSc with builtin Hypre is not
  # detected successfully.  This is useful for app codes which require
  # both PETSc and Hypre (like MOOSE-based apps), since it prevents
  # libmesh from being accidentally built without PETSc and Hypre
  # support.
  AC_ARG_ENABLE(petsc-hypre-required,
                AC_HELP_STRING([--enable-petsc-hypre-required],
                               [Error if a PETSc with Hypre is not detected by configure]),
                [AS_CASE("${enableval}",
                         [yes], [petschyprerequired=yes
                                 petscrequired=yes],
                         [no],  [petschyprerequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-petsc-hypre-required)])],
                     [petschyprerequired=no])

  AS_IF([test "$enablepetsc" !=  no],
        [
    # AC_REQUIRE:
    # If the M4 macro AC_PROG_F77 has not already been called, call
    # it (without any arguments). Make sure to quote AC_PROG_F77 with
    # square brackets. AC_PROG_F77 must have been defined using
    # AC_DEFUN or else contain a call to AC_PROVIDE to indicate
    # that it has been called.
    AC_REQUIRE([AC_PROG_F77])


    # Let's use a C compiler for the AC_CHECK_HEADER test, although this is
    # not strictly necessary...
    AC_LANG_PUSH(C)
    AC_CHECK_HEADER(${PETSC_DIR}/${PETSC_ARCH}/include/petscversion.h,
                    [enablepetsc=yes; enablepetsc_mpi=yes
                     petsc_version_h_file=${PETSC_DIR}/${PETSC_ARCH}/include/petscversion.h],
                    [
      AC_CHECK_HEADER(${PETSC_DIR}/include/petscversion.h,
                      [enablepetsc=yes; enablepetsc_mpi=yes
                       petsc_version_h_file=${PETSC_DIR}/include/petscversion.h],
                      [enablepetsc=no; enablepetsc_mpi=no])
                    ])
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
            petscmajor=`grep "define PETSC_VERSION_MAJOR" $petsc_version_h_file | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
            petscminor=`grep "define PETSC_VERSION_MINOR" $petsc_version_h_file | sed -e "s/#define PETSC_VERSION_MINOR[ ]*//g"`
            petscsubminor=`grep "define PETSC_VERSION_SUBMINOR" $petsc_version_h_file | sed -e "s/#define PETSC_VERSION_SUBMINOR[ ]*//g"`
            petscrelease=`grep "define PETSC_VERSION_RELEASE" $petsc_version_h_file | sed -e "s/#define PETSC_VERSION_RELEASE[ ]*//g"`
            petscversion=$petscmajor.$petscminor.$petscsubminor
            petscmajorminor=$petscmajor.$petscminor.x

            AS_IF([test "$petscmajor" = "2" && test "x${PETSC_ARCH}" = "x"],
                  [
                    dnl PETSc config failed.  We will try MPI at the end of this function.
                    enablepetsc=no; enablepetsc_mpi=no
                    AC_MSG_RESULT([<<< PETSc 2.x detected and "\$PETSC_ARCH" not set.  PETSc disabled. >>>])
                  ])

            dnl We look for petscconf.h in both $PETSC_DIR/include and
            dnl $PETSC_DIR/$PETSC_ARCH/include, since it can appear in either.
            petsc_use_debug=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_USE_DEBUG`
            petsc_have_superlu_dist=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_SUPERLU_DIST`
            petsc_have_mumps=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_MUMPS`
            petsc_have_metis=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_METIS`
            petsc_have_chaco=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_CHACO`
            petsc_have_party=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_PARTY`
            petsc_have_ptscotch=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_PTSCOTCH`
            petsc_have_parmetis=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_PARMETIS`
            petsc_have_hypre=`cat ${PETSC_DIR}/include/petscconf.h ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null | grep -c PETSC_HAVE_HYPRE`
          ],
          [enablepetsc=no; enablepetsc_mpi=no])

    # If we haven't been disabled yet, carry on!
    AS_IF([test $enablepetsc != no],
          [

        # Print informative message about the version of PETSc we detected
        AC_MSG_RESULT([<<< Found PETSc $petscversion installation in ${PETSC_DIR} ... >>>])

        # We can skip the rest of the tests because they aren't going to pass
        AS_IF([test $enablepetsc != no],
              [
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
            enablepetsc=no; enablepetsc_mpi=no
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
        ])
    ])
  ])

  AC_SUBST(enablepetsc)

  # If PETSc is still enabled at this point, do the required AC_SUBST
  # and AC_DEFINE commands.  This prevents libmesh_config.h from having
  # confusing information in it if the test compilation steps fail.
  AS_IF([test $enablepetsc != no],
        [
          AC_SUBST(petscversion)
          AC_SUBST(petscmajor)
          AC_SUBST(petscminor)
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
          AS_IF([test $petsc_use_debug -gt 0],
                [AC_DEFINE(PETSC_USE_DEBUG, 1, [Flag indicating whether or not PETSc was configured with debugging enabled])])

          # Set a #define if PETSc was built with SuperLU_dist support
          AS_IF([test $petsc_have_superlu_dist -gt 0],
                [AC_DEFINE(PETSC_HAVE_SUPERLU_DIST, 1, [Flag indicating whether or not PETSc was configured with SuperLU_dist support])])

          # Set a #define if PETSc was built with MUMPS support
          AS_IF([test $petsc_have_mumps -gt 0],
                [AC_DEFINE(PETSC_HAVE_MUMPS, 1, [Flag indicating whether or not PETSc was configured with MUMPS support])])

          # Set a #define if PETSc was built with Chaco support
          AS_IF([test $petsc_have_chaco -gt 0],
                [AC_DEFINE(PETSC_HAVE_CHACO, 1, [Flag indicating whether or not PETSc was configured with CHACO support])])

          # Set a #define if PETSc was built with Party support
          AS_IF([test $petsc_have_party -gt 0],
                [AC_DEFINE(PETSC_HAVE_PARTY, 1, [Flag indicating whether or not PETSc was configured with PARTY support])])

          # Set a #define if PETSc was built with PTScotch support
          AS_IF([test $petsc_have_ptscotch -gt 0],
                [AC_DEFINE(PETSC_HAVE_PTSCOTCH, 1, [Flag indicating whether or not PETSc was configured with PTSCOTCH support])])

          # Set a #define if PETSc was built with ParMETIS support
          AS_IF([test $petsc_have_parmetis -gt 0],
                [AC_DEFINE(PETSC_HAVE_PARMETIS, 1, [Flag indicating whether or not PETSc was configured with ParMETIS support])])

          AC_SUBST(PETSC_ARCH) # Note: may be empty...
          AC_SUBST(PETSC_DIR)

          AC_SUBST(PETSCLINKLIBS)
          AC_SUBST(PETSCINCLUDEDIRS)
          AC_SUBST(PETSC_CC_INCLUDES)
          AC_SUBST(PETSC_FC_INCLUDES)

          AS_IF([test $petsc_have_hypre -gt 0],
                [AC_DEFINE(HAVE_PETSC_HYPRE, 1, [Flag indicating whether or not PETSc was compiled with Hypre support])])

          AC_DEFINE(HAVE_PETSC, 1, [Flag indicating whether or not PETSc is available])

          AS_IF([test $petsc_have_tao != no],
                [AC_DEFINE(HAVE_PETSC_TAO, 1, [Flag indicating whether or not the Toolkit for Advanced Optimization (TAO) is available via PETSc])])
        ])

  # If PETSc is not enabled, but it *was* required, error out now
  # instead of compiling libmesh in an invalid configuration.
  AS_IF([test "$enablepetsc" = "no" && test "$petscrequired" = "yes"],
        dnl We return error code 3 here, since 0 means success and 1 is
        dnl indistinguishable from other errors.  Ideally, all of the
        dnl AC_MSG_ERROR calls in our m4 files would return a different
        dnl error code, but currently this is not implemented.
        [AC_MSG_ERROR([*** PETSc was not found, but --enable-petsc-required was specified.], 3)])

  # If PETSc + Hypre is required, throw an error if we don't have it.
  AS_IF([test "x$petschyprerequired" = "xyes" && test $petsc_have_hypre -eq 0],
        [AC_MSG_ERROR([*** PETSc with Hypre was not found, but --enable-petsc-hypre-required was specified.], 4)])
])
