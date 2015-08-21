# -------------------------------------------------------------
# SLEPc
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SLEPC],
[
  AC_ARG_ENABLE(slepc,
                AS_HELP_STRING([--disable-slepc],
                               [build without SLEPc eigen solver support]),
                [case "${enableval}" in
                  yes)  enableslepc=yes ;;
                  no)  enableslepc=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-slepc) ;;
                esac],
                [enableslepc=$enablepetsc]) # if unspecified, infer from what's up with PETSc

  AC_ARG_VAR([SLEPC_DIR], [path to SLEPc installation])

  if (test "$enableslepc" !=  no) ; then

    # Although PETSc 3.x does not require PETSC_ARCH, it appears that
    # SLEPc 3.x may (based on patches sent in by Jed Brown)?
    # If that's the case, we should probably uncomment
    # this check and disable SLEPc when PETSC_ARCH is not set...
    # if test "x$PETSC_ARCH" = x ; then
    #   enableslepc=no
    #   AC_MSG_RESULT([<<< SLEPc disabled.  Please set your "\$PETSC_ARCH" environment variable correctly. >>>])
    # fi

    # Test to see if SLEPC_DIR set by user.  If not set, then
    # try to autodetect in a default directory
    if test "x$SLEPC_DIR" = x ; then
      export SLEPC_DIR=/usr/lib/slepc
    fi

    # Test to see if SLEPC_DIR and petscversion were set by user or
    # autodetection.  If not set, then disable slepc, print a message.
    if (test "x$SLEPC_DIR" = x || test "x$petscversion" = x); then
      enableslepc=no
      AC_MSG_RESULT(<<< SLEPc disabled.  Please set your "\$SLEPC_DIR" environment variable correctly to enable SLEPc. >>>)
    else
      AC_CHECK_HEADER([$SLEPC_DIR/include/slepcversion.h],
                      [SLEPC_INCLUDE="-I$SLEPC_DIR/include"],
                      [
                      AC_MSG_RESULT(<<< Invalid "\$SLEPC_DIR" detected (slepcversion.h not found). SLEPc disabled. >>>)
                      unset SLEPC_DIR
                      enableslepc=no
                      ])

      # I didn't check this code branch since I don't use a
      # PETSC_ARCH, but running AC_CHECK_HEADER on slepcconf.h
      # should be pretty safe... from what I can tell it just
      # #defines a few things.
      if (test "x$PETSC_ARCH" != "x"); then
        AC_CHECK_HEADER([$SLEPC_DIR/$PETSC_ARCH/include/slepcconf.h],
                        [SLEPC_INCLUDE="$SLEPC_INCLUDE -I$SLEPC_DIR/$PETSC_ARCH/include"])
      fi
    fi

    if (test "x$enableslepc" = "xyes") ; then

      # Similar to petsc, we need the slepc version number.
      # Note slepc will most likely only work with the corresponding version of petsc
      slepcmajor=`grep "define SLEPC_VERSION_MAJOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MAJOR[ ]*//g"`
      slepcminor=`grep "define SLEPC_VERSION_MINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MINOR[ ]*//g"`
      slepcsubminor=`grep "define SLEPC_VERSION_SUBMINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_SUBMINOR[ ]*//g"`
      slepcversion=$slepcmajor.$slepcminor.$slepcsubminor

      # Older PETSc (3.3 and earlier) might not work if the version numbers don't match exactly.
      oldpetsc=no

      if (test $petscmajor -lt 3) ; then
        oldpetsc=yes
      fi

      if (test $petscmajor == 3 -a $petscminor -lt 4) ; then
        oldpetsc=yes
      fi

      # If PETSc is old, warn if the full version numbers (including the subminor version number) don't match exactly.
      if (test "$oldpetsc" = "yes" -a $slepcversion != $petscversion) ; then
        AC_MSG_RESULT([<<< WARNING: PETSc version $petscversion does not match SLEPc version $slepcversion >>>])
      fi

      # For any PETSc version, warn if the major and minor version numbers don't match.
      if (test $slepcmajor != $petscmajor -o $slepcminor != $petscminor) ; then
        AC_MSG_RESULT([<<< WARNING: PETSc version $petscmajor.$petscminor does not match SLEPc version $slepcmajor.$slepcminor >>>])
      fi

      # OK, now we will create a temporary makefile to query SLEPc libs
      includefile=""

      # These files are not headers, so they will never pass an AC_CHECK_HEADER test.
      # If AC_CHECK_FILE is forbidden while cross-compiling, perhaps we can just fall back on
      # 'if test -r' calls?
      if (test -r ${SLEPC_DIR}/bmake/${PETSC_ARCH}/slepcconf) ; then
         includefile=${SLEPC_DIR}/bmake/${PETSC_ARCH}/slepcconf
      elif (test -r ${SLEPC_DIR}/conf/slepc_common_variables) ; then
         includefile=${SLEPC_DIR}/conf/slepc_common_variables
      elif (test -r ${SLEPC_DIR}/conf/slepc_variables) ; then
         includefile=${SLEPC_DIR}/conf/slepc_variables
      # 3.6.0
      elif (test -r ${SLEPC_DIR}/lib/slepc/conf/slepc_variables) ; then
        includefile=${SLEPC_DIR}/lib/slepc/conf/slepc_variables
      else
         AC_MSG_RESULT([<<< SLEPc configuration failed.  Could not find slepcconf/slepc_variables file. >>>])
         enableslepc=no
      fi

      if (test "x$enableslepc" = "xyes" -a "x$includefile" != "x"); then
        AC_MSG_RESULT(<<< Querying SLEPc configuration from $includefile >>>)
        printf '%s\n' "include $includefile" > Makefile_config_slepc
        printf '%s\n' "getSLEPC_LIBS:" >> Makefile_config_slepc
        printf '\t%s\n' "echo \$(SLEPC_LIB) \$(ARPACK_LIB)" >> Makefile_config_slepc

        SLEPC_LIBS=`make -s -f Makefile_config_slepc getSLEPC_LIBS`

        # If calling make worked, try to compile a trivial SLEPc
        # program to check our configuration...
        if (test x$? = x0); then
          AC_MSG_CHECKING(whether we can compile a trivial SLEPc program)
          AC_LANG_PUSH([C++])

          # Save the original CXXFLAGS contents
          saveCXXFLAGS="$CXXFLAGS"

          # Append both PETSc and SLEPc include paths to the CXXFLAGS variable.
          CXXFLAGS="$saveCXXFLAGS $PETSCINCLUDEDIRS $SLEPC_INCLUDE"

          AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
          @%:@include <slepc.h>
          static char help[]="";

          int main(int argc, char **argv)
          {
            SlepcInitialize(&argc, &argv, (char*)0, help);
            SlepcFinalize();
            return 0;
          }
          ]])],[
            AC_MSG_RESULT(yes)
            AC_MSG_RESULT(<<< Configuring library with SLEPc version $slepcversion support >>>)

            AC_DEFINE(HAVE_SLEPC, 1, [Flag indicating whether or not SLEPc is available])
            AC_SUBST(SLEPC_INCLUDE)
            AC_SUBST(SLEPC_LIBS)
          ],[
            AC_MSG_RESULT(no)
            AC_MSG_RESULT(<<< Compiling trivial SLEPc program failed. SLEPc disabled. >>>)
            enableslepc=no
          ])

          # Return CXXFLAGS and LANG to their original state.
          CXXFLAGS="$saveCXXFLAGS"
          AC_LANG_POP([C++])

        else
          enableslepc=no
          AC_MSG_RESULT(<<< SLEPc configuration query failed. SLEPc disabled. >>>)
        fi

        # rm temporary Makefile
        rm -f Makefile_config_slepc
      fi
    fi
  fi
])
