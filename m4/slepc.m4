# -------------------------------------------------------------
# SLEPc
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SLEPC],
[
  AC_ARG_ENABLE(slepc,
                AS_HELP_STRING([--disable-slepc],
                               [build without SLEPc eigen solver support]),
                [AS_CASE("${enableval}",
                         [yes], [enableslepc=yes],
                         [no],  [enableslepc=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-slepc)])],
                [enableslepc=$enablepetsc]) # if unspecified, infer from what's up with PETSc

  AC_ARG_VAR([SLEPC_DIR], [path to SLEPc installation])

  AS_IF([test "$enableslepc" !=  no],
        [
          dnl Test to see if SLEPC_DIR set by user.  If not set, then
          dnl try to autodetect in a default directory
          AS_IF([test "x$SLEPC_DIR" = x],
                [export SLEPC_DIR=/usr/lib/slepc])

          dnl Test to see if SLEPC_DIR and petscversion were set by user or
          dnl autodetection.  If not set, then disable slepc, print a message.
          AS_IF([test "x$SLEPC_DIR" = "x" || test "x$petscversion" = "x"],
                [
                  enableslepc=no
                  AC_MSG_RESULT(<<< SLEPc disabled.  Please set your "\$SLEPC_DIR" environment variable correctly to enable SLEPc. >>>)
                ],
                [
                  AC_CHECK_HEADER([$SLEPC_DIR/include/slepcversion.h],
                                  [SLEPC_INCLUDE="-I$SLEPC_DIR/include"],
                                  [
                                    AC_MSG_RESULT(<<< Invalid "\$SLEPC_DIR" detected (slepcversion.h not found). SLEPc disabled. >>>)
                                    unset SLEPC_DIR
                                    enableslepc=no
                                  ])
                ])

            dnl I didn't check this code branch since I don't use a
            dnl PETSC_ARCH, but running AC_CHECK_HEADER on slepcconf.h
            dnl should be pretty safe... from what I can tell it just
            dnl #defines a few things.
            AS_IF([test "x$PETSC_ARCH" != "x"],
                  [AC_CHECK_HEADER([$SLEPC_DIR/$PETSC_ARCH/include/slepcconf.h],
                                   [SLEPC_INCLUDE="$SLEPC_INCLUDE -I$SLEPC_DIR/$PETSC_ARCH/include"])])
        ])

    AS_IF([test "x$enableslepc" = "xyes"],
          [
            dnl Similar to petsc, we need the slepc version number.
            dnl Note slepc will most likely only work with the corresponding version of petsc
            slepcmajor=`grep "define SLEPC_VERSION_MAJOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MAJOR[ ]*//g"`
            slepcminor=`grep "define SLEPC_VERSION_MINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MINOR[ ]*//g"`
            slepcsubminor=`grep "define SLEPC_VERSION_SUBMINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_SUBMINOR[ ]*//g"`
            slepcversion=$slepcmajor.$slepcminor.$slepcsubminor

            dnl Older PETSc (3.3 and earlier) might not work if the version numbers don't match exactly.
            oldpetsc=no

            AS_IF([test $petscmajor -lt 3], [oldpetsc=yes])
            AS_IF([test "$petscmajor" = "3" && test $petscminor -lt 4], [oldpetsc=yes])

            dnl If PETSc is old, warn if the full version numbers (including the subminor version number) don't match exactly.
            AS_IF([test "$oldpetsc" = "yes" && test "$slepcversion" != "$petscversion"],
                  [AC_MSG_RESULT([<<< WARNING: PETSc version $petscversion does not match SLEPc version $slepcversion >>>])])

            dnl For any PETSc version, warn if the major and minor version numbers don't match.
            AS_IF([test "$slepcmajor" != "$petscmajor" || test "$slepcminor" != "$petscminor"],
                  [AC_MSG_RESULT([<<< WARNING: PETSc version $petscmajor.$petscminor does not match SLEPc version $slepcmajor.$slepcminor >>>])])

            dnl OK, now we will create a temporary makefile to query SLEPc libs
            includefile=""

            dnl Figure out which file to include for SLEPc configuration variables.
            AS_IF([test -r ${SLEPC_DIR}/bmake/${PETSC_ARCH}/slepcconf],  [includefile=${SLEPC_DIR}/bmake/${PETSC_ARCH}/slepcconf],
                  [test -r ${SLEPC_DIR}/conf/slepc_common_variables],    [includefile=${SLEPC_DIR}/conf/slepc_common_variables],
                  [test -r ${SLEPC_DIR}/conf/slepc_variables],           [includefile=${SLEPC_DIR}/conf/slepc_variables],
                  [test -r ${SLEPC_DIR}/lib/slepc/conf/slepc_variables], [includefile=${SLEPC_DIR}/lib/slepc/conf/slepc_variables], dnl 3.6.0
                  [AC_MSG_RESULT([<<< SLEPc configuration failed.  Could not find slepcconf/slepc_variables file. >>>])
                   enableslepc=no])

            AS_IF([test "x$enableslepc" = "xyes" && test "x$includefile" != "x"],
                  [
                    AC_MSG_RESULT(<<< Querying SLEPc configuration from $includefile >>>)
                    printf '%s\n' "include $includefile" > Makefile_config_slepc
                    printf '%s\n' "getSLEPC_LIBS:" >> Makefile_config_slepc
                    printf '\t%s\n' "echo \$(SLEPC_LIB) \$(ARPACK_LIB)" >> Makefile_config_slepc

                    SLEPC_LIBS=`make -s -f Makefile_config_slepc getSLEPC_LIBS`

                    dnl If calling make worked, try to compile a trivial SLEPc
                    dnl program to check our configuration...
                    AS_IF([test x$? = x0],
                          [
                            AC_MSG_CHECKING(whether we can compile a trivial SLEPc program)
                            AC_LANG_PUSH([C++])

                            dnl Save the original CXXFLAGS contents
                            saveCXXFLAGS="$CXXFLAGS"

                            dnl Append both PETSc and SLEPc include paths to the CXXFLAGS variable.
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

                              AC_DEFINE_UNQUOTED(DETECTED_SLEPC_VERSION_MAJOR, [$slepcmajor],
                                [SLEPc's major version number, as detected by slepc.m4])

                              AC_DEFINE_UNQUOTED(DETECTED_SLEPC_VERSION_MINOR, [$slepcminor],
                                [SLEPc's minor version number, as detected by slepc.m4])

                              AC_DEFINE_UNQUOTED(DETECTED_SLEPC_VERSION_SUBMINOR, [$slepcsubminor],
                                [SLEPc's subminor version number, as detected by slepc.m4])
                            ],[
                              AC_MSG_RESULT(no)
                              AC_MSG_RESULT(<<< Compiling trivial SLEPc program failed. SLEPc disabled. >>>)
                              enableslepc=no
                            ])

                            dnl Return CXXFLAGS and LANG to their original state.
                            CXXFLAGS="$saveCXXFLAGS"
                            AC_LANG_POP([C++])
                          ],
                          [
                            enableslepc=no
                            AC_MSG_RESULT(<<< SLEPc configuration query failed. SLEPc disabled. >>>)
                          ])

                    dnl rm temporary Makefile
                    rm -f Makefile_config_slepc
                  ])
         ])
])
