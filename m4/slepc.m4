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
          dnl If SLEPC_DIR is not set, but $PETSC_DIR is, it's possible
          dnl that there is a SLEPc within PETSc.  Not sure about using
          dnl "export" here but the same thing is done in petsc.m4 so it
          dnl probably doesn't hurt anything.
          AS_IF([test "x$SLEPC_DIR" = "x" && test "x$PETSC_DIR" != "x"],
                [export SLEPC_DIR=$PETSC_DIR])

          # It doesn't seem likely that we would find a /usr/lib/slepc
          # when we already didn't find a /usr/lib/petsc (otherwise,
          # PETSC_DIR would be set!)  so let's rule that case out just
          # to be safe.
          AS_IF([test "x$SLEPC_DIR" = "x"],
                [
                  dnl SLEPC_DIR (and hence PETSC_DIR) was not set. Even *if* we found a SLEPc
                  dnl in /usr/lib/slepc, I don't see how we would be able to configure it reliably
                  dnl without PETSC_DIR? This seems like it shouldn't happen, but we'll print a
                  dnl message just in case so people know what's going on.
                  enableslepc=no
                  AC_MSG_RESULT([<<< Could not infer SLEPC_DIR from PETSC_DIR, SLEPc disabled. >>>])
                ])

          dnl Check for different flavors of PETSc installs.
          AS_IF([test "x$SLEPC_DIR" != "x"],
                [
                  dnl 1.) A standalone (make install) SLEPc build
                  dnl 2.) SLEPc in a PETSC_ARCH directory
                  dnl 3.) SLEPc in /usr/lib/slepc
                  slepc_standalone=no
                  slepc_in_petsc_arch=no
                  slepc_in_usr_lib=no

                  dnl 1.)
                  AC_CHECK_HEADER([$SLEPC_DIR/include/slepcversion.h],
                                  [SLEPC_INCLUDE="-I$SLEPC_DIR/include"
                                   slepc_standalone=yes],[])

                  dnl 2.)
                  AS_IF([test "x$slepc_standalone" != "xyes" && test "x$PETSC_ARCH" != "x"],
                        [
                          AC_CHECK_HEADER([$SLEPC_DIR/$PETSC_ARCH/include/slepcversion.h],
                                          [SLEPC_INCLUDE="-I$SLEPC_DIR/$PETSC_ARCH/include"
                                          slepc_in_petsc_arch=yes],[])
                        ])

                  dnl 3.)
                  AS_IF([test "x$slepc_standalone" != "xyes" && test "x$slepc_in_petsc_arch" != "xyes"],
                        [
                          AC_CHECK_HEADER([/usr/lib/slepc/slepcversion.h],
                                          [
                                            export SLEPC_DIR=/usr/lib/slepc
                                            SLEPC_INCLUDE="-I$SLEPC_DIR/include"
                                          ],
                                          [
                                            AC_MSG_RESULT([<<< No valid SLEPC installation found, SLEPc disabled. >>>])
                                            unset SLEPC_DIR
                                            enableslepc=no
                                          ])
                        ])
                ])
        ])

    AS_IF([test "x$enableslepc" = "xyes"],
          [
            dnl Similar to petsc, we need the slepc version number.
            dnl Note slepc will most likely only work with the corresponding version of petsc
            slepc_version_header=$SLEPC_DIR/$PETSC_ARCH/include/slepcversion.h
            AS_IF([test "x$slepc_in_petsc_arch" != "xyes"],
                  [slepc_version_header=$SLEPC_DIR/include/slepcversion.h])

            slepcmajor=`grep "define SLEPC_VERSION_MAJOR" $slepc_version_header | sed -e "s/#define SLEPC_VERSION_MAJOR[ ]*//g"`
            slepcminor=`grep "define SLEPC_VERSION_MINOR" $slepc_version_header | sed -e "s/#define SLEPC_VERSION_MINOR[ ]*//g"`
            slepcsubminor=`grep "define SLEPC_VERSION_SUBMINOR" $slepc_version_header | sed -e "s/#define SLEPC_VERSION_SUBMINOR[ ]*//g"`

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
                  [test -r ${SLEPC_DIR}/${PETSC_ARCH}/lib/slepc/conf/slepc_variables], [includefile=${SLEPC_DIR}/${PETSC_ARCH}/lib/slepc/conf/slepc_variables],
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
