#---------------------------------------------------------
# Scrape configuration information from PETSc, but don't
# try to verify the installation yet. We need to establish
# a compiler first
#---------------------------------------------------------
AC_DEFUN([SCRAPE_PETSC_CONFIGURE],
[
  dnl We need to verify that we've done AC_ARG_ENABLE(petsc)
  dnl which occurs in COMPILER_CONTROL_ARGS
  AC_REQUIRE([COMPILER_CONTROL_ARGS])

  # Trump --enable-petsc with --disable-mpi
  AS_IF([test "x$enablempi" = xno],
        [enablepetsc=no;enablepetsc_mpi=no])

  AC_ARG_VAR([PETSC_DIR],  [path to PETSc installation])
  AC_ARG_VAR([PETSC_ARCH], [PETSc build architecture])

  AS_IF([test "$enablepetsc" !=  no],
        [
          # If the user doesn't have any PETSC directory specified, let's check to
          # see if it's installed via Ubuntu module
          AS_IF([test "x${PETSC_DIR}" = x],
                [
                  AC_PATH_PROG(PETSCARCH, petscarch)
                  AS_IF([test "x$PETSCARCH" != x],
                        [
                          export PETSC_DIR=/usr/lib/petsc
                          export PETSC_ARCH=`$PETSCARCH`
                          AS_IF([test -d ${PETSC_DIR}],
                                [
                                  AC_MSG_RESULT([using system-provided PETSC_DIR ${PETSC_DIR}])
                                  AC_MSG_RESULT([using system-provided PETSC_ARCH ${PETSC_ARCH}])
                                ])
                        ])
                ])

          AS_IF([test x"$PETSC_DIR" = x], [enablepetsc=no;enablepetsc_mpi=no])
        ])

  AS_IF([test "$enablepetsc" != no],
        [
          dnl Check for snoopable MPI
          AS_IF([test -r ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf], dnl 2.3.x
                [PETSC_MPI=`grep MPIEXEC ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf | grep -v mpiexec.uni`],
                [test -r ${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables], dnl 3.0.x
                [PETSC_MPI=`grep MPIEXEC ${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables | grep -v mpiexec.uni`],
                [test -r ${PETSC_DIR}/conf/petscvariables], dnl 3.0.x
                [PETSC_MPI=`grep MPIEXEC ${PETSC_DIR}/conf/petscvariables | grep -v mpiexec.uni`],
                [test -r ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables], dnl 3.6.x
                [PETSC_MPI=`grep MPIEXEC ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables | grep -v mpiexec.uni`],
                [test -r ${PETSC_DIR}/lib/petsc/conf/petscvariables], dnl 3.6.x
                [PETSC_MPI=`grep MPIEXEC ${PETSC_DIR}/lib/petsc/conf/petscvariables | grep -v mpiexec.uni`])

          dnl If we couldn't snoop MPI from PETSc, fall back on ACX_MPI.
          AS_IF([test "x$PETSC_MPI" != x],
                [
                ],
                [
                  enablepetsc_mpi=no
                ])

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

                dnl Ubuntu PETSc dpkg
                [test -r ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/variables],
                [PREFIX_INSTALLED_PETSC=yes
                 PETSC_VARS_FILE=${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/variables],

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
                [
                  AC_MSG_RESULT([<<< Could not find a viable PETSc Makefile to determine PETSC_CC_INCLUDES, etc. >>>])
                   enablepetsc=no; enablepetsc_mpi=no
                ]) dnl determining PETSC_VARS_FILE
        ]) dnl AS_IF(enablepetsc)

  dnl If we haven't been disabled yet...
  AS_IF([test "$enablepetsc" != no],
        [
          dnl Set some include and link variables by building and running temporary Makefiles.
          AS_IF([test "$PREFIX_INSTALLED_PETSC" = "no"],
                [
                  PETSCLINKLIBS=`make -s -C ${PETSC_DIR} getlinklibs`
                  PETSCINCLUDEDIRS=`make -s -C ${PETSC_DIR} getincludedirs`
                  PETSC_CXX=`make -s -C $PETSC_DIR getcxxcompiler`
                  PETSC_MPI_INCLUDE_DIRS=`make -s -C $PETSC_DIR getmpiincludedirs`
                  PETSC_MPI_LINK_LIBS=`make -s -C $PETSC_DIR getmpilinklibs`
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
                  printf '\t%s' "echo " >> Makefile_config_petsc
                  AS_IF([test -d ${PETSC_DIR}/include],
                        [printf '%s ' "-I\$(PETSC_DIR)/include" >> Makefile_config_petsc])
                  AS_IF([test -d ${PETSC_DIR}/${PETSC_ARCH}/include],
                        [printf '%s ' "-I\$(PETSC_DIR)/\$(PETSC_ARCH)/include" >> Makefile_config_petsc])
                  printf '%s\n' "\$(BLOCKSOLVE_INCLUDE) \$(HYPRE_INCLUDE) \$(PACKAGES_INCLUDES)" >> Makefile_config_petsc
                  printf '%s\n'  "getPETSC_CC_INCLUDES:" >> Makefile_config_petsc
                  printf '\t%s\n' "echo \$(PETSC_CC_INCLUDES)" >> Makefile_config_petsc
                  printf '%s\n' "getPETSC_FC_INCLUDES:" >> Makefile_config_petsc
                  printf '\t%s\n' "echo \$(PETSC_FC_INCLUDES)" >> Makefile_config_petsc
                  printf '%s\n' "getlinklibs:" >> Makefile_config_petsc
                  printf '\t%s\n' "echo \$(PETSC_SNES_LIB)" >> Makefile_config_petsc
                  printf '%s\n' "getcxxcompiler:" >> Makefile_config_petsc
                  printf '\t%s\n' "echo \$(CXX)" >> Makefile_config_petsc
                  printf '%s\n' "getmpiincludedirs:" >> Makefile_config_petsc
                  printf '\t%s\n' "echo \$(MPI_INCLUDE)" >> Makefile_config_petsc
                  printf '%s\n' "getmpilinklibs:" >> Makefile_config_petsc
                  printf '\t%s\n' "echo \$(MPI_LIB)" >> Makefile_config_petsc
                  PETSCLINKLIBS=`make -s -f Makefile_config_petsc getlinklibs`
                  PETSCINCLUDEDIRS=`make -s -f Makefile_config_petsc getincludedirs`
                  PETSC_CC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_CC_INCLUDES`
                  PETSC_FC_INCLUDES=`make -s -f Makefile_config_petsc getPETSC_FC_INCLUDES`
                  PETSC_CXX=`make -s -f Makefile_config_petsc getcxxcompiler`
                  PETSC_MPI_INCLUDE_DIRS=`make -s -f Makefile_config_petsc getmpiincludedirs`
                  PETSC_MPI_LINK_LIBS=`make -s -f Makefile_config_petsc getmpilinklibs`
                  rm -f Makefile_config_petsc
                ]) dnl scrape petsc cxx, includes, and libs

          dnl Debugging: see what actually got set for PETSCINCLUDEDIRS
          dnl echo ""
          dnl echo "PETSCLINKLIBS=$PETSCLINKLIBS"
          dnl echo "PETSCINCLUDEDIRS=$PETSCINCLUDEDIRS"
          dnl echo ""

          # We sometimes need the full CC_INCLUDES to access a
          # PETSc-snooped MPI
          PETSCINCLUDEDIRS="$PETSCINCLUDEDIRS $PETSC_CC_INCLUDES"

    ]) dnl AS_IF(enable_petsc)

  AS_IF([test "$enablepetsc" = no && test "$enablepetsc_mpi" != no],
        [
          AC_MSG_ERROR([petsc was disabled but petscs mpi was not disabled])
          AC_MSG_ERROR([something wrong must have happened during the configure process])
          AC_MSG_ERROR([please contact the libmesh-users mailing list for support])
        ])

  AS_IF([test "$enablepetsc_mpi" != no], [PETSC_HAVE_MPI=1])
])
