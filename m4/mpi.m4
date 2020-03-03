dnl ---------------------------------------------------------------------------
dnl check for the required MPI library
dnl ---------------------------------------------------------------------------
AC_DEFUN([ACX_MPI], [
dnl We need to verify that we've done AC_ARG_ENABLE(mpi) and AC_ARG_WITH(mpi)
dnl which occur in LIBMESH_COMPILER_CONTROL_ARGS
AS_IF([test x"CALLED_COMPILER_CONTROL_ARGS" != x1],
      [AC_MSG_ERROR([the compiler control args macro must be called before the mpi macro])])


AS_IF(
  [test x"$MPI_USING_WRAPPERS" = x1],
  [
    AC_LANG_PUSH([C++])
    AC_TRY_LINK([@%:@include <mpi.h>],
                [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                [
                  enablempi=yes
                  MPI_IMPL="mpi-wrapper-built-in"
                  AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                ],
                [
                  AC_MSG_RESULT([MPI wrapper $CXX cannot compile an MPI program. Disabling MPI...])
                  enablempi=no
                ])
    AC_LANG_POP([C++])
  ],
  [test -n "$MPI"],
  [
    MPI_LIBS_PATH="$MPI/lib"
    MPI_INCLUDES_PATH="$MPI/include"

    MPI_LIBS_TO_TEST=""
    AS_IF([test -e $MPI_LIBS_PATH/libmpi.a], [MPI_LIBS_TO_TEST="$MPI_LIBS_PATH/libmpi.a $MPI_LIBS_TO_TEST"])
    AS_IF([test -e $MPI_LIBS_PATH/libmpi.so], [MPI_LIBS_TO_TEST="$MPI_LIBS_PATH/libmpi.so $MPI_LIBS_TO_TEST"])
    AS_IF([test -e $MPI_LIBS_PATH/libmpi.dylib], [MPI_LIBS_TO_TEST="$MPI_LIBS_PATH/libmpi.dylib $MPI_LIBS_TO_TEST"])


    dnl look for LAM or other MPI implementation
    AS_IF([test x"$MPI_LIBS_TO_TEST" != x],
          [
            AS_ECHO(["note: testing $MPI_LIBS_PATH/libmpi(.a/.so)"])

            dnl Ensure the compiler finds the library...
            tmpLIBS=$LIBS
            AC_LANG_SAVE
            AC_LANG_CPLUSPLUS

            LIBS="-L$MPI_LIBS_PATH $LIBS"

            dnl look for lam_version_show in liblam.(a/so)
            dnl (this is needed in addition to libmpi.(a/so) for
            dnl LAM MPI
            AC_CHECK_LIB([lam],
                         [lam_show_version],
                         [
                           LIBS="-llam $LIBS"
                           MPI_LIBS="-llam $MPI_LIBS"
                         ],
                         [])

            dnl Quadrics MPI requires the elan library to be included too
            for lib in $MPI_LIBS_TO_TEST
            do
              if (nm $lib | grep elan > /dev/null); then
                AS_ECHO(["note: MPI found to use Quadrics switch, looking for elan library"])
                AC_CHECK_LIB([elan],
                             [elan_init],
                             [
                               LIBS="-lelan $LIBS"
                               MPI_LIBS="-lelan $MPI_LIBS"
                             ],
                             [AC_MSG_ERROR([Could not find elan library... exiting])])
                break
              fi
            done

            AC_CHECK_LIB([mpi],
                         [MPI_Init],
                         [
                           MPI_LIBS="-lmpi $MPI_LIBS"
                           MPI_LDFLAGS="-L$MPI_LIBS_PATH"
                           MPI_IMPL="mpi"
                           AC_MSG_RESULT([Found valid MPI installation...])
                         ],
                         [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])

            AC_LANG_RESTORE
            LIBS=$tmpLIBS
          ])

    dnl have we not found an implementation yet?
    AS_IF([test x"$MPI_IMPL" = x],
          [
            AS_IF([test -e $MPI_LIBS_PATH/libmpich.a || test -e $MPI_LIBS_PATH/libmpich.so],
                  [
                    AS_ECHO(["note: using $MPI_LIBS_PATH/libmpich(.a/.so)"])

                    dnl Ensure the compiler finds the library...
                    tmpLIBS=$LIBS
                    AC_LANG_SAVE
                    AC_LANG_CPLUSPLUS
                    LIBS="-L$MPI_LIBS_PATH $LIBS"

                    dnl Myricomm MPICH requires the gm library to be included too
                    if (nm $MPI_LIBS_PATH/libmpich.* | grep gm_open > /dev/null); then
                      AS_ECHO(["note: MPICH found to use Myricomm's Myrinet, looking for gm library"])

                            AS_IF([test "x$GMHOME" = "x"], [GMHOME="/usr"])
                            AC_ARG_WITH([gm],
                                        AS_HELP_STRING([--with-gm=PATH],
                                                       [Prefix where GM is installed (GMHOME)]),
                                        [GM="$withval"],
                                        [
                                          AS_ECHO(["note: GM library path not given... trying prefix=$MPIHOME"])
                                          GM=$GMHOME
                                        ])

                            LIBS="-L$GM/lib $LIBS"
                            MPI_LIBS="-L$GM/lib $MPI_LIBS"

                            AC_CHECK_LIB([gm],
                                         [gm_open],
                                         [
                                           LIBS="$LIBS -lgm"
                                           MPI_LIBS="$MPI_LIBS -lgm"
                                         ],
                                         [AC_MSG_ERROR( [Could not find gm library... exiting])])
                    fi

                    dnl look for MPI_Init in libmich.(a/so)
                    dnl try adding libmpl if we see it there; some MPICH2 versions
                    dnl require it.
                    AS_IF([test -e $MPI_LIBS_PATH/libmpl.a || test -e $MPI_LIBS_PATH/libmpl.so],
                          [
                            LIBS="-L$MPI_LIBS_PATH -lmpl $tmpLIBS"
                            AC_CHECK_LIB([mpich],
                                         [MPI_Init],
                                         [
                                           MPI_LIBS="-lmpich -lmpl $MPI_LIBS"
                                           MPI_LDFLAGS="-L$MPI_LIBS_PATH"
                                           MPI_IMPL="mpich"
                                           AC_MSG_RESULT([Found valid MPICH installation with libmpl...])
                                         ],
                                         [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])
                          ],
                          [
                            AC_CHECK_LIB([mpich],
                                         [MPI_Init],
                                         [
                                           MPI_LIBS="-lmpich $MPI_LIBS"
                                           MPI_LDFLAGS="-L$MPI_LIBS_PATH"
                                           MPI_IMPL="mpich"
                                           AC_MSG_RESULT([Found valid MPICH installation...])
                                         ],
                                         [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])
                          ])

                    AC_LANG_RESTORE
                    LIBS=$tmpLIBS
                  ])
          ])

    dnl ok we've done all the library work. Now let's check the header
    AS_IF([test "x$MPI_IMPL" != x],
          [
            dnl Ensure the compiler finds the header file...
            AS_IF([test -e $MPI_INCLUDES_PATH/mpi.h],
                  [
                    AS_ECHO(["note: using $MPI_INCLUDES_PATH/mpi.h"])
                    tmpCPPFLAGS=$CPPFLAGS
                    AC_LANG_SAVE
                    AC_LANG_CPLUSPLUS
                    CPPFLAGS="-I$MPI_INCLUDES_PATH $CPPFLAGS"
                    AC_CHECK_HEADER([mpi.h],
                                    [
                                      AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                                      enablempi=yes
                                    ],
                                    [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no])
                    MPI_INCLUDES="-I$MPI_INCLUDES_PATH"
                    AC_LANG_RESTORE
                    CPPFLAGS=$tmpCPPFLAGS
                  ],
                  [test -e $MPI_INCLUDES_PATH/mpi/mpi.h],
                  [
                    MPI_INCLUDES_PATH=$MPI_INCLUDES_PATH/mpi
                    AS_ECHO(["note: using $MPI_INCLUDES_PATH/mpi.h"])
                    tmpCPPFLAGS=$CPPFLAGS
                    AC_LANG_SAVE
                    AC_LANG_CPLUSPLUS
                    CPPFLAGS="-I$MPI_INCLUDES_PATH $CPPFLAGS"
                    AC_CHECK_HEADER([mpi.h],
                                    [
                                      AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                                      enablempi=yes
                                    ],
                                    [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no] )
                    MPI_INCLUDES="-I$MPI_INCLUDES_PATH"
                    AC_LANG_RESTORE
                    CPPFLAGS=$tmpCPPFLAGS
                  ],
                  [
                    AC_MSG_RESULT([Could not find MPI header <mpi.h>...])
                    enablempi=no
                  ])
          ],
          [
            dnl no MPI install found, see if the compiler "natively" supports it by
            dnl attempting to link a test application without any special flags.
            AC_LANG_SAVE
            AC_LANG_CPLUSPLUS
            AC_TRY_LINK([@%:@include <mpi.h>],
                        [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                        [
                           MPI_IMPL="CXX-built-in"
                           AC_MSG_RESULT( [$CXX Compiler Supports MPI] )
                           AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                           enablempi=yes
                        ],
                        [AC_MSG_RESULT([$CXX Compiler Does NOT Support MPI...]); enablempi=no])
            AC_LANG_RESTORE
          ])
  ],
  dnl We weren't supplied --with-mpi and no MPI environment variables were set
  dnl Let's check to see if we got anything from PETSc
  [test x"$PETSC_HAVE_MPI" = x1],
  [
    tmpCPPFLAGS=$CPPFLAGS
    tmpLIBS=$LIBS
    AC_LANG_PUSH([C++])
    LIBS="$PETSC_MPI_LINK_LIBS $LIBS"
    CPPFLAGS="$PETSC_MPI_INCLUDE_DIRS $CPPFLAGS"
    AC_TRY_LINK([@%:@include <mpi.h>],
                [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                [
                  enablempi=yes
                  MPI_IMPL="petsc"
                  AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                  AS_IF([test x"$PETSC_MPI_INCLUDE_DIRS" != x], [MPI_INCLUDES="$PETSC_MPI_INCLUDE_DIRS"])
                  AS_IF([test x"$PETSC_MPI_LINK_LIBS" != x], [MPI_LIBS="$PETSC_MPI_LINK_LIBS"])
                ],
                [
                  AC_MSG_RESULT([PETSc MPI cannot compile an MPI program. Disabling MPI...])
                  enablempi=no
                ])
    AC_LANG_POP([C++])
    LIBS=$tmpCPPFLAGS
    CPPFLAGS=$tmpLIBS
  ],
  dnl We weren't supplied --with-mpi and no MPI environment variables were set. Also we didn't
  dnl get any MPI info from a PETSc install. Let's check to see if CXX natively supports MPI
  [
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    AC_TRY_LINK([@%:@include <mpi.h>],
                [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                [
                   MPI_IMPL="CXX-built-in"
                   AC_MSG_RESULT( [$CXX Compiler Supports MPI] )
                   AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                   enablempi=yes
                ],
                [AC_MSG_RESULT([$CXX Compiler Does NOT Support MPI...]); enablempi=no])
    AC_LANG_RESTORE
  ])

  dnl Save variables...
  AC_SUBST(MPI_IMPL)
  AC_SUBST(MPI_LIBS)
  AC_SUBST(MPI_LDFLAGS)
  AC_SUBST(MPI_INCLUDES)
])
