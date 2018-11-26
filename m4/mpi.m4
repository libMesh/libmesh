dnl ---------------------------------------------------------------------------
dnl check for the required MPI library
dnl ---------------------------------------------------------------------------
AC_DEFUN([ACX_MPI], [

dnl if MPIHOME is empty, set it to an invalid directory and try to find the correct location.
AS_IF([test "x$MPIHOME" = "x"], [MPIHOME="/.."])

AC_ARG_WITH([mpi],
            AS_HELP_STRING([--with-mpi=PATH],
                           [Prefix where MPI is installed (MPIHOME)]),
            [MPI="$withval"],
            [
              AS_ECHO(["note: MPI library path not given... trying to find it myself"])
              MPI=$MPIHOME
            ])

AS_IF([test -z "$MPI"], [MPI="/.."])

dnl if MPI is not set in the config, search it via the $INCLUDE and $LD_LIBRARY_PATH variables;
dnl     with this we can hope to get a configure that is consistent with the environment.
AS_IF([test "$MPI" = "/.."],
 [
    for i in $(echo $LD_LIBRARY_PATH | tr ':' '\n')
    do
       AS_IF([ test -e $i/libmpi.a || test -e $MPI_LIBS_PATH/libmpi.so ],
       [MPI=$(dirname $i)
       break])
    done
 ])
AS_IF([test "$MPI" = "/.."],
 [
    for i in $(echo $INCLUDE | tr ':' '\n')
    do
       AS_IF([ test -e $i/mpi.h ],
       [MPI=$(dirname $i)
       break ])
    done
 ])

dnl if we don't find it, default to standard-values and proceed as before:
AS_IF([test "$MPI" = "/.."], [MPI="/usr"])

MPI_LIBS_PATH="$MPI/lib"
MPI_INCLUDES_PATH="$MPI/include"

dnl Check that the compiler uses the library we specified...

dnl look for LAM or other MPI implementation
AS_IF([test -e $MPI_LIBS_PATH/libmpi.a || test -e $MPI_LIBS_PATH/libmpi.so],
      [
        AS_ECHO(["note: using $MPI_LIBS_PATH/libmpi(.a/.so)"])

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
        if (nm $MPI_LIBS_PATH/libmpi.* | grep elan > /dev/null); then
          AS_ECHO(["note: MPI found to use Quadrics switch, looking for elan library"])
          AC_CHECK_LIB([elan],
                       [elan_init],
                       [
                         LIBS="-lelan $LIBS"
                         MPI_LIBS="-lelan $MPI_LIBS"
                       ],
                       [AC_MSG_ERROR([Could not find elan library... exiting])])
        fi

        AC_CHECK_LIB([mpi],
                     [MPI_Init],
                     [
                       MPI_LIBS="-lmpi $MPI_LIBS"
                       MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
                       MPI_IMPL="mpi"
                       AC_MSG_RESULT([Found valid MPI installation...])
                     ],
                     [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])

        AC_LANG_RESTORE
        LIBS=$tmpLIBS
      ])

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
                               MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
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
                               MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
                               MPI_IMPL="mpich"
                               AC_MSG_RESULT([Found valid MPICH installation...])
                             ],
                             [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])
              ])

        AC_LANG_RESTORE
        LIBS=$tmpLIBS
      ])

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
                                [AC_MSG_RESULT([MPI-headers are working as expected.])],
                                [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no])
                MPI_INCLUDES_PATHS="-I$MPI_INCLUDES_PATH"
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
                                [AC_MSG_RESULT([MPI-headers are working as expected.])],
                                [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no] ) MPI_INCLUDES_PATHS="-I$MPI_INCLUDES_PATH"
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
        AC_TRY_LINK([@%:@include <mpi.h>],
                    [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                    [
                       MPI_IMPL="built-in"
                       AC_MSG_RESULT( [$CXX Compiler Supports MPI] )
                    ],
                    [AC_MSG_RESULT([$CXX Compiler Does NOT Support MPI...]); enablempi=no])
      ])
     dnl Save variables...
     AC_SUBST(MPI)
     AC_SUBST(MPI_IMPL)
     AC_SUBST(MPI_LIBS)
     AC_SUBST(MPI_LIBS_PATH)
     AC_SUBST(MPI_LIBS_PATHS)
     AC_SUBST(MPI_INCLUDES_PATH)
     AC_SUBST(MPI_INCLUDES_PATHS)

     VERSION_MPI
])


AC_DEFUN([VERSION_MPI], [

  AC_MSG_RESULT([Checking Version of MPI now:])
  AS_IF([test "x$enablempi" != xno],
  [
        dnl check that MPI_VERSION and MPI_SUBVERSION are defined.
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        dnl first catch undefined (evaluated to 0) and 1.X cases
        AC_TRY_LINK([@%:@include <mpi.h>],
                    [
                     @%:@if MPI_VERSION < 2
                     @%:@error "MPI-version 1.X is not supported."
                     @%:@endif
                    ],
                    [ dnl true (MPI_VERSION >= 2): Check which version it is:
                     AC_TRY_LINK([@%:@include <mpi.h>],
                                  [@%:@if MPI_VERSION > 4
                                  @%:@error "MPI-version too high. There is some error."
                                  @%:@endif
                                  ],
                                  [
                                   dnl if MPI_VERSION is 2, through a warning
                                   AC_TRY_LINK([@%:@include <mpi.h>],
                                               [@%:@if MPI_VERSION == 2
                                               @%:@error "deprecated MPI-version"
                                               @%:@endif
                                               ],
                                               [
                                                 AC_MSG_RESULT([ MPI was found with version >= 3.0. ]);
                                               ],
                                               [
                                                 AC_MSG_WARN([MPI-version is 2.X. This feature is deprecated. ]);
                                               ])
                                  dnl set necessary variables etc: Here, all enablempi-cases are in one place.
                                  AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                                  ],
                                  [
                                  AC_MSG_WARN(["ERROR: MPI-version seems to be unreasonably high. Please check your library. Disable MPI now..."]); enablempi=no
                                  ])
                    ],
                    [ dnl MPI_VERSION is not defined or has unexpected value
                    dnl Lets look more closely to be able to tell the user what is actually wrong.
                    dnl First, test if it defined at all
                    AC_TRY_LINK([@%:@include <mpi.h>],
                                [
                                 @%:@if MPI_VERSION == 0
                                 @%:@error "MPI-version 1.X is not supported."
                                 @%:@endif
                                 ],
                                 [
                                 dnl test if it defined at all
                                 AC_TRY_LINK([@%:@include <mpi.h>],
                                             [
                                              @%:@if MPI_VERSION == 1
                                              @%:@error "MPI-version 1.X is not supported."
                                              @%:@endif
                                             ],
                                             [
                                              AC_MSG_WARN(["ERROR: MPI-version is 1.X. Only MPI 2.X or compatible is supported. Disable MPI now..."])
                                             ],
                                             [ dnl MPI_VERSION is < 0 or not an integer. This should make you worry!
                                              AC_MSG_WARN(["ERROR: The MPI_VERSION-string has a weird value. Please check your library. Disable MPI now..."])
                                             ])
                                 ],
                                 [ dnl MPI_VERSION is 0.
                                  AC_MSG_WARN(["ERROR: MPI_VERSION is not defined: You seem to have MPI< 1.5 but need MPI 2.X or compatible. Disable MPI now..."])
                                 ])
                     dnl in any case, we need to disable MPI now.
                     enablempi=no
                    ])
        dnl this is just for understanding the Civet-problem:
        AC_TRY_LINK([@%:@include <mpi.h>],
                    [
                     @%:@ifndef MPICH
                     @%:@error "MPICH is undefined"
                     @%:@endif
                     ],
                     [AC_MSG_WARN(["MPICH is defined."])],
                     [AC_MSG_WARN(["MPICH is not defined."])]),
        AC_TRY_LINK([@%:@include <mpi.h>],
                    [
                     @%:@if MPICH_NAME == 3
                     @%:@error "MPICH-version is 3."
                     @%:@endif
                     ],
                     [AC_MSG_WARN(["MPICH is not version 3."])],
                     [AC_MSG_WARN(["MPICH is version 3."])])
        AC_TRY_LINK([@%:@include <mpi.h>],
                    [
                     @%:@ifndef MPI_T_ERR_INVALID
                     @%:@error "this variable belongs to MPI 3.1 and should be supported."
                     @%:@endif
                     ],
                     [AC_MSG_WARN(["MPI 3 is supported."])],
                     [AC_MSG_WARN(["MPI 3 is not supported."])])
        AC_LANG_RESTORE
  ])
])
