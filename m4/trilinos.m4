dnl -------------------------------------------------------------
dnl Trilinos 10
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS_10],
[
  AS_IF([test "x$TRILINOS_DIR" = "x"],
        [
          dnl Ubuntu trilinos package?
          AS_IF([test -d /usr/include/trilinos],
                [TRILINOS_DIR=/usr/include/trilinos],
                [TRILINOS_DIR=no])
        ])

  AC_ARG_WITH(trilinos,
              AS_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  AS_IF([test "x$withtrilinosdir" != "xno"],
        [
          AS_IF([test -r $withtrilinosdir/include/Makefile.export.Trilinos],
                [TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.Trilinos],
                [test -r $withtrilinosdir/Makefile.export.Trilinos],
                [TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/Makefile.export.Trilinos],
                [enabletrilinos10=no])

          AS_IF([test "x$enabletrilinos10" != "xno"],
                [
                  dnl This actually implies Trilinos 10+, as this configure
                  dnl function also works with Trilinos 11.
                  enabletrilinos10=yes

                  dnl Try to detect the full Trilinos version string by grepping through Trilinos_version.h.
                  dnl The version string is surrounded by quotes, so we use a slight variation on our usual
                  dnl grep command, and pass the result through 'tr' to remove those.
                  trilinosversionstring=`grep "define TRILINOS_VERSION_STRING" $withtrilinosdir/include/Trilinos_version.h | sed -e "s/#define TRILINOS_VERSION_STRING[ ]*//g" | tr -d '"'`

                  dnl If we couldn't extract anything for the full version string, try to just get the major version number.
                  AS_IF([test "x$trilinosversionstring" = "x"],
                        [
                          trilinosversionstring=`grep "define TRILINOS_MAJOR_VERSION" $withtrilinosdir/include/Trilinos_version.h | sed -e "s/#define TRILINOS_MAJOR_VERSION[ ]*//g"`
                        ])

                  dnl If we couldn't extract the actual major version number, go with "10+"
                  AS_IF([test "x$trilinosversionstring" = "x"],
                        [trilinosversionstring=10+])

                  AC_DEFINE(HAVE_TRILINOS, 1,
                            [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
                  AC_MSG_RESULT([<<< Configuring library with Trilinos $trilinosversionstring support >>>])

                  dnl ------------------------------------------------------
                  dnl AztecOO
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/AztecOO_config.h],
                                  [enableaztecoo=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/AztecOO_config.h],
                                                   [enableaztecoo=yes],
                                                   [AC_CHECK_HEADER([$withtrilinosdir/packages/aztecoo/src/AztecOO_config.h],
                                                                    [enableaztecoo=yes],
                                                                    [enableaztecoo=no])
                                                                    ])])

                  AS_IF([test "x$enableaztecoo" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_AZTECOO, 1, [Flag indicating whether the library shall be compiled to use the Trilinos AztecOO linear solver])
                          AC_MSG_RESULT(<<< Configuring library with AztecOO support >>>)
                        ])

                  dnl ------------------------------------------------------
                  dnl NOX
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/NOX_Config.h],
                                  [enablenox=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/NOX_Config.h],
                                                   [enablenox=yes],
                                                   [AC_CHECK_HEADER([$withtrilinosdir/packages/nox/src/NOX_Config.h],
                                                                    [enablenox=yes],
                                                                    [enablenox=no])])])

                  AS_IF([test "x$enablenox" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_NOX, 1, [Flag indicating whether the library shall be compiled to use the Trilinos NOX nonlinear solver])
                          AC_MSG_RESULT(<<< Configuring library with NOX support >>>)
                        ])

                  dnl ------------------------------------------------------
                  dnl ML - prevent ML from keying on our 'HAVE_PETSC' by
                  dnl undefining it using the fourth argument to AC_CHECK_HEADER.
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/ml_include.h],
                                  [enableml=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/ml_include.h],
                                                   [enableml=yes],
                                                   [AC_CHECK_HEADER([$withtrilinosdir/packages/ml/src/ml_config.h],
                                                                    [enableml=yes],
                                                                    [enableml=no],
                                                                    [
                                                                    @%:@ifdef HAVE_PETSC
                                                                    @%:@undef HAVE_PETSC
                                                                    @%:@endif
                                                                    ])],
                                                   [
                                                   @%:@ifdef HAVE_PETSC
                                                   @%:@undef HAVE_PETSC
                                                   @%:@endif
                                                   ])],
                                  [
                                  @%:@ifdef HAVE_PETSC
                                  @%:@undef HAVE_PETSC
                                  @%:@endif
                                  ])

                  AS_IF([test "x$enableml" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_ML, 1, [Flag indicating whether the library shall be compiled to use the Trilinos ML package])
                          AC_MSG_RESULT(<<< Configuring library with ML support >>>)
                        ])

                  dnl ------------------------------------------------------
                  dnl TPetra
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/Tpetra_config.h],
                                  [enabletpetra=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/Tpetra_config.h],
                                                   [enabletpetra=yes],
                                                   [AC_CHECK_HEADER([$withtrilinosdir/packages/tpetra/src/Tpetra_config.h],
                                                                    [enabletpetra=yes],
                                                                    [AC_CHECK_HEADER([$withtrilinosdir/packages/tpetra/core/src/TpetraCore_config.h],
                                                                                     [enabletpetra=yes],
                                                                                     [AC_CHECK_HEADER([$withtrilinosdir/include/TpetraCore_config.h],
                                                                                                      [enabletpetra=yes],
                                                                                                      [enabletpetra=no])])])])])

                  AS_IF([test "x$enabletpetra" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_TPETRA, 1, [Flag indicating whether the library shall be compiled to use the Trilinos TPetra package])
                          AC_MSG_RESULT(<<< Configuring library with TPetra support >>>)
                        ])

                  dnl ------------------------------------------------------
                  dnl DTK
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/DataTransferKit_config.hpp],
                                  [enabledtk=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/DataTransferKit_config.hpp],
                                                   [enabledtk=yes],
                                                   [AC_CHECK_HEADER([$withtrilinosdir/DataTransferKit/src/DataTransferKit_config.hpp],
                                                                    [enabledtk=yes],
                                                                    [AC_CHECK_HEADER([$withtrilinosdir/DataTransferKit/packages/Utils/src/DataTransferKitUtils_config.hpp],
                                                                                     [enabledtk=yes],
                                                                                     [AC_CHECK_HEADER([$withtrilinosdir/include/DataTransferKitUtils_config.hpp],
                                                                                                      [enabledtk=yes],
                                                                                                      [enabledtk=no])])])])])

                  AS_IF([test "x$enabledtk" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_DTK, 1, [Flag indicating whether the library shall be compiled to use the Trilinos DTK interfaces])
                          AC_MSG_RESULT(<<< Configuring library with DTK support >>>)
                        ])


                  dnl ------------------------------------------------------
                  dnl Ifpack - We are only going to look in one place for this,
                  dnl as I don't have access to lots of older versions of
                  dnl Trilinos to know where else it might be...
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/Ifpack_config.h],
                                  [enableifpack=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/packages/ifpack/src/Ifpack_config.h],
                                                   [enableifpack=yes],
                                                   [enableifpack=no])])

                  AS_IF([test "x$enableifpack" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_IFPACK, 1, [Flag indicating whether the library shall be compiled to use the Trilinos Ifpack interfaces])
                          AC_MSG_RESULT([<<< Configuring library with Trilinos Ifpack support >>>])
                        ])

                  dnl ------------------------------------------------------
                  dnl EpetraExt - We are only going to look in one place for this,
                  dnl as I don't have access to lots of older versions of
                  dnl Trilinos to know where else it might be...
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/EpetraExt_config.h],
                                  [enableepetraext=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/packages/epetraext/src/EpetraExt_config.h],
                                                   [enableepetraext=yes],
                                                   [enableepetraext=no])])

                  AS_IF([test "x$enableepetraext" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_EPETRAEXT, 1, [Flag indicating whether the library shall be compiled to use the Trilinos EpetraExt interfaces])
                          AC_MSG_RESULT([<<< Configuring library with Trilinos EpetraExt support >>>])
                        ])

                  dnl ------------------------------------------------------
                  dnl Epetra - Trilinos can be built without Epetra, but
                  dnl libmesh can't do much with such a build.  There are several
                  dnl header files whose absence indicates Epetra is not available,
                  dnl we just choose one here that libmesh actually includes.
                  dnl ------------------------------------------------------
                  AC_CHECK_HEADER([$withtrilinosdir/include/Epetra_config.h],
                                  [enableepetra=yes],
                                  [AC_CHECK_HEADER([$withtrilinosdir/packages/epetra/src/Epetra_config.h],
                                                   [enableepetra=yes],
                                                   [enableepetra=no])])

                  AS_IF([test "x$enableepetra" != "xno"],
                        [
                          AC_DEFINE(TRILINOS_HAVE_EPETRA, 1, [Flag indicating whether the library shall be compiled to use Epetra interface in Trilinos])
                          AC_MSG_RESULT([<<< Configuring library with Trilinos Epetra support >>>])
                        ])
                ])
        ],
        [
          enabletrilinos10=no
        ])

  AC_SUBST(TRILINOS_MAKEFILE_EXPORT)

  dnl get requisite include and library variables by snarfing
  dnl them from the exported makefiles
  AS_IF([test "x$enabletrilinos10" != "xno"],
        [
          printf '%s\n' "include $TRILINOS_MAKEFILE_EXPORT" > Makefile_config_trilinos
          printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(Trilinos_LIBRARIES) \$(Trilinos_LIBRARY_DIRS) \$(Trilinos_TPL_LIBRARIES) \$(Trilinos_TPL_LIBRARY_DIRS)" >> Makefile_config_trilinos
          printf '%s\n' "echo_include:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(Trilinos_INCLUDE_DIRS) \$(Trilinos_TPL_INCLUDE_DIRS)" >> Makefile_config_trilinos

          TRILINOS_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
          TRILINOS_LIBS=`make -sf Makefile_config_trilinos echo_libs`
          rm -f Makefile_config_trilinos
       ])

  dnl Add an rpath for $withtrilinosdir/lib to the link line.  You don't
  dnl need this if Trilinos is built with static libs or if you can rely
  dnl on {DYLD,LD}_LIBRARY_PATH, but we don't want to assume either.
  AS_IF([test "x$RPATHFLAG" != "x" && test -d ${withtrilinosdir}/lib],
        [TRILINOS_LIBS="${RPATHFLAG}${withtrilinosdir}/lib $TRILINOS_LIBS"])

  AC_SUBST(TRILINOS_LIBS)
  AC_SUBST(TRILINOS_INCLUDES)
])
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl Trilinos 9
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS_9],
[
  AS_IF([test "x$TRILINOS_DIR" = "x"],
        [TRILINOS_DIR=no])

  AC_ARG_WITH(trilinos,
              AS_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  dnl AztecOO
  AS_IF([test "x$withtrilinosdir" != "xno"],
        [
          AS_IF([test -r $withtrilinosdir/include/Makefile.export.aztecoo],
                [AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.aztecoo],
                [test -r $withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo],
                [AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo],
                [enableaztecoo=no])

          AS_IF([test "x$enableaztecoo" != "xno"],
                [
                  enableaztecoo=yes
                  AC_DEFINE(TRILINOS_HAVE_AZTECOO, 1, [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
                  AC_DEFINE(HAVE_TRILINOS, 1, [])
                  AC_MSG_RESULT(<<< Configuring library with Trilinos 9 support >>>)
                ])
        ],
        [enableaztecoo=no])

  dnl Nox
  AC_ARG_WITH(nox,
              AS_HELP_STRING([--with-nox=PATH],[Specify the path to Nox installation]),
              withnoxdir=$withval,
              withnoxdir=$TRILINOS_DIR)

  AS_IF([test "x$withnoxdir" != "xno"],
        [
          AS_IF([test -r $withnoxdir/include/Makefile.export.nox],
                [NOX_MAKEFILE_EXPORT=$withnoxdir/include/Makefile.export.nox],
                [test -r $withnoxdir/packages/nox/Makefile.export.nox],
                [NOX_MAKEFILE_EXPORT=$withnoxdir/packages/nox/Makefile.export.nox],
                [enablenox=no])

          AS_IF([test "$enablenox" != "xno"],
                [
                  enablenox=yes
                  AC_DEFINE(TRILINOS_HAVE_NOX, 1, [Flag indicating whether the library shall be compiled to use the Nox solver collection])
                  AC_MSG_RESULT(<<< Configuring library with Nox support >>>)
                ])
        ],
        [enablenox=no])

  dnl ML
  AC_ARG_WITH(ml,
              AS_HELP_STRING([--with-ml=PATH],[Specify the path to ML installation]),
              withmldir=$withval,
              withmldir=$TRILINOS_DIR)

  AS_IF([test "x$withmldir" != "xno"],
        [
          AS_IF([test -r $withmldir/include/Makefile.export.ml],
                [ML_MAKEFILE_EXPORT=$withmldir/include/Makefile.export.ml],
                [test -r $withmldir/packages/nox/Makefile.export.ml],
                [ML_MAKEFILE_EXPORT=$withmldir/packages/nox/Makefile.export.ml],
                [enableml=no])

          AS_IF([test "x$enableml" != "xno"],
                [
                  enableml=yes
                  AC_DEFINE(TRILINOS_HAVE_ML, 1, [Flag indicating whether the library shall be compiled to use the ML package])
                  AC_MSG_RESULT(<<< Configuring library with ML support >>>)
                ])
        ],
        [enableml=no])

  dnl Tpetra
  AC_ARG_WITH(tpetra,
              AS_HELP_STRING([--with-tpetra=PATH],[Specify the path to Tpetra installation]),
              withtpetradir=$withval,
              withtpetradir=$TRILINOS_DIR)

  AS_IF([test "x$withtpetradir" != "xno"],
        [
          AS_IF([test -r $withtpetradir/include/Makefile.export.Tpetra],
                [TPETRA_MAKEFILE_EXPORT=$withtpetradir/include/Makefile.export.Tpetra],
                [test -r $withtpetradir/packages/tpetra/Makefile.export.Tpetra],
                [TPETRA_MAKEFILE_EXPORT=$withtpetradir/packages/tpetra/Makefile.export.Tpetra],
                [enabletpetra=no])

          AS_IF([test "x$enabletpetra" != "xno"],
                [
                  enabletpetra=yes
                  AC_DEFINE(TRILINOS_HAVE_TPETRA, 1, [Flag indicating whether the library shall be compiled to use the Tpetra solver collection])
                  AC_MSG_RESULT(<<< Configuring library with Tpetra support >>>)
                ])
       ],
       [enabletpetra=no])



  dnl DTK
  AC_ARG_WITH(dtk,
              AS_HELP_STRING([--with-dtk=PATH],[Specify the path to Dtk installation]),
              withdtkdir=$withval,
              withdtkdir=$TRILINOS_DIR)

  AS_IF([test "x$withdtkdir" != "xno"],
        [
          AS_IF([test -r $withdtkdir/include/Makefile.export.DataTransferKit],
                [DTK_MAKEFILE_EXPORT=$withdtkdir/include/Makefile.export.Makefile.export.DataTransferKit],
                [test -r $withdtkdir/DataTransferKit/Makefile.export.Makefile.export.DataTransferKit],
                [DTK_MAKEFILE_EXPORT=$withdtkdir/packages/dtk/Makefile.export.Makefile.export.DataTransferKit],
                [enabledtk=no])

          AS_IF([test "$enabledtk" != "xno"],
                [
                  enabledtk=yes
                  AC_DEFINE(TRILINOS_HAVE_DTK, 1, [Flag indicating whether the library shall be compiled to use the DataTransferKit])
                  AC_MSG_RESULT(<<< Configuring library with DTK support >>>)
                ])
        ],
        [enabledtk=no])



  dnl Ifpack - rather than try and guess how this would have worked in
  dnl Trilinos 9, we're just going to assume we don't have it.  If
  dnl anyone ever wants to go back and write a configure test with
  dnl an older Trilinos, that would be great!
  enableifpack=no

  dnl EpetraEXT - rather than try and guess how this would have worked in
  dnl Trilinos 9, we're just going to assume we don't have it.  If
  dnl anyone ever wants to go back and write a configure test with
  dnl an older Trilinos, that would be great!
  enableepetraext=no


  AC_SUBST(AZTECOO_MAKEFILE_EXPORT)

  AC_SUBST(NOX_MAKEFILE_EXPORT)

  AC_SUBST(ML_MAKEFILE_EXPORT)

  AC_SUBST(TPETRA_MAKEFILE_EXPORT)

  AC_SUBST(DTK_MAKEFILE_EXPORT)

  AS_IF([test "x$enableml" = xyes || test "x$enableaztecoo" = xyes || test "x$enablenox" = xyes],
        [enabletrilinos9=yes],
        [enabletrilinos9=no])

  dnl get requisite include and library variables by snarfing
  dnl them from the exported makefiles
  dnl
  dnl AztecOO
  AS_IF([test "x$enableaztecoo" != "xno"],
        [
          printf '%s\n' "include $AZTECOO_MAKEFILE_EXPORT" > Makefile_config_trilinos
          printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(AZTECOO_LIBS)" >> Makefile_config_trilinos
          printf '%s\n' "echo_include:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(AZTECOO_INCLUDES)" >> Makefile_config_trilinos

          AZTECOO_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
          AZTECOO_LIBS=`make -sf Makefile_config_trilinos echo_libs`
          rm -f Makefile_config_trilinos
        ])

  dnl Nox
  AS_IF([test "x$enablenox" != "xno"],
        [
          printf '%s\n' "include $NOX_MAKEFILE_EXPORT" > Makefile_config_trilinos
          printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(NOX_LIBS)" >> Makefile_config_trilinos
          printf '%s\n' "echo_include:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(NOX_INCLUDES)" >> Makefile_config_trilinos

          NOX_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
          NOX_LIBS=`make -sf Makefile_config_trilinos echo_libs`
          rm -f Makefile_config_trilinos
        ])

  dnl ML
  AS_IF([test "x$enableml" != "xno"],
        [
          printf '%s\n' "include $ML_MAKEFILE_EXPORT" > Makefile_config_trilinos
          printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(ML_LIBS)" >> Makefile_config_trilinos
          printf '%s\n' "echo_include:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(ML_INCLUDES)" >> Makefile_config_trilinos

          ML_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
          ML_LIBS=`make -sf Makefile_config_trilinos echo_libs`
          rm -f Makefile_config_trilinos
       ])

  dnl Tpetra
  AS_IF([test "x$enabletpetra" != "xno"],
        [
          printf '%s\n' "include $TPETRA_MAKEFILE_EXPORT" > Makefile_config_trilinos
          printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(TPETRA_LIBS)" >> Makefile_config_trilinos
          printf '%s\n' "echo_include:" >> Makefile_config_trilinos
          printf '\t%s\n' "@echo \$(TPETRA_INCLUDES)" >> Makefile_config_trilinos

          TPETRA_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
          TPETRA_LIBS=`make -sf Makefile_config_trilinos echo_libs`
          rm -f Makefile_config_trilinos
       ])

  AC_SUBST(AZTECOO_LIBS)
  AC_SUBST(AZTECOO_INCLUDES)
  AC_SUBST(NOX_LIBS)
  AC_SUBST(NOX_INCLUDES)
  AC_SUBST(ML_LIBS)
  AC_SUBST(ML_INCLUDES)
  AC_SUBST(TPETRA_LIBS)
  AC_SUBST(TPETRA_INCLUDES)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Trilinos -- wrapper for v9 and v10
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS],
[
  AC_ARG_ENABLE(trilinos,
                AS_HELP_STRING([--disable-trilinos],
                               [build without Trilinos support]),
                               [AS_CASE("${enableval}",
                                 [yes], [enabletrilinos=yes],
                                 [no],  [enabletrilinos=no],
                                 [AC_MSG_ERROR(bad value ${enableval} for --enable-trilinos)])],
                               [enabletrilinos=$enableoptional])

  AS_IF([test "x$enabletrilinos" = xyes],
        [
          dnl Trump --enable-trilinos with --disable-mpi
          AC_MSG_CHECKING([Whether MPI is available for Trilinos])
          AS_IF([test "x$enablempi" = "xno"],
                [
                  enabletrilinos=no
                  AC_MSG_RESULT(no)
                ],
                [
                  AC_MSG_RESULT(yes)
                ])

          dnl Trump --enable-trilinos with --enable anything except double
          dnl precision - the modules we use are hardcoded to double
          AC_MSG_CHECKING([Whether Real precision is compatible with Trilinos])
          AS_IF([test "x$enablesingleprecision" != "xno" || test "x$enabletripleprecision" != "xno" || test "x$enablequadrupleprecision" != "xno"],
                [
                  enabletrilinos=no
                  AC_MSG_RESULT(no)
                ],
                [
                  AC_MSG_RESULT(yes)
                ])
        ])

  AC_ARG_VAR([TRILINOS_DIR],  [path to Trilinos installation])

  AS_IF([test "x$enablecomplex" = "xno"],
        [
          AS_IF([test "x$enabletrilinos" != "xno"],
                [
                  dnl -- try Trilinos 10 first
                  CONFIGURE_TRILINOS_10
                  dnl -- then Trilinos 9
                  AS_IF([test "x$enabletrilinos10" = "xno"],
                        [
                          CONFIGURE_TRILINOS_9
                          AS_IF([test "x$enabletrilinos9" = "xno"], [enabletrilinos=no])
                        ])
                ])
        ])
])
