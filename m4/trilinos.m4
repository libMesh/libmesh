dnl -------------------------------------------------------------
dnl Trilinos 10
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS_10],
[
  if (test "x$TRILINOS_DIR" = "x"); then
    # Ubuntu trilinos package?
    if (test -d /usr/include/trilinos); then
      TRILINOS_DIR=/usr/include/trilinos
    else
      TRILINOS_DIR=no
    fi
  fi

  AC_ARG_WITH(trilinos,
              AS_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  if test "$withtrilinosdir" != no ; then
    if (test -r $withtrilinosdir/include/Makefile.export.Trilinos) ; then
      TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.Trilinos
    elif (test -r $withtrilinosdir/Makefile.export.Trilinos) ; then
      TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/Makefile.export.Trilinos
    else
      enabletrilinos10=no
    fi

    if test "$enabletrilinos10" != no ; then
       dnl This actually implies Trilinos 10+, as this configure
       dnl function also works with Trilinos 11.
       enabletrilinos10=yes

       dnl Try to detect the full Trilinos version string by grepping through Trilinos_version.h.
       dnl The version string is surrounded by quotes, so we use a slight variation on our usual
       dnl grep command, and pass the result through 'tr' to remove those.
       trilinosversionstring=`grep "define TRILINOS_VERSION_STRING" $withtrilinosdir/include/Trilinos_version.h | sed -e "s/#define TRILINOS_VERSION_STRING[ ]*//g" | tr -d '"'`

       dnl If we couldn't extract anything for the full version string, try to just get the major version number.
       if test "x$trilinosversionstring" = "x" ; then
         trilinosversionstring=`grep "define TRILINOS_MAJOR_VERSION" $withtrilinosdir/include/Trilinos_version.h | sed -e "s/#define TRILINOS_MAJOR_VERSION[ ]*//g"`
       fi

       dnl If we couldn't extract the actual major version number, go with "10+"
       if test "x$trilinosversionstring" = "x" ; then
         trilinosversionstring=10+
       fi

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

       if test "$enableaztecoo" != no ; then
          AC_DEFINE(TRILINOS_HAVE_AZTECOO, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos AztecOO linear solver])
          AC_MSG_RESULT(<<< Configuring library with AztecOO support >>>)
       fi

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

       if test "$enablenox" != no ; then
          AC_DEFINE(TRILINOS_HAVE_NOX, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos NOX nonlinear solver])
          AC_MSG_RESULT(<<< Configuring library with NOX support >>>)
       fi

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

       if test "$enableml" != no ; then
          AC_DEFINE(TRILINOS_HAVE_ML, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos ML package])
          AC_MSG_RESULT(<<< Configuring library with ML support >>>)
       fi

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

       if test "$enabletpetra" != no ; then
          AC_DEFINE(TRILINOS_HAVE_TPETRA, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos TPetra package])
          AC_MSG_RESULT(<<< Configuring library with TPetra support >>>)
       fi

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

       if test "$enabledtk" != no ; then
          AC_DEFINE(TRILINOS_HAVE_DTK, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos DTK interfaces])
          AC_MSG_RESULT(<<< Configuring library with DTK support >>>)
       fi


       dnl ------------------------------------------------------
       dnl Ifpack - We are only going to look in one place for this,
       dnl as I don't have access to lots of older versions of
       dnl Trilinos to know where else it might be...
       dnl ------------------------------------------------------
       AC_CHECK_HEADER([$withtrilinosdir/include/Ifpack_config.h],
                       [enableifpack=yes],
                       [enableifpack=no])

       if test "$enableifpack" != no ; then
          AC_DEFINE(TRILINOS_HAVE_IFPACK, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos Ifpack interfaces])
          AC_MSG_RESULT([<<< Configuring library with Trilinos Ifpack support >>>])
       fi

       dnl ------------------------------------------------------
       dnl EpetraExt - We are only going to look in one place for this,
       dnl as I don't have access to lots of older versions of
       dnl Trilinos to know where else it might be...
       dnl ------------------------------------------------------
       AC_CHECK_HEADER([$withtrilinosdir/include/EpetraExt_config.h],
                       [enableepetraext=yes],
                       [enableepetraext=no])

       if test "$enableepetraext" != no ; then
          AC_DEFINE(TRILINOS_HAVE_EPETRAEXT, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos EpetraExt interfaces])
          AC_MSG_RESULT([<<< Configuring library with Trilinos EpetraExt support >>>])
       fi

    fi
  else
    enabletrilinos10=no
  fi

  AC_SUBST(TRILINOS_MAKEFILE_EXPORT)

  #########################################################
  # get requisite include and library variables by snarfing
  # them from the exported makefiles
  if (test $enabletrilinos10 != no); then
    printf '%s\n' "include $TRILINOS_MAKEFILE_EXPORT" > Makefile_config_trilinos
    printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(Trilinos_LIBRARIES) \$(Trilinos_LIBRARY_DIRS) \$(Trilinos_TPL_LIBRARIES) \$(Trilinos_TPL_LIBRARY_DIRS)" >> Makefile_config_trilinos
    printf '%s\n' "echo_include:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(Trilinos_INCLUDE_DIRS) \$(Trilinos_TPL_INCLUDE_DIRS)" >> Makefile_config_trilinos

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    TRILINOS_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    TRILINOS_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo TRILINOS_LIBS=$TRILINOS_LIBS
    #echo TRILINOS_INCLUDES=$TRILINOS_INCLUDES

    rm -f Makefile_config_trilinos
  fi

  dnl Add an rpath for $withtrilinosdir/lib to the link line.  You don't
  dnl need this if Trilinos is built with static libs or if you can rely
  dnl on {DYLD,LD}_LIBRARY_PATH, but we don't want to assume either.
  if (test "x$RPATHFLAG" != "x" -a -d ${withtrilinosdir}/lib); then
    TRILINOS_LIBS="${RPATHFLAG}${withtrilinosdir}/lib $TRILINOS_LIBS"
  fi

  AC_SUBST(TRILINOS_LIBS)
  AC_SUBST(TRILINOS_INCLUDES)
])
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl Trilinos 9
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS_9],
[
  if test "x$TRILINOS_DIR" = "x"; then
    TRILINOS_DIR=no
  fi

  AC_ARG_WITH(trilinos,
              AS_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  dnl AztecOO
  if test "$withtrilinosdir" != no ; then
    if (test -r $withtrilinosdir/include/Makefile.export.aztecoo) ; then
      AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.aztecoo
    elif (test -r $withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo) ; then
      AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo
    else
      enableaztecoo=no
    fi

    if test "$enableaztecoo" != no ; then
       enableaztecoo=yes
       AC_DEFINE(TRILINOS_HAVE_AZTECOO, 1,
                 [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
       AC_DEFINE(HAVE_TRILINOS, 1,
                 [])
       AC_MSG_RESULT(<<< Configuring library with Trilinos 9 support >>>)
    fi
  else
    enableaztecoo=no
  fi

  dnl Nox
  AC_ARG_WITH(nox,
              AS_HELP_STRING([--with-nox=PATH],[Specify the path to Nox installation]),
              withnoxdir=$withval,
              withnoxdir=$TRILINOS_DIR)

  if test "$withnoxdir" != no ; then
    if (test -r $withnoxdir/include/Makefile.export.nox) ; then
      NOX_MAKEFILE_EXPORT=$withnoxdir/include/Makefile.export.nox
    elif (test -r $withnoxdir/packages/nox/Makefile.export.nox) ; then
      NOX_MAKEFILE_EXPORT=$withnoxdir/packages/nox/Makefile.export.nox
    else
      enablenox=no
    fi

    if test "$enablenox" != no ; then
       enablenox=yes
       AC_DEFINE(TRILINOS_HAVE_NOX, 1,
                 [Flag indicating whether the library shall be compiled to use the Nox solver collection])
       AC_MSG_RESULT(<<< Configuring library with Nox support >>>)
    fi
  else
    enablenox=no
  fi



  dnl ML
  AC_ARG_WITH(ml,
              AS_HELP_STRING([--with-ml=PATH],[Specify the path to ML installation]),
              withmldir=$withval,
              withmldir=$TRILINOS_DIR)

  if test "$withmldir" != no ; then
    if (test -r $withmldir/include/Makefile.export.ml) ; then
      ML_MAKEFILE_EXPORT=$withmldir/include/Makefile.export.ml
    elif (test -r $withmldir/packages/nox/Makefile.export.ml) ; then
      ML_MAKEFILE_EXPORT=$withmldir/packages/nox/Makefile.export.ml
    else
      enableml=no
    fi

    if test "$enableml" != no ; then
       enableml=yes
       AC_DEFINE(TRILINOS_HAVE_ML, 1,
                 [Flag indicating whether the library shall be compiled to use the ML package])
       AC_MSG_RESULT(<<< Configuring library with ML support >>>)
    fi
  else
    enableml=no
  fi

  dnl Tpetra
  AC_ARG_WITH(tpetra,
              AS_HELP_STRING([--with-tpetra=PATH],[Specify the path to Tpetra installation]),
              withtpetradir=$withval,
              withtpetradir=$TRILINOS_DIR)

  if test "$withtpetradir" != no ; then
    if (test -r $withtpetradir/include/Makefile.export.Tpetra) ; then
      TPETRA_MAKEFILE_EXPORT=$withtpetradir/include/Makefile.export.Tpetra
    elif (test -r $withtpetradir/packages/tpetra/Makefile.export.Tpetra) ; then
      TPETRA_MAKEFILE_EXPORT=$withtpetradir/packages/tpetra/Makefile.export.Tpetra
    else
      enabletpetra=no
    fi

    if test "$enabletpetra" != no ; then
       enabletpetra=yes
       AC_DEFINE(TRILINOS_HAVE_TPETRA, 1,
                 [Flag indicating whether the library shall be compiled to use the Tpetra solver collection])
       AC_MSG_RESULT(<<< Configuring library with Tpetra support >>>)
    fi
  else
    enabletpetra=no
  fi



  dnl DTK
  AC_ARG_WITH(dtk,
              AS_HELP_STRING([--with-dtk=PATH],[Specify the path to Dtk installation]),
              withdtkdir=$withval,
              withdtkdir=$TRILINOS_DIR)

  if test "$withdtkdir" != no ; then
    if (test -r $withdtkdir/include/Makefile.export.DataTransferKit) ; then
      DTK_MAKEFILE_EXPORT=$withdtkdir/include/Makefile.export.Makefile.export.DataTransferKit
    elif (test -r $withdtkdir/DataTransferKit/Makefile.export.Makefile.export.DataTransferKit) ; then
      DTK_MAKEFILE_EXPORT=$withdtkdir/packages/dtk/Makefile.export.Makefile.export.DataTransferKit
    else
      enabledtk=no
    fi

    if test "$enabledtk" != no ; then
       enabledtk=yes
       AC_DEFINE(TRILINOS_HAVE_DTK, 1,
                 [Flag indicating whether the library shall be compiled to use the DataTransferKit])
       AC_MSG_RESULT(<<< Configuring library with DTK support >>>)
    fi
  else
    enabledtk=no
  fi



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

  if test "x$enableml" = xyes -o "x$enableaztecoo" = xyes -o "x$enablenox" = xyes; then
    enabletrilinos9=yes
  else
    enabletrilinos9=no
  fi

  #########################################################
  # get requisite include and library variables by snarfing
  # them from the exported makefiles
  #
  # AztecOO
  if (test $enableaztecoo != no); then
    printf '%s\n' "include $AZTECOO_MAKEFILE_EXPORT" > Makefile_config_trilinos
    printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(AZTECOO_LIBS)" >> Makefile_config_trilinos
    printf '%s\n' "echo_include:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(AZTECOO_INCLUDES)" >> Makefile_config_trilinos

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    AZTECOO_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    AZTECOO_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo AZTECOO_LIBS=$AZTECOO_LIBS
    #echo AZTECOO_INCLUDES=$AZTECOO_INCLUDES

    rm -f Makefile_config_trilinos
  fi

  #
  # Nox
  if (test $enablenox != no); then
    printf '%s\n' "include $NOX_MAKEFILE_EXPORT" > Makefile_config_trilinos
    printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(NOX_LIBS)" >> Makefile_config_trilinos
    printf '%s\n' "echo_include:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(NOX_INCLUDES)" >> Makefile_config_trilinos

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    NOX_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    NOX_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo NOX_LIBS=$NOX_LIBS
    #echo NOX_INCLUDES=$NOX_INCLUDES

    rm -f Makefile_config_trilinos
  fi

  #
  # ML
  if (test $enableml != no); then
    printf '%s\n' "include $ML_MAKEFILE_EXPORT" > Makefile_config_trilinos
    printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(ML_LIBS)" >> Makefile_config_trilinos
    printf '%s\n' "echo_include:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(ML_INCLUDES)" >> Makefile_config_trilinos

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    ML_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    ML_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo ML_LIBS=$ML_LIBS
    #echo ML_INCLUDES=$ML_INCLUDES

    rm -f Makefile_config_trilinos
  fi



  #
  # Tpetra
  if (test $enabletpetra != no); then
    printf '%s\n' "include $TPETRA_MAKEFILE_EXPORT" > Makefile_config_trilinos
    printf '%s\n' "echo_libs:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(TPETRA_LIBS)" >> Makefile_config_trilinos
    printf '%s\n' "echo_include:" >> Makefile_config_trilinos
    printf '\t%s\n' "@echo \$(TPETRA_INCLUDES)" >> Makefile_config_trilinos

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    TPETRA_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    TPETRA_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo TPETRA_LIBS=$TPETRA_LIBS
    #echo TPETRA_INCLUDES=$TPETRA_INCLUDES

    rm -f Makefile_config_trilinos
  fi

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
                               [case "${enableval}" in
                                 yes)  enabletrilinos=yes ;;
                                 no)  enabletrilinos=no ;;
                                 *)  AC_MSG_ERROR(bad value ${enableval} for --enable-trilinos) ;;
                               esac],
                               [enabletrilinos=$enableoptional])

  # Trump --enable-trilinos with --disable-mpi
  if (test "x$enablempi" = xno); then
    enabletrilinos=no
  fi

  AC_ARG_VAR([TRILINOS_DIR],  [path to Trilinos installation])

  if test "$enablecomplex" = no ; then
    if test "$enabletrilinos" != no ; then
      # -- try Trilinos 10 first
      CONFIGURE_TRILINOS_10
      # -- then Trilinos 9
      if test "$enabletrilinos10" = no ; then
        CONFIGURE_TRILINOS_9
        if test "$enabletrilinos9" = no; then
          enabletrilinos=no
        fi
      fi
    fi
  fi
])
