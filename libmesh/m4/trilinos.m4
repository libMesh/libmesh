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
              AC_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  if test "$withtrilinosdir" != no ; then
    AC_CHECK_FILE($withtrilinosdir/include/Makefile.export.Trilinos,
                  TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.Trilinos,
                  AC_CHECK_FILE($withtrilinosdir/Makefile.export.Trilinos,
                                TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/Makefile.export.Trilinos,
	 	                enabletrilinos10=no))

    if test "$enabletrilinos10" != no ; then
       enabletrilinos10=yes
       AC_DEFINE(HAVE_TRILINOS, 1,
                 [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
       AC_MSG_RESULT(<<< Configuring library with Trilinos 10 support >>>)
       
       dnl ------------------------------------------------------
       dnl AztecOO
       dnl ------------------------------------------------------
       AC_CHECK_FILE($withtrilinosdir/include/AztecOO_config.h,
                     enableaztecoo=yes,
                     AC_CHECK_FILE($withtrilinosdir/AztecOO_config.h,
                                   enableaztecoo=yes,
                                   AC_CHECK_FILE($withtrilinosdir/packages/aztecoo/src/AztecOO_config.h,
                                                 enableaztecoo=yes,
                                                 enableaztecoo=no)))
                     
       if test "$enableaztecoo" != no ; then
          AC_DEFINE(HAVE_AZTECOO, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos AztecOO linear solver])
          AC_MSG_RESULT(<<< Configuring library with AztecOO support >>>)
       fi
       
       dnl ------------------------------------------------------
       dnl NOX
       dnl ------------------------------------------------------
       AC_CHECK_FILE($withtrilinosdir/include/NOX_Config.h,
                     enablenox=yes,
                     AC_CHECK_FILE($withtrilinosdir/NOX_Config.h,
                                   enablenox=yes,
                                   AC_CHECK_FILE($withtrilinosdir/packages/nox/src/NOX_Config.h,
                                                 enablenox=yes,
                                                 enablenox=no)))
                     
       if test "$enablenox" != no ; then
          AC_DEFINE(HAVE_NOX, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos NOX nonlinear solver])
          AC_MSG_RESULT(<<< Configuring library with NOX support >>>)
       fi
       
       dnl ------------------------------------------------------
       dnl ML
       dnl ------------------------------------------------------
       AC_CHECK_FILE($withtrilinosdir/include/ml_include.h,
                     enableml=yes,
                     AC_CHECK_FILE($withtrilinosdir/ml_include.h,
                                   enableml=yes,
                                   AC_CHECK_FILE($withtrilinosdir/packages/ml/src/Include/ml_include.h,
                                                enableml=yes,
                                                enableml=no)))
                     
       if test "$enableml" != no ; then
          AC_DEFINE(HAVE_ML, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos ML package])
          AC_MSG_RESULT(<<< Configuring library with ML support >>>)
       fi
    fi
  else
    enabletrilinos10=no
  fi

  AC_SUBST(TRILINOS_MAKEFILE_EXPORT)
  AC_SUBST(enabletrilinos10)
  
  #########################################################
  # get requisite include and library variables by snarfing
  # them from the exported makefiles
  if (test $enabletrilinos10 != no); then
    cat <<EOF >Makefile_config_trilinos
include $TRILINOS_MAKEFILE_EXPORT
echo_libs:
	@echo \$(Trilinos_LIBRARIES) \$(Trilinos_LIBRARY_DIRS) \$(Trilinos_TPL_LIBRARIES) \$(Trilinos_TPL_LIBRARY_DIRS)

echo_include:
	@echo \$(Trilinos_INCLUDE_DIRS) \$(Trilinos_TPL_INCLUDE_DIRS)
EOF

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    TRILINOS_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    TRILINOS_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo TRILINOS_LIBS=$TRILINOS_LIBS
    #echo TRILINOS_INCLUDES=$TRILINOS_INCLUDES

    libmesh_optional_INCLUDES="$TRILINOS_INCLUDES $libmesh_optional_INCLUDES"
    libmesh_optional_LIBS="$TRILINOS_LIBS $libmesh_optional_LIBS"

    rm -f Makefile_config_trilinos
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
              AC_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  dnl AztecOO
  if test "$withtrilinosdir" != no ; then
    AC_CHECK_FILE($withtrilinosdir/include/Makefile.export.aztecoo,
                  AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.aztecoo,
                  AC_CHECK_FILE($withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo,
                                AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo,
	                        enableaztecoo=no))

    if test "$enableaztecoo" != no ; then
       enableaztecoo=yes
       AC_DEFINE(HAVE_AZTECOO, 1,
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
              AC_HELP_STRING([--with-nox=PATH],[Specify the path to Nox installation]),
              withnoxdir=$withval,
              withnoxdir=$TRILINOS_DIR)
  if test "$withnoxdir" != no ; then
    AC_CHECK_FILE($withnoxdir/include/Makefile.export.nox,
                  NOX_MAKEFILE_EXPORT=$withnoxdir/include/Makefile.export.nox,
                  AC_CHECK_FILE($withnoxdir/packages/nox/Makefile.export.nox,
                                NOX_MAKEFILE_EXPORT=$withnoxdir/packages/nox/Makefile.export.nox,
	 	                enablenox=no))

    if test "$enablenox" != no ; then
       enablenox=yes
       AC_DEFINE(HAVE_NOX, 1,
                 [Flag indicating whether the library shall be compiled to use the Nox solver collection])
       AC_MSG_RESULT(<<< Configuring library with Nox support >>>)
    fi
  else
    enablenox=no
  fi
  
  dnl ML
  AC_ARG_WITH(ml,
              AC_HELP_STRING([--with-ml=PATH],[Specify the path to ML installation]),
              withmldir=$withval,
              withmldir=$TRILINOS_DIR)
  if test "$withmldir" != no ; then
    AC_CHECK_FILE($withmldir/include/Makefile.export.ml,
                  ML_MAKEFILE_EXPORT=$withmldir/include/Makefile.export.ml,
                  AC_CHECK_FILE($withmldir/packages/nox/Makefile.export.ml,
                                ML_MAKEFILE_EXPORT=$withmldir/packages/nox/Makefile.export.ml,
	 	                enableml=no))

    if test "$enableml" != no ; then
       enableml=yes
       AC_DEFINE(HAVE_ML, 1,
                 [Flag indicating whether the library shall be compiled to use the ML package])
       AC_MSG_RESULT(<<< Configuring library with ML support >>>)
    fi
  else
    enableml=no
  fi

  AC_SUBST(AZTECOO_MAKEFILE_EXPORT)
  AC_SUBST(enableaztecoo)

  AC_SUBST(NOX_MAKEFILE_EXPORT)
  AC_SUBST(enablenox)

  AC_SUBST(ML_MAKEFILE_EXPORT)
  AC_SUBST(enableml)

  #########################################################
  # get requisite include and library variables by snarfing
  # them from the exported makefiles
  # 
  # AztecOO
  if (test $enableaztecoo != no); then
    cat <<EOF >Makefile_config_trilinos
include $AZTECOO_MAKEFILE_EXPORT
echo_libs:
	@echo \$(AZTECOO_LIBS)

echo_include:
	@echo \$(AZTECOO_INCLUDES)
EOF

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    AZTECOO_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    AZTECOO_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo AZTECOO_LIBS=$AZTECOO_LIBS
    #echo AZTECOO_INCLUDES=$AZTECOO_INCLUDES

    libmesh_optional_INCLUDES="$AZTECOO_INCLUDES $libmesh_optional_INCLUDES"
    libmesh_optional_LIBS="$AZTECOO_LIBS $libmesh_optional_LIBS"

    rm -f Makefile_config_trilinos
  fi

  # 
  # Nox
  if (test $enablenox != no); then
    cat <<EOF >Makefile_config_trilinos
include $NOX_MAKEFILE_EXPORT
echo_libs:
	@echo \$(NOX_LIBS)

echo_include:
	@echo \$(NOX_INCLUDES)
EOF

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    NOX_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    NOX_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo NOX_LIBS=$NOX_LIBS
    #echo NOX_INCLUDES=$NOX_INCLUDES

    libmesh_optional_INCLUDES="$NOX_INCLUDES $libmesh_optional_INCLUDES"
    libmesh_optional_LIBS="$NOX_LIBS $libmesh_optional_LIBS"

    rm -f Makefile_config_trilinos
  fi

  # 
  # ML
  if (test $enableml != no); then
    cat <<EOF >Makefile_config_trilinos
include $ML_MAKEFILE_EXPORT
echo_libs:
	@echo \$(ML_LIBS)

echo_include:
	@echo \$(ML_INCLUDES)
EOF

    #echo "Makefile_config_trilinos="
    #cat Makefile_config_trilinos
    ML_INCLUDES=`make -sf Makefile_config_trilinos echo_include`
    ML_LIBS=`make -sf Makefile_config_trilinos echo_libs`

    #echo ML_LIBS=$ML_LIBS
    #echo ML_INCLUDES=$ML_INCLUDES

    libmesh_optional_INCLUDES="$ML_INCLUDES $libmesh_optional_INCLUDES"
    libmesh_optional_LIBS="$ML_LIBS $libmesh_optional_LIBS"

    rm -f Makefile_config_trilinos
  fi


  AC_SUBST(AZTECOO_LIBS)
  AC_SUBST(AZTECOO_INCLUDES)
  AC_SUBST(NOX_LIBS)
  AC_SUBST(NOX_INCLUDES)
  AC_SUBST(ML_LIBS)
  AC_SUBST(ML_INCLUDES)

])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Trilinos -- wrapper for v9 and v10
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS], 
[
  AC_ARG_ENABLE(trilinos,
                AC_HELP_STRING([--enable-trilinos],
                               [build with Trilinos support]),
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
          # -- then Trlinos 9
	  if test "$enabletrilinos10" = no ; then
              CONFIGURE_TRILINOS_9
             if test "$enabletrilinos9" = no; then
	       enabletrilinos=no
	     fi
	  fi
      fi
  fi

])
dnl -------------------------------------------------------------
