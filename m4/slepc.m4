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

      if (test $slepcversion != $petscversion) ; then
        AC_MSG_RESULT(WARNING:)
        AC_MSG_RESULT(<<< Different version numbers for SLEPc and PETSc >>>)
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
      else
         enableslepc=no
      fi

      if (test "x$enableslepc" = "xyes" -a "x$includefile" != "x"); then
	  AC_MSG_RESULT(<<< Querying SLEPc configuration from $includefile >>>)
	  cat <<EOF >Makefile_config_slepc
include $includefile
getSLEPC_LIBS:
	echo \$(SLEPC_LIB) \$(ARPACK_LIB)
EOF
	SLEPC_LIBS=`make -s -f Makefile_config_slepc getSLEPC_LIBS`
        if (test x$? = x0); then
	  rm -f Makefile_config_slepc

	  #echo " "
	  #echo "SLEPC_INCLUDE=$SLEPC_INCLUDE"
	  #echo "SLEPC_LIBS=$SLEPC_LIBS"
	  #echo " "

	  AC_DEFINE(HAVE_SLEPC, 1,
                   [Flag indicating whether or not SLEPc is available])
	  AC_MSG_RESULT(<<< Configuring library with SLEPc version $slepcversion support >>>)

	  AC_SUBST(SLEPC_INCLUDE)
	  AC_SUBST(SLEPC_LIBS)
        else
          enableslepc=no
          AC_MSG_RESULT(<<< SLEPc configuration query failed. SLEPc disabled. >>>)
        fi
      fi
    fi
  fi
])
