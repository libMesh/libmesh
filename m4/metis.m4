dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS],
[
  AC_ARG_ENABLE(metis,
                AS_HELP_STRING([--disable-metis],
                               [build without Metis graph partitioning suppport]),
                [case "${enableval}" in
                  yes)  enablemetis=yes ;;
                  no)  enablemetis=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
                esac],
                [enablemetis=$enableoptional])

  AC_ARG_WITH(metis,
             AS_HELP_STRING([--with-metis=<internal,PETSc>],
                            [metis to use. interal: build from contrib, PETSc: rely on PETSc]),
             [case "${withval}" in
               internal)
                 build_metis=yes ;;
               PETSc)
                 build_metis=no ;;
               *)
                 AC_MSG_ERROR(bad value ${withval} for --with-metis) ;;
              esac],
              [build_metis=yes])

  # If PETSc has its own METIS, default to using that one regardless
  # of what the user specified (if anything) in --with-metis.
  if (test "x$petsc_have_metis" != "x") ; then
    if (test $petsc_have_metis -gt 0) ; then
      build_metis=no
    fi
  fi

  # Conversely, if:
  # .) METIS is enabled in libmesh,
  # .) PETSc does not have a METIS or we aren't using PETSc, and
  # .) build_metis=no because user said --with-metis=PETSc,
  # then we need to make sure that libmesh builds its own METIS!
  if (test $enablemetis = yes) ; then
    if (test $build_metis = no) ; then
      if (test "x$petsc_have_metis" = "x0") ; then
        build_metis=yes
      fi
      if (test $enablepetsc = no) ; then
        build_metis=yes
      fi
    fi
  fi

  dnl The METIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablemetis = yes); then
    dnl Only put Metis contrib directories into the compiler include paths if we're building our own Metis.
    if (test $build_metis = yes); then
      METIS_INCLUDE="-I\$(top_srcdir)/contrib/metis/include"
      METIS_LIB="\$(EXTERNAL_LIBDIR)/libmetis\$(libext) \$(EXTERNAL_LIBDIR)/libGK\$(libext)"
    fi
    AC_DEFINE(HAVE_METIS, 1, [Flag indicating whether the library will be compiled with Metis support])
    AC_MSG_RESULT(<<< Configuring library with Metis support >>>)

    dnl look for thread-local storage
    AX_TLS
 else
     METIS_INCLUDE=""
     METIS_LIB=""
     enablemetis=no
  fi
  AM_CONDITIONAL(BUILD_METIS, test x$build_metis = xyes)

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(METIS_LIB)
])
