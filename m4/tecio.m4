dnl -------------------------------------------------------------
dnl tecio
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECIO],
[
  AC_ARG_ENABLE(tecio,
                AS_HELP_STRING([--disable-tecio],
                               [build without Tecplot TecIO API support (from source)]),
                [case "${enableval}" in
                  yes)  enabletecio=yes ;;
                  no)  enabletecio=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tecio) ;;
                esac],
                [enabletecio=$enableoptional])

  AC_ARG_WITH(tecio-x11-include,
              AS_HELP_STRING([--with-tecio-x11-include=PATH],
                             [Path to X11 header files. E.g. /usr/include but _not_ /usr/include/X11]),
              withteciox11inc=$withval,
              withteciox11inc=no)

  # Start off with an empty TECIO_CPPFLAGS and possibly append to it.
  TECIO_CPPFLAGS=""

  # Look for an X11 installation
  if (test $enabletecio = yes); then
    # If the user specified a path to look for X11 headers in, honor
    # that and don't try to look elsewhere.  If the compilation fails
    # then they must figure out why on their own.
    if (test $withteciox11inc != no); then
      TECIO_CPPFLAGS="-I$withteciox11inc"

    # The user did not specify where to look, so see if the file
    # exists in the usual Linux location...
    elif (test -r /usr/include/X11/Intrinsic.h); then
      TECIO_CPPFLAGS="-I/usr/include"

    # ... and if not there, try the Mac (XQuartz) location.
    elif (test -r /opt/X11/include/X11/Intrinsic.h); then
      TECIO_CPPFLAGS="-I/opt/X11/include"
    fi

    # Print a status message
    if (test "x$TECIO_CPPFLAGS" != "x"); then
      AC_MSG_RESULT(<<< Testing X11 headers with $TECIO_CPPFLAGS >>>)
    fi

    # Run AC_CHECK_HEADER to see if we can actually compile a test
    # program with the TECIO_CPPFLAGS we determined.  If this doesn't
    # work, we assume the X11 installation is not working and disable
    # the tecio interface.
    saved_CPPFLAGS=$CPPFLAGS
    CPPFLAGS="$TECIO_CPPFLAGS $CPPFLAGS"

    AC_CHECK_HEADER([X11/Intrinsic.h],
                    [],
                    [enabletecio=no])

    # Restore original value of CPPFLAGS
    CPPFLAGS=$saved_CPPFLAGS
  fi

  # The TECIO API is distributed with libmesh, so we don't have to guess
  # where it might be installed...
  if (test $enabletecio = yes); then

    # Set tecio platform-specific compiler flags.
    case "${host_os}" in
      *linux*)
        TECIO_CPPFLAGS="-DLINUX $TECIO_CPPFLAGS"
        AC_CHECK_SIZEOF([void *])
        if (test $ac_cv_sizeof_void_p = 8); then
          TECIO_CPPFLAGS="-DLINUX64 $TECIO_CPPFLAGS"
        fi
        ;;

      *darwin*)
        TECIO_CPPFLAGS="-DDARWIN -DLONGIS64 $TECIO_CPPFLAGS"
        ;;

        *)
          AC_MSG_RESULT([>>> Unrecognized TecIO platform, see contrib/tecplot/tecio/Runmake for hints on how to extend <<<])
          ;;
    esac


     TECIO_INCLUDE="-I\$(top_srcdir)/contrib/tecplot/tecio/include"
     AC_DEFINE(HAVE_TECPLOT_API, 1, [Flag indicating whether the library will be compiled with Tecplot TecIO API support])
     AC_DEFINE(HAVE_TECPLOT_API_112, 1, [Flag indicating tecplot API understands newer features])
     AC_MSG_RESULT(<<< Configuring library with Tecplot TecIO support >>>)
     have_tecio=yes
  else
     TECIO_INCLUDE=""
     enabletecio=no
     have_tecio=no
  fi

  AC_SUBST(TECIO_INCLUDE)
  AC_SUBST(TECIO_CPPFLAGS)
])
