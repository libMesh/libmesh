# ----------------------------------------------------------------
# The Reduced Basis code uses Cap'n Proto to write training data to
# disk.  By default, we check for the Capn'n Proto installation files
# in the directory provided to configure via --with-capnproto=xxx,
# which is superseded by $CAPNPROTO_DIR if it is set.
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_CAPNPROTO],
[
  AC_ARG_ENABLE(capnproto,
                AS_HELP_STRING([--disable-capnproto],
                               [build without Cap'n Proto support]),
                [case "${enableval}" in
                  yes) enablecapnproto=yes ;;
                   no) enablecapnproto=no ;;
                    *) AC_MSG_ERROR(bad value ${enableval} for --enable-capnproto) ;;
                 esac],
                [enablecapnproto=$enableoptional])

  # Cap'n Proto code uses 'auto', and therefore requires C++11 support.
  if (test "x$HAVE_CXX11" = "x" -o "x$HAVE_CXX11" = "x0"); then
    enablecapnproto=no
    AC_MSG_RESULT([<<< Cap'n Proto disabled -- C++11 support (auto keyword) is required. >>>])
  fi

  if (test $enablecapnproto = yes); then
    # The Cap'n Proto include and library paths
    CAPNPROTO_INC=""
    CAPNPROTO_LIB=""

    # The user can specify the location of their capnproto installation by configuring with
    # --with-capnproto=/path/to/capnproto/installation
    AC_ARG_WITH(capnproto,
                AS_HELP_STRING([--with-capnproto=PATH],[Specify location of the CAPNPROTO installation]),
                [
                  CAPNPROTO_INC="$withval/include"
                  CAPNPROTO_LIB="$withval/lib"
                ],
                [])

    # Let CAPNPROTO_DIR possibly override the user's --with-capnproto setting.
    if test "x$CAPNPROTO_DIR" != x; then
      CAPNPROTO_INC="$CAPNPROTO_DIR/include"
      CAPNPROTO_LIB="$CAPNPROTO_DIR/lib"
    fi

    # Initialize Makefile/config.h substitution variables.  These will
    # include flags and library names where necessary.
    CAPNPROTO_INCLUDE=""
    CAPNPROTO_LIBRARY=""

    # Check for existence of a header file in the specified location
    if (test -r $CAPNPROTO_INC/capnp/common.h) ; then
      enablecapnproto=yes
    else
      enablecapnproto=no
      AC_MSG_RESULT([<<< Required header files not found, Cap'n Proto support disabled. >>>])
    fi

    # If the header file was found, continue.
    if (test x$enablecapnproto = xyes); then
      #             Var,           look for,   name if found,   name if not, where
      AC_CHECK_PROG(CAPNP_BINARY,  capnp,      capnp,           none,        $PATH)
      if test "$CAPNP_BINARY" = capnp; then
        enablecapnproto=yes
      else
        enablecapnproto=no
        AC_MSG_RESULT([<<< The 'capnp' utility is not in your PATH, Cap'n Proto support disabled. >>>])
      fi
    fi

    # If the required programs were found in $PATH, continue.
    if (test x$enablecapnproto = xyes); then
      # Create a simple capnp schema file, make sure we can process it
      # with the capnp utility.  This will only work if $CAPNPROTO_DIR/bin is
      # in the user's PATH.
      printf '%s\n' "@0xdbb9ad1f14bf0b36;" > example.capnp
      printf '%s\n' "struct Person {" >> example.capnp
      printf '%s\n' "name @0 :Text;" >> example.capnp
      printf '%s\n' "birthdate @3 :Date;" >> example.capnp
      printf '%s\n' "email @1 :Text;" >> example.capnp
      printf '%s\n' "phones @2 :List(PhoneNumber);" >> example.capnp
      printf '%s\n' "struct PhoneNumber {" >> example.capnp
      printf '%s\n' "number @0 :Text;" >> example.capnp
      printf '%s\n' "type @1 :Type;" >> example.capnp
      printf '%s\n' "enum Type {" >> example.capnp
      printf '%s\n' "mobile @0;" >> example.capnp
      printf '%s\n' "home @1;" >> example.capnp
      printf '%s\n' "work @2;" >> example.capnp
      printf '%s\n' "} } }" >> example.capnp
      printf '%s\n' "struct Date {" >> example.capnp
      printf '%s\n' "year @0 :Int16;" >> example.capnp
      printf '%s\n' "month @1 :UInt8;" >> example.capnp
      printf '%s\n' "day @2 :UInt8;" >> example.capnp
      printf '%s\n' "}" >> example.capnp

      # Call the capnp utility
      capnp compile -oc++ example.capnp

      if (test "x$?" = "x0"); then
        enablecapnproto=yes
      else
        enablecapnproto=no
        AC_MSG_RESULT([<<< The 'capnp' utility function failed, Cap'n Proto support disabled. >>>])
      fi

      # Remove temporary files
      rm -f example.capnp example.capnp.h example.capnp.c++
    fi

    # If calling the capnp utility worked, continue.
    if (test x$enablecapnproto = xyes); then
      CAPNPROTO_INCLUDE="-I$CAPNPROTO_INC"
      CAPNPROTO_LIBRARY="-L$CAPNPROTO_LIB -lcapnp -lkj"

      # add the CAPNPROTO_LIB dir to the linker run path if it is a directory
      if (test "x$RPATHFLAG" != "x" -a -d $CAPNPROTO_LIB); then
        CAPNPROTO_LIBRARY="${RPATHFLAG}${CAPNPROTO_LIB} $CAPNPROTO_LIBRARY"
      fi

      AC_DEFINE(HAVE_CAPNPROTO, 1, [Flag indicating whether the library will be compiled with CAPNPROTO support])
      AC_MSG_RESULT([<<< Configuring library with CAPNPROTO support >>>])
    fi
  fi

  # Substitute the substitution variables
  AC_SUBST(CAPNPROTO_INCLUDE)
  AC_SUBST(CAPNPROTO_LIBRARY)
])
