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
                [AS_CASE("${enableval}",
                         [yes], [enablecapnproto=yes],
                         [no], [enablecapnproto=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-capnproto)])],
                [enablecapnproto=$enableoptional])

  dnl Cap'n Proto code uses 'auto', and therefore requires C++11 support.
  AS_IF([test "x$HAVE_CXX11" = "x" || test "x$HAVE_CXX11" = "x0"],
        [
          enablecapnproto=no
          AC_MSG_RESULT([<<< Cap'n Proto disabled -- C++11 support (auto keyword) is required. >>>])
        ])

  AS_IF([test "x$enablecapnproto" = "xyes"],
        [
          dnl The Cap'n Proto include and library paths
          CAPNPROTO_INC=""
          CAPNPROTO_LIB=""

          dnl The user can specify the location of their capnproto installation by configuring with
          dnl --with-capnproto=/path/to/capnproto/installation
          AC_ARG_WITH(capnproto,
                      AS_HELP_STRING([--with-capnproto=PATH],[Specify location of the CAPNPROTO installation]),
                      [
                        CAPNPROTO_INC="$withval/include"
                        CAPNPROTO_LIB="$withval/lib"
                      ],
                      [])

          dnl Let CAPNPROTO_DIR possibly override the user's --with-capnproto setting.
          AS_IF([test "x$CAPNPROTO_DIR" != x],
                [
                  CAPNPROTO_INC="$CAPNPROTO_DIR/include"
                  CAPNPROTO_LIB="$CAPNPROTO_DIR/lib"
                ])

          dnl Initialize Makefile/config.h substitution variables.  These will
          dnl include flags and library names where necessary.
          CAPNPROTO_INCLUDE=""
          CAPNPROTO_LIBRARY=""

          dnl Check for existence of a header file in the specified location
          AS_IF([test -r $CAPNPROTO_INC/capnp/common.h],
                [enablecapnproto=yes],
                [
                  enablecapnproto=no
                  AC_MSG_RESULT([<<< Required header files not found, Cap'n Proto support disabled. >>>])
                ])

          dnl If the header file was found, continue.
          AS_IF([test "x$enablecapnproto" = "xyes"],
                [
                  AC_CHECK_PROG(CAPNP_BINARY, capnp, capnp, none, $PATH)
                  AS_IF([test "x$CAPNP_BINARY" = "xcapnp"],
                        [enablecapnproto=yes],
                        [
                          enablecapnproto=no
                          AC_MSG_RESULT([<<< The 'capnp' utility is not in your PATH, Cap'n Proto support disabled. >>>])
                        ])
                ])

          # If the required programs were found in $PATH, continue.
          AS_IF([test "x$enablecapnproto" = "xyes"],
                [
                  dnl Create a simple capnp schema file, make sure we can process it
                  dnl with the capnp utility.  This will only work if $CAPNPROTO_DIR/bin is
                  dnl in the user's PATH.
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

                  AS_IF([test "x$?" = "x0"],
                        [enablecapnproto=yes],
                        [
                          enablecapnproto=no
                          AC_MSG_RESULT([<<< The 'capnp' utility function failed, Cap'n Proto support disabled. >>>])
                        ])

                  dnl Remove temporary files
                  rm -f example.capnp example.capnp.h example.capnp.c++
                ])

          dnl If calling the capnp utility worked, continue.
          AS_IF([test "x$enablecapnproto" = "xyes"],
                [
                  CAPNPROTO_INCLUDE="-I$CAPNPROTO_INC"
                  CAPNPROTO_LIBRARY="-L$CAPNPROTO_LIB -lcapnp -lkj"

                  dnl add the CAPNPROTO_LIB dir to the linker run path if it is a directory
                  AS_IF([test "x$RPATHFLAG" != "x" && test -d $CAPNPROTO_LIB],
                        [CAPNPROTO_LIBRARY="${RPATHFLAG}${CAPNPROTO_LIB} $CAPNPROTO_LIBRARY"])

                  AC_DEFINE(HAVE_CAPNPROTO, 1, [Flag indicating whether the library will be compiled with CAPNPROTO support])
                  AC_MSG_RESULT([<<< Configuring library with CAPNPROTO support >>>])
                ])
        ])

  dnl Substitute the substitution variables
  AC_SUBST(CAPNPROTO_INCLUDE)
  AC_SUBST(CAPNPROTO_LIBRARY)
])
