# --------------------------------------------------------------
# Make sure that the Boost installation we found can actually compile
# Howard Hinnant's unique_ptr implementation.
# --------------------------------------------------------------
AC_DEFUN([CONFIGURE_HINNANT_UNIQUE_PTR],
[
  # Howard Hinnant's unique_ptr is enabled by default, but we allow it to be explicitly disabled.
  AC_ARG_ENABLE(hinnant-unique-ptr,
                AS_HELP_STRING([--disable-hinnant-unique-ptr],
                               [build without Howard Hinnant's unique_ptr implementation]),
                [case "${enableval}" in
                  yes)  enablehinnant=yes ;;
                   no)  enablehinnant=no ;;
                    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-boost) ;;
                esac],
                enablehinnant=yes)

  # Shell variable that will eventually be used to set the AM_CONDITIONAL
  install_hinnant_unique_ptr=no

  # If the user did not explicitly disable it, do some more testing
  if (test x$enablehinnant = xyes); then

    # We require boost to be enabled for Hinnant's unique pointer to be installed.
    if (test x$enableboost = xyes); then

      # Try to compile a test program that uses the unique_ptr.hpp header file
      AC_LANG_PUSH([C++])

      # Save any original value that CXXFLAGS had
      saveCXXFLAGS="$CXXFLAGS"

      # Add location of contrib header file to CXXFLAGS
      CXXFLAGS="$saveCXXFLAGS -I$top_srcdir/contrib/unique_ptr";

      # Depending on whether the boost installation found is external
      # or internal, also add that path (see also boost.m4)
      if (test x$install_internal_boost = xyes); then
        CXXFLAGS="$CXXFLAGS -I$top_srcdir/contrib/boost/include"
      else
        CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS"
      fi

      AC_MSG_CHECKING([if Howard Hinnant's C++03 unique_ptr implementation can be compiled])

      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <iostream>
      @%:@include "unique_ptr.hpp"
      struct Foo
      {
        Foo()      { std::cout << "Foo::Foo\n";  }
        ~Foo()     { std::cout << "Foo::~Foo\n"; }
      };
      ]], [[
      {
        // up now owns a Foo
        boost::unique_ptr<Foo> up(new Foo);
      } // Foo deleted when up goes out of scope
      ]])],[
        install_hinnant_unique_ptr=yes
        AC_MSG_RESULT(yes)
      ],[
        AC_MSG_RESULT(no)
      ])

      # Restore the original flags, whatever they were.
      CXXFLAGS="$saveCXXFLAGS"

      AC_LANG_POP([C++])

      # Set the header file define and print a message on compilation success
      if (test x$install_hinnant_unique_ptr = xyes); then
        AC_DEFINE(HAVE_HINNANT_UNIQUE_PTR, 1, [Flag indicating whether the library will be compiled with Howard Hinnant's C++03 unique_ptr implementation])
        AC_MSG_RESULT([<<< Installing Howard Hinnant's unique_ptr implementation >>>])
      fi

    else
      AC_MSG_RESULT([<<< No boost, Howard Hinnant's unique_ptr implementation not available >>>])
    fi

  else
    AC_MSG_RESULT([<<< Howard Hinnant's unique_ptr implementation explicitly disabled >>>])
  fi

  # Set Makefile.am conditional variable.  This call must be outside of all if-statements
  AM_CONDITIONAL(LIBMESH_INSTALL_HINNANT_UNIQUE_PTR, test x$install_hinnant_unique_ptr = xyes)
])
