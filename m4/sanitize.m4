# Test that a simple test program can be compiled with the sanitizer
# options.
#
# So far, the sanitizer options have worked with all the GCC 4.8
# installations I've tried, but not all versions of clang are created
# equal.  For example, Apple's Mavericks/Yosemite clang compiler
# reports:
#
# $ /usr/bin/clang++ --version
# Apple LLVM version 6.0 (clang-600.0.51) (based on LLVM 3.5svn)
#
# i.e. LLVM 3.5, but it does not support address sanitizer.
# On the other hand, if you build clang from source with the
# address sanitizer enabled, you get:
#
# $ clang++ --version
# clang version 3.5.0 (tags/RELEASE_350/final 219211)
#
# so there is no way to tell if it supports address
# sanitizer just from the version info.
AC_DEFUN([LIBMESH_TEST_SANITIZE_FLAGS],
  [
    # Save current lang and CXXFLAGS values
    AC_LANG_PUSH([C++])
    saveCXXFLAGS="$CXXFLAGS"

    # Set CXXFLAGS to be used for the test compilation.  Note:
    # we don't test these flags in conjunction with any
    # others...  that doesn't seem to be necessary yet.
    CXXFLAGS="$1"

    # Tell the user what we are doing
    AC_MSG_CHECKING([if compiler has support for address sanitizer])

    # Try compiling and running a simple main program with
    # sanitizer flags.  Since address sanitizer requires
    # library-level support, we want to be sure that a compiled
    # executable can run, not just that the compiler accepts the
    # sanitizer flags.  This program does not actually have a
    # memory error, otherwise configure would consider the test
    # to have failed...
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
    [],
    [[
        int *array = new int[100];
        delete [] array;
    ]])],
    [
        # Result if program succeeds
        have_address_sanitizer=yes
        AC_MSG_RESULT(yes)
    ],
    [
        # Result if program fails
        have_address_sanitizer=no
        AC_MSG_RESULT(no)
    ],
    [
        # Result if cross-compiling
        have_address_sanitizer=no
        AC_MSG_RESULT(no)
    ])

    # Restore the original lang and flags
    CXXFLAGS="$saveCXXFLAGS"
    AC_LANG_POP([C++])
  ])
