# -------------------------------------------------------------
# Set C++ compiler flags to their default values. They will be
# modified according to other options in later steps of
# configuration
#
# CXXFLAGS_OPT    : flags for optimized mode
# CXXFLAGS_DEVEL  : flags for development mode
# CXXFLAGS_DBG    : flags for debug mode
# CPPFLAGS_OPT    : preprocessor flags for optimized mode
# CPPFLAGS_DEVEL  : preprocessor flags for development mode
# CPPFLAGS_DBG    : preprocessor flags for debug mode
# PROFILING_FLAGS : flags to enable code profiling
# ASSEMBLY_FLAGS  : flags to enable assembly language output
# WERROR_FLAGS    : flags to turn compiler warnings into errors
# PARANOID_FLAGS  : flags to turn on many more compiler warnings
#
# Usage: SET_CXX_FLAGS
#
# (Note the CXXFLAGS and the CPPFLAGS used for further tests may
#  be augmented)
# -------------------------------------------------------------
AC_DEFUN([PASSMESS_SET_CXX_FLAGS],
[
  # method-specific preprocessor flags, independent of compiler.
  CPPFLAGS_OPT="-DNDEBUG"
  CPPFLAGS_DBG="-DDEBUG"
  CPPFLAGS_DEVEL=""

  # Flag to add directories to the dynamic library search path; can
  # be changed at a later stage
  RPATHFLAG="-Wl,-rpath,"

  # Flag for profiling mode; can be modified at a later stage
  PROFILING_FLAGS="-pg"

  # Flag for assembly-output mode; can be modified at a later stage
  ASSEMBLY_FLAGS="-S"

  # Flag to turn warnings into errors; can be modified at a later stage
  WERROR_FLAGS="-Werror"

  # Flags to add every additional warning we expect the library itself
  # to not trigger
  #
  # These can be fairly compiler-specific so the default is blank: we
  # only add warnings within specific compiler version tests.
  PARANOID_FLAGS=""

  # The -g flag is necessary for OProfile to produce annotations
  # -fno-omit-frame-pointer flag turns off an optimization that
  # interferes with OProfile callgraphs
  OPROFILE_FLAGS="-g -fno-omit-frame-pointer"

  # For compilers that support it (clang >= 3.5.0 and GCC >= 4.8) the
  # user can selectively enable "sanitize" flags for different METHODs
  # with e.g.
  #
  # --enable-sanitize="dbg opt"
  #
  # These flags generally slow down code execution, so you don't
  # necessarily want to turn them on all the time or in all METHODs.

  # Declaring something AC_ARG_VAR does several things, see
  # http://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/Setting-Output-Variables.html
  # for more information... not sure we need it in this case.
  # AC_ARG_VAR([SANITIZE_METHODS], [methods we apply sanitizer flags to, e.g. "opt dbg devel"])

  AC_ARG_ENABLE([sanitize],
              AS_HELP_STRING([--enable-sanitize="opt dbg devel prof oprof"],
                             [turn on sanitizer flags for the given methods]),
              [for method in ${enableval} ; do
                 dnl make sure each method specified makes sense
                 AS_CASE("${method}",
                         [optimized|opt|debug|dbg|devel|profiling|pro|prof|oprofile|oprof], [],
                         [AC_MSG_ERROR(bad value ${method} for --enable-sanitize)])
               done
               dnl If we made it here, the case statement didn't detect any errors
               SANITIZE_METHODS=${enableval}],
               [])

  AS_IF([test "x$SANITIZE_METHODS" != x],
        [
          AC_MSG_RESULT([<<< Testing sanitizer flags for method(s) "$SANITIZE_METHODS" >>>])

          dnl Both Clang and GCC docs suggest using "-fsanitize=address -fno-omit-frame-pointer".
          dnl The Clang documentation further suggests using "-O1 -g -fno-optimize-sibling-calls".
          dnl Since these flags also work in GCC, we'll use them there as well...
          COMMON_SANITIZE_OPTIONS="-fsanitize=address -fno-omit-frame-pointer -O1 -g -fno-optimize-sibling-calls"

          dnl Test the sanitizer flags.  Currently Clang and GCC are the only
          dnl compilers that support the address sanitizer, and they use the
          dnl same set of flags.  If the set of flags used by Clang and GCC ever
          dnl diverges, we'll need to set up separate flags and test them in the
          dnl case blocks below...  The TEST_SANITIZE_FLAGS function sets
          dnl the variable have_address_sanitizer to either "no" or "yes"
          TEST_SANITIZE_FLAGS([$COMMON_SANITIZE_OPTIONS])

          dnl Enable the address sanitizer stuff if the test code compiled
          AS_IF([test "x$have_address_sanitizer" = xyes],
                [
                  dnl As of clang 3.9.0 or so, we also need to pass the sanitize flag to the linker
                  dnl if it is being used during compiling. It seems that we do not need to pass all
                  dnl the flags above, just the sanitize flag itself.
                  libmesh_LDFLAGS="$libmesh_LDFLAGS -Wc,-fsanitize=address"

                  for method in ${SANITIZE_METHODS}; do
                      AS_CASE("${method}",
                              [optimized|opt],      [SANITIZE_OPT_FLAGS=$COMMON_SANITIZE_OPTIONS],
                              [debug|dbg],          [SANITIZE_DBG_FLAGS=$COMMON_SANITIZE_OPTIONS],
                              [devel],              [SANITIZE_DEVEL_FLAGS=$COMMON_SANITIZE_OPTIONS],
                              [profiling|pro|prof], [SANITIZE_PROF_FLAGS=$COMMON_SANITIZE_OPTIONS],
                              [oprofile|oprof],     [SANITIZE_OPROF_FLAGS=$COMMON_SANITIZE_OPTIONS],
                              [AC_MSG_ERROR(bad value ${method} for --enable-sanitize)])
                  done
                ])
        ])

  # in the case blocks below we may add GLIBCXX-specific pedantic debugging preprocessor
  # definitions. Should only be used by a knowing user, because these options break
  # ABI compatibility.
  AC_ARG_ENABLE(glibcxx-debugging,
                [AS_HELP_STRING([--enable-glibcxx-debugging],
                                [enable -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC in dbg mode])],
                [AS_CASE("${enableval}",
                         [yes], [enableglibcxxdebugging=yes],
                         [no],  [enableglibcxxdebugging=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-glibcxx-debugging)])],
                [enableglibcxxdebugging=no])

  # GLIBCXX debugging causes untold woes on mac machines - so disable it
  AS_IF([test `uname` = "Darwin"],
        [
          AC_MSG_RESULT(<<< Disabling GLIBCXX debugging on Darwin >>>)
          enableglibcxxdebugging=no
        ])
  AM_CONDITIONAL(PASSMESS_ENABLE_GLIBCXX_DEBUGGING, test x$enableglibcxxdebugging = xyes)


  # First the flags for gcc compilers
  AS_IF([test "$GXX" = "yes" && test "x$REAL_GXX" != "x"],
        [
          CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization"
          dnl devel flags are added on two lines since there are so many
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -O2 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused"
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -Wpointer-arith -Wformat -Wparentheses -Wuninitialized -funroll-loops -fstrict-aliasing -Woverloaded-virtual -Wdisabled-optimization"
          CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Woverloaded-virtual"
          NODEPRECATEDFLAG="-Wno-deprecated"

          CFLAGS_OPT="-O2 -funroll-loops -fstrict-aliasing"
          CFLAGS_DEVEL="$CFLAGS_OPT -g -Wimplicit -funroll-loops -fstrict-aliasing"
          CFLAGS_DBG="-g -Wimplicit"
          ASSEMBLY_FLAGS="$ASSEMBLY_FLAGS -fverbose-asm"

          dnl Tested on gcc 4.8.5; hopefully the other 4.8.x and all
          dnl later versions support these too:
          PARANOID_FLAGS="-Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2"
          PARANOID_FLAGS="$PARANOID_FLAGS -Wformat-nonliteral -Wformat-security -Wformat-y2k"
          PARANOID_FLAGS="$PARANOID_FLAGS -Winvalid-pch -Wlogical-op -Wmissing-field-initializers"
          PARANOID_FLAGS="$PARANOID_FLAGS -Wmissing-include-dirs -Wpacked -Wredundant-decls"
          PARANOID_FLAGS="$PARANOID_FLAGS -Wshadow -Wstack-protector -Wstrict-aliasing -Wswitch-default"
          PARANOID_FLAGS="$PARANOID_FLAGS -Wtrigraphs -Wunreachable-code -Wunused-label"
          PARANOID_FLAGS="$PARANOID_FLAGS -Wunused-parameter -Wunused-value -Wvariadic-macros"
          PARANOID_FLAGS="$PARANOID_FLAGS -Wvolatile-register-var -Wwrite-strings"

          AS_IF([test "x$enableglibcxxdebugging" = "xyes"],
                [CPPFLAGS_DBG="$CPPFLAGS_DBG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"])

          # GCC 4.6.3 warns about variadic macros but supports them just
          # fine, so let's turn off that warning.
          AS_CASE("$GXX_VERSION",
                  [gcc4.6 | gcc5], [CXXFLAGS_OPT="$CXXFLAGS_OPT -Wno-variadic-macros"
                                    CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -Wno-variadic-macros"
                                    CXXFLAGS_DBG="$CXXFLAGS_DBG -Wno-variadic-macros"])

          dnl Set OS-specific flags for linkers & other stuff
          dnl For Solaris we need to pass a different flag to the linker for specifying the
          dnl dynamic library search path and add -lrpcsvc to use XDR
          AS_CASE("$target",
                  [*solaris*], [RPATHFLAG="-Wl,-R,"
                                LIBS="-lrpcsvc $LIBS"])
        ],
        [
    dnl Non-gcc compilers
    AS_CASE("$GXX_VERSION",
            [ibm_xlc], [
                          CXXFLAGS_OPT="-O3 -qmaxmem=-1 -w -qansialias -Q=10 -qrtti=all -qstaticinline"
                          CXXFLAGS_DBG="-qmaxmem=-1 -qansialias -qrtti=all -g -qstaticinline"
                          CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
                          NODEPRECATEDFLAG=""
                          CFLAGS_OPT="-O3 -qmaxmem=-1 -w -qansialias -Q=10"
                          CFLAGS_DBG="-qansialias -g"
                          CFLAGS_DEVEL="$CFLAGS_DBG"
                        ],

            dnl All Intel ICC/ECC flavors
            [intel_*], [
                          dnl Intel understands the gcc-like no-deprecated flag
                          NODEPRECATEDFLAG="-Wno-deprecated"

                          dnl Intel compilers use -qp for profiling
                          PROFILING_FLAGS="-qp"

                          dnl Intel options for annotated assembly
                          ASSEMBLY_FLAGS="$ASSEMBLY_FLAGS -fverbose-asm -fsource-asm"

                          dnl The -g flag is all OProfile needs to produce annotations
                          OPROFILE_FLAGS="-g"

                          dnl A dozen or so g++-supported warnings aren't supported on
                          dnl all icpc versions
                          PARANOID_FLAGS="-Wall -Wextra -Wdisabled-optimization -Wformat=2"
                          PARANOID_FLAGS="$PARANOID_FLAGS -Wformat-security"
                          PARANOID_FLAGS="$PARANOID_FLAGS -Winvalid-pch"
                          PARANOID_FLAGS="$PARANOID_FLAGS -Wmissing-include-dirs"
                          PARANOID_FLAGS="$PARANOID_FLAGS -Wtrigraphs"
                          PARANOID_FLAGS="$PARANOID_FLAGS -Wunused-parameter"
                          PARANOID_FLAGS="$PARANOID_FLAGS -Wwrite-strings"

                          dnl Disable some warning messages on Intel compilers:
                          dnl 161:  unrecognized pragma GCC diagnostic warning "-Wdeprecated-declarations"
                          dnl 175:  subscript out of range
                          dnl       FIN-S application code causes many false
                          dnl       positives with this
                          dnl 266:  function declared implicitly
                          dnl       Metis function "GKfree" caused this error
                          dnl       in almost every file.
                          dnl 488:  template parameter "Scalar1" is not used in declaring the
                          dnl       parameter types of function template
                          dnl       This warning was generated from one of the type_vector.h
                          dnl       constructors that uses some SFINAE tricks.
                          dnl 1476: field uses tail padding of a base class
                          dnl 1505: size of class is affected by tail padding
                          dnl       simply warns of a possible incompatibility with
                          dnl       the g++ ABI for this case
                          dnl 1572: floating-point equality and inequality comparisons are unreliable
                          dnl       Well, duh, when the tested value is computed...  OK when it
                          dnl       was from an assignment.
                          AS_CASE("$GXX_VERSION",
                                  [intel_icc_v13.x | intel_icc_v14.x | intel_icc_v15.x | intel_icc_v16.x | intel_icc_v17.x | intel_icc_v18.x],
                                  [
                                    PROFILING_FLAGS="-p"
                                    CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -g -wd175 -wd1476 -wd1505 -wd1572 -wd488 -wd161"
                                    CXXFLAGS_OPT="$CXXFLAGS_OPT -O3 -unroll -w0 -ftz"
                                    CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -w1 -g -wd175 -wd1476 -wd1505 -wd1572 -wd488 -wd161"
                                    CFLAGS_DBG="$CFLAGS_DBG -w1 -g -wd266 -wd1572 -wd488 -wd161"
                                    CFLAGS_OPT="$CFLAGS_OPT -O3 -unroll -w0 -ftz"
                                    CFLAGS_DEVEL="$CFLAGS_DBG"
                                  ],
                                  [AC_MSG_RESULT(Unknown Intel compiler, "$GXX_VERSION")])
                       ],

            [portland_group], [
                                CXXFLAGS_DBG="$CXXFLAGS_DBG -g --no_using_std"
                                CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 --no_using_std -fast -Minform=severe"
                                CXXFLAGS_DEVEL="$CXXFLAGS_DBG"

                                dnl PG C++ definitely doesnt understand -Wno-deprecated...
                                NODEPRECATEDFLAG=""

                                CFLAGS_DBG="$CFLAGS_DBG -g"
                                CFLAGS_OPT="$CFLAGS_OPT -O2"
                                CFLAGS_DEVEL="$CFLAGS_DBG"

                                dnl Disable exception handling if we dont use it
                                AS_IF([test "$enableexceptions" = no],
                                      [
                                        CXXFLAGS_DBG="$CXXFLAGS_DBG --no_exceptions"
                                        CXXFLAGS_OPT="$CXXFLAGS_OPT --no_exceptions"
                                      ])
                              ],

            [cray_cc], [
                         CXXFLAGS_DBG="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
                         CXXFLAGS_OPT="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
                         CXXFLAGS_DEVEL="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
                         NODEPRECATEDFLAG=""
                         CFLAGS_DBG="-G n"
                         CFLAGS_OPT="-G n"
                         CFLAGS_DEVEL="-G n"
                       ],

            [clang], [
                       CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors -Qunused-arguments -Wunused-parameter -Wunused"
                       dnl devel flags are added on two lines since there are so many
                       CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -O2 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long"
                       CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -Wunused-parameter -Wunused -Wpointer-arith -Wformat -Wparentheses -Wuninitialized -Qunused-arguments -Woverloaded-virtual -fno-limit-debug-info"
                       dnl dbg flags are added on two lines since there are so many
                       CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long"
                       CXXFLAGS_DBG="$CXXFLAGS_DBG -Wunused-parameter -Wunused -Wpointer-arith -Wformat -Wparentheses -Qunused-arguments -Woverloaded-virtual -fno-limit-debug-info"
                       NODEPRECATEDFLAG="-Wno-deprecated"

                       dnl Tested on clang 3.4.2
                       PARANOID_FLAGS="-Wall -Wextra -Wcast-align -Wdisabled-optimization -Wformat=2"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Wformat-nonliteral -Wformat-security -Wformat-y2k"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Winvalid-pch -Wmissing-field-initializers"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Wmissing-include-dirs -Wpacked"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Wstack-protector -Wtrigraphs"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Wunreachable-code -Wunused-label"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Wunused-parameter -Wunused-value -Wvariadic-macros"
                       PARANOID_FLAGS="$PARANOID_FLAGS -Wvolatile-register-var -Wwrite-strings"

                       CFLAGS_OPT="-O2 -Qunused-arguments -Wunused"
                       CFLAGS_DEVEL="$CFLAGS_OPT -g -Wimplicit -fno-limit-debug-info -Wunused"
                       CFLAGS_DBG="-g -Wimplicit -Qunused-arguments -fno-limit-debug-info -Wunused"
                     ],

            dnl default case
            [
              AC_MSG_RESULT(No specific options for this C++ compiler known)
              CXXFLAGS_DBG="$CXXFLAGS"
              CXXFLAGS_OPT="$CXXFLAGS"
              CXXFLAGS_DEVEL="$CXXFLAGS"
              NODEPRECATEDFLAG=""

              CFLAGS_DBG="$CFLAGS"
              CFLAGS_OPT="$CFLAGS"
              CFLAGS_DEVEL="$CFLAGS"
            ])
  ])
])
