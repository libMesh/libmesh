AC_DEFUN([LIBMESH_SET_METHODS],
[
 AC_ARG_VAR([METHODS], [methods used to build libMesh, e.g. "opt dbg devel". Possibilities include: (opt,dbg,devel,prof,oprof)])

 dnl accept --with-methods=METHODS.  but default to $METHODS, which is either set
 dnl by the user already or defaulted above
 AC_ARG_WITH(methods,
             AS_HELP_STRING([--with-methods=METHODS],
                            [methods used to build libMesh (opt,dbg,devel,prof,oprof)]),
             [for method in ${withval} ; do
                dnl make sure each method specified makes sense
                AS_CASE("${method}",
                  [optimized|opt], [],
                  [debug|dbg], [],
                  [devel], [],
                  [profiling|pro|prof], [],
                  [oprofile|oprof], [],
                  [AC_MSG_ERROR(bad value ${method} for --with-methods)])
              done
              METHODS=${withval}],
             [
               dnl default METHOD is opt if not specified.
               AS_IF([test "x${METHODS}" = x],
                     [
                       METHODS="dbg devel opt"
                       AC_MSG_RESULT([No build methods specified, defaulting to "$METHODS"])
                     ])
             ])

 AC_MSG_RESULT([<<< Configuring libMesh with methods "$METHODS" >>>])

 AC_ARG_VAR([libmesh_CPPFLAGS], [User-specified C/C++ preprocessor flags])
 AC_ARG_VAR([libmesh_CXXFLAGS], [User-specified C++ compilation flags])
 AC_ARG_VAR([libmesh_CFLAGS],   [User-specified C compilation flags])

 build_opt=no
 build_dbg=no
 build_devel=no
 build_prof=no
 build_oprof=no

 dnl define compiler flags for all methods
 CPPFLAGS_OPT="$CPPFLAGS_OPT $libmesh_CPPFLAGS"
 CXXFLAGS_OPT="$CXXFLAGS_OPT $libmesh_CXXFLAGS"
 CFLAGS_OPT="$CFLAGS_OPT $libmesh_CFLAGS"

 CPPFLAGS_DBG="$CPPFLAGS_DBG $libmesh_CPPFLAGS"
 CXXFLAGS_DBG="$CXXFLAGS_DBG $libmesh_CXXFLAGS"
 CFLAGS_DBG="$CFLAGS_DBG $libmesh_CFLAGS"

 CPPFLAGS_DEVEL="$CPPFLAGS_DEVEL $libmesh_CPPFLAGS"
 CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $libmesh_CXXFLAGS"
 CFLAGS_DEVEL="$CFLAGS_DEVEL $libmesh_CFLAGS"

 dnl profiling-specific flags are derivatives of optimized flags
 CPPFLAGS_PROF="$CPPFLAGS_OPT"
 CXXFLAGS_PROF="$CXXFLAGS_OPT $PROFILING_FLAGS"
 CFLAGS_PROF="$CFLAGS_OPT $PROFILING_FLAGS"

 CPPFLAGS_OPROF="$CPPFLAGS_OPT"
 CXXFLAGS_OPROF="$CXXFLAGS_OPT $OPROFILE_FLAGS"
 CFLAGS_OPROF="$CFLAGS_OPT $OPROFILE_FLAGS"

 dnl The PROF and OPROF methods are derived from opt, but we don't want
 dnl to turn on address sanitizer stuff in prof/oprof mode just because
 dnl it got turned on in opt mode.  Hence this second round of setting
 dnl all the flags vars...
 CXXFLAGS_OPT="$CXXFLAGS_OPT $SANITIZE_OPT_FLAGS"
 CFLAGS_OPT="$CFLAGS_OPT $SANITIZE_OPT_FLAGS"

 CXXFLAGS_DBG="$CXXFLAGS_DBG $SANITIZE_DBG_FLAGS"
 CFLAGS_DBG="$CFLAGS_DBG $SANITIZE_DBG_FLAGS"

 CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $SANITIZE_DEVEL_FLAGS"
 CFLAGS_DEVEL="$CFLAGS_DEVEL $SANITIZE_DEVEL_FLAGS"

 CXXFLAGS_PROF="$CXXFLAGS_PROF $SANITIZE_PROF_FLAGS"
 CFLAGS_PROF="$CFLAGS_PROF $SANITIZE_PROF_FLAGS"

 CXXFLAGS_OPROF="$CXXFLAGS_OPROF $SANITIZE_OPROF_FLAGS"
 CFLAGS_OPROF="$CFLAGS_OPROF $SANITIZE_OPROF_FLAGS"

 dnl conditionally compile selected methods
 for method in ${METHODS}; do
     AS_CASE("${method}",
             [optimized|opt],      [build_opt=yes],
             [debug|dbg],          [build_dbg=yes],
             [devel],              [build_devel=yes],
             [profiling|pro|prof], [build_prof=yes],
             [oprofile|oprof],     [build_oprof=yes],
             [AC_MSG_ERROR(bad value ${method} for --with-methods)])
 done

 AM_CONDITIONAL(LIBMESH_OPT_MODE,   test x$build_opt   = xyes)
 AM_CONDITIONAL(LIBMESH_DBG_MODE,   test x$build_dbg   = xyes)
 AM_CONDITIONAL(LIBMESH_DEVEL_MODE, test x$build_devel = xyes)
 AM_CONDITIONAL(LIBMESH_PROF_MODE,  test x$build_prof = xyes)
 AM_CONDITIONAL(LIBMESH_OPROF_MODE, test x$build_oprof = xyes)

 LIBMESH_PC_IN=""
  # set the configuration input file for libmesh.pc
  AS_IF([test x$build_opt = xyes], [LIBMESH_PC_IN="contrib/utils/libmesh-opt.pc.in"],
        [test x$build_dbg = xyes], [LIBMESH_PC_IN="contrib/utils/libmesh-dbg.pc.in"],
        [test x$build_devel = xyes], [LIBMESH_PC_IN="contrib/utils/libmesh-devel.pc.in"],
        [test x$build_prof = xyes], [LIBMESH_PC_IN="contrib/utils/libmesh-prof.pc.in"],
        [test x$build_oprof = xyes], [LIBMESH_PC_IN="contrib/utils/libmesh-oprof.pc.in"])

  dnl substitute all methods
  AC_SUBST(CPPFLAGS_DBG)
  AC_SUBST(CXXFLAGS_DBG)
  AC_SUBST(CFLAGS_DBG)

  AC_SUBST(CPPFLAGS_DEVEL)
  AC_SUBST(CXXFLAGS_DEVEL)
  AC_SUBST(CFLAGS_DEVEL)

  AC_SUBST(CPPFLAGS_OPT)
  AC_SUBST(CXXFLAGS_OPT)
  AC_SUBST(CFLAGS_OPT)

  AC_SUBST(CPPFLAGS_PROF)
  AC_SUBST(CXXFLAGS_PROF)
  AC_SUBST(CFLAGS_PROF)

  AC_SUBST(CPPFLAGS_OPROF)
  AC_SUBST(CXXFLAGS_OPROF)
  AC_SUBST(CFLAGS_OPROF)

])
