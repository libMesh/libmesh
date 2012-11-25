# # -------------------------------------------------------------
# AC_DEFUN([LIBMESH_SET_METHOD],
# [
#  AC_ARG_VAR([METHOD], [method used to build libMesh (opt,dbg,devel,oprofile)])

#  # default METHOD is opt if not specified.
#  if test "x${METHOD}" = x; then
#    AC_MSG_RESULT([No build method specified, defaulting to opt])
#    METHOD=opt
#  fi 

#  # accept --with-method=METHOD.  but default to $METHOD, which is either set
#  # by the user already or defaulted to opt above
#  AC_ARG_WITH(method,
#              AC_HELP_STRING([--with-method=METHOD],
#                             [method used to build libMesh (opt,dbg,devel,oprofile)]),
#              [case "${withval}" in
#                   opt)   METHOD=opt   ;;
#                   dbg)   METHOD=dbg   ;;
#                   devel) METHOD=devel ;;
#                       *) AC_MSG_ERROR(bad value ${withval} for --with-method) ;;
#                   esac],
#                  [METHOD=$METHOD])

#  AC_MSG_RESULT([<<< Configuring libMesh in $METHOD mode >>>])

#  AM_CONDITIONAL(LIBMESH_OPT_MODE,   test x$METHOD = xopt)
#  AM_CONDITIONAL(LIBMESH_DBG_MODE,   test x$METHOD = xdbg)
#  AM_CONDITIONAL(LIBMESH_DEVEL_MODE, test x$METHOD = xdevel)

#  AC_ARG_VAR([libmesh_CPPFLAGS], [User-specified C/C++ preprocessor flags])
#  AC_ARG_VAR([libmesh_CXXFLAGS], [User-specified C++ compliation flags])
#  AC_ARG_VAR([libmesh_CFLAGS],   [User-specified C compliation flags])

#  case "${METHOD}" in
#     opt)
#        CPPFLAGS_METHOD="-DNDEBUG $libmesh_CPPFLAGS"
#        CXXFLAGS_METHOD="$CXXFLAGS_OPT $libmesh_CXXFLAGS"
#        CFLAGS_METHOD="$CFLAGS_OPT $libmesh_CFLAGS"
#        AC_CONFIG_FILES(Make.common.opt:Make.common.in)
#        ;;
#     dbg)
#        CPPFLAGS_METHOD="-DDEBUG $libmesh_CPPFLAGS"
#        CXXFLAGS_METHOD="$CXXFLAGS_DBG $libmesh_CXXFLAGS"
#        CFLAGS_METHOD="$CFLAGS_DBG $libmesh_CFLAGS"
#        AC_CONFIG_FILES(Make.common.dbg:Make.common.in)
#        ;;
#     devel)
#        CPPFLAGS_METHOD="$libmesh_CPPFLAGS"
#        CXXFLAGS_METHOD="$CXXFLAGS_DEVEL $libmesh_CXXFLAGS"
#        CFLAGS_METHOD="$CFLAGS_DEVEL $libmesh_CFLAGS"
#        AC_CONFIG_FILES(Make.common.devel:Make.common.in)
#        ;;
#     oprofile)
#        CPPFLAGS_METHOD="-DNDEBUG  $libmesh_CPPFLAGS"
#        CXXFLAGS_METHOD="$CXXFLAGS_OPT $OPROFILE_FLAGS $libmesh_CXXFLAGS"
#        CFLAGS_METHOD="$CFLAGS_OPT $OPROFILE_FLAGS $libmesh_CFLAGS"
#        AC_CONFIG_FILES(Make.common.oprof:Make.common.in)
#        ;;       

#     *)
#        AC_MSG_ERROR([bad value for METHOD: $METHOD])
#        ;;
#    esac

#   # substitute the requested method's flags
#   AC_SUBST(CPPFLAGS_METHOD)
#   AC_SUBST(CXXFLAGS_METHOD)
#   AC_SUBST(CFLAGS_METHOD)

#   # substitute all the others just in case.
#   CPPFLAGS_DBG="-DDEBUG $libmesh_CPPFLAGS"
#   AC_SUBST(CPPFLAGS_DBG)
#   AC_SUBST(CXXFLAGS_DBG)
#   AC_SUBST(CFLAGS_DBG)
  
#   CPPFLAGS_DEVEL="$libmesh_CPPFLAGS"
#   AC_SUBST(CPPFLAGS_DEVEL)
#   AC_SUBST(CXXFLAGS_DEVEL)
#   AC_SUBST(CFLAGS_DEVEL)
  
#   CPPFLAGS_OPT="-DNDEBUG $libmesh_CPPFLAGS"
#   AC_SUBST(CPPFLAGS_OPT)
#   AC_SUBST(CXXFLAGS_OPT)
#   AC_SUBST(CFLAGS_OPT)
  
# ])
# # -------------------------------------------------------------


# -------------------------------------------------------------
AC_DEFUN([LIBMESH_SET_METHODS],
[
 AC_ARG_VAR([METHODS], [methods used to build libMesh (opt,dbg,devel)])


 # accept --with-methods=METHODS.  but default to $METHODS, which is either set
 # by the user already or defaulted above
 AC_ARG_WITH(methods,
             AC_HELP_STRING([--with-methods=METHOD],
                            [methods used to build libMesh (opt,dbg,devel)]),
             [for method in ${withval} ; do
                # make sure each method specified makes sense
	        case "${method}" in
                  opt)   ;;
                  dbg)   ;;
                  devel) ;;
                      *) AC_MSG_ERROR(bad value ${method} for --with-methods) ;;
                esac
              done
	      METHODS=${withval}],
              [
	       # default METHOD is opt if not specified.
	        if test "x${METHODS}" = x; then
   		  AC_MSG_RESULT([No build methods specified, defaulting to "dbg devel opt"])
  		  METHODS="dbg devel opt"
 		fi 
	      ])

 AC_MSG_RESULT([<<< Configuring libMesh with methods "$METHODS" >>>])

 AC_ARG_VAR([libmesh_CPPFLAGS], [User-specified C/C++ preprocessor flags])
 AC_ARG_VAR([libmesh_CXXFLAGS], [User-specified C++ compliation flags])
 AC_ARG_VAR([libmesh_CFLAGS],   [User-specified C compliation flags])

 build_opt=no
 build_dbg=no
 build_devel=no
 
 # define compiler flags for all methods
 CPPFLAGS_OPT="-DNDEBUG $libmesh_CPPFLAGS"
 CXXFLAGS_OPT="$CXXFLAGS_OPT $libmesh_CXXFLAGS"
 CFLAGS_OPT="$CFLAGS_OPT $libmesh_CFLAGS"
 
 CPPFLAGS_DBG="-DDEBUG $libmesh_CPPFLAGS"
 CXXFLAGS_DBG="$CXXFLAGS_DBG $libmesh_CXXFLAGS"
 CFLAGS_DBG="$CFLAGS_DBG $libmesh_CFLAGS"
 
 CPPFLAGS_DEVEL="$libmesh_CPPFLAGS"
 CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $libmesh_CXXFLAGS"
 CFLAGS_DEVEL="$CFLAGS_DEVEL $libmesh_CFLAGS"

 AC_CONFIG_FILES(Make.common:Make.common.in)

 # conditionally compile selected methods
 for method in ${METHODS}; do
     case "${method}" in
         opt)   build_opt=yes   ;;  
         dbg)   build_dbg=yes   ;;  
         devel) build_devel=yes ;;
         *)
	     AC_MSG_ERROR(bad value ${method} for --with-methods) 
	     ;;
     esac
 done

 AM_CONDITIONAL(LIBMESH_OPT_MODE,   test x$build_opt   = xyes)
 AM_CONDITIONAL(LIBMESH_DBG_MODE,   test x$build_dbg   = xyes)
 AM_CONDITIONAL(LIBMESH_DEVEL_MODE, test x$build_devel = xyes)

  # substitute all methods
  AC_SUBST(CPPFLAGS_DBG)
  AC_SUBST(CXXFLAGS_DBG)
  AC_SUBST(CFLAGS_DBG)
  
  CPPFLAGS_DEVEL="$libmesh_CPPFLAGS"
  AC_SUBST(CPPFLAGS_DEVEL)
  AC_SUBST(CXXFLAGS_DEVEL)
  AC_SUBST(CFLAGS_DEVEL)
  
  CPPFLAGS_OPT="-DNDEBUG $libmesh_CPPFLAGS"
  AC_SUBST(CPPFLAGS_OPT)
  AC_SUBST(CXXFLAGS_OPT)
  AC_SUBST(CFLAGS_OPT)
  
])
# -------------------------------------------------------------


