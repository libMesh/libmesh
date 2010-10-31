dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl ----------------------------------------------------------------
dnl VTK Mesh I/O API (by Wout Ruijter) requires VTK headers and lib
dnl ----------------------------------------------------------------
AC_DEFUN([CONFIGURE_VTK], 
[
  dnl Look for VTK location in the environment, then default paths
  VTK_LS_CHECK=$(dirname $(ls -d /usr/include/vtk*/vtkConfigure.h 2>/dev/null | tail -n 1) 2>/dev/null)
  if test "x$VTK_INC" = x; then
    if test "x$VTK_INCLUDE" != x; then
      VTK_INC=$VTK_INCLUDE
    elif test "x$VTK_DIR" != x; then
      VTK_INC=$VTK_DIR/include
    elif (test -d $VTK_LS_CHECK); then
      VTK_INC=$VTK_LS_CHECK
    else
      VTK_INC="/usr/include/vtk"
    fi
  fi

  if test "x$VTK_LIB" = x; then
    if test "x$VTK_DIR" != x; then
      VTK_LIB=$VTK_DIR/lib
    else
      VTK_LIB="/usr/lib"
    fi
  fi

  dnl User-specific include path
  AC_ARG_WITH(vtk-include,
              AC_HELP_STRING([--with-vtk-include=PATH],[Specify the path for VTK header files]),
              withvtkinc=$withval,
              withvtkinc=no)
	      
  dnl User-specific library path
  AC_ARG_WITH(vtk-lib,
              AC_HELP_STRING([--with-vtk-lib=PATH],[Specify the path for VTK libs]),
              withvtklib=$withval,
              withvtklib=no)

  if (test $withvtkinc != no); then
    VTK_INC="$withvtkinc"
  fi

  if (test $withvtklib != no); then
    VTK_LIB="$withvtklib"
  fi

  dnl Initialize Makefile/config.h substitution variables
  VTK_INCLUDE=""
  VTK_LIBRARY=""

  dnl Properly let the substitution variables
  if (test $enablevtk = yes); then
  
     dnl Check for existence of a header file in the specified location
     dnl AC_CHECK_FILE([$VTK_INC/vtkConfigure.h], [vtkincFound="OK"], [vtkincFound="FAIL"])
     vtkincFound=no;
     AC_CHECK_HEADERS($VTK_INC/vtkConfigure.h, vtkincFound=yes)

     if (test $vtkincFound = no); then
       AC_MSG_RESULT(VTK header files not found!)
       enablevtk=no;
     fi

     dnl Discover the major, minor, and build versions of VTK by looking in
     dnl vtkConfigure.h.  This may eventually be useful for linking against
     dnl different subsets of libraries.
     if (test $enablevtk = yes); then
       vtkmajor=`grep "define VTK_MAJOR_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/#define VTK_MAJOR_VERSION[ ]*//g"`
       vtkminor=`grep "define VTK_MINOR_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/#define VTK_MINOR_VERSION[ ]*//g"`
       vtkbuild=`grep "define VTK_BUILD_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/#define VTK_BUILD_VERSION[ ]*//g"`
       vtkversion=$vtkmajor.$vtkminor.$vtkbuild
       AC_MSG_RESULT(<<< Configuring library with VTK version $vtkversion support >>>)
     fi

     if (test $enablevtk = yes); then
       dnl Also Check for existence of required libraries.  This is not really
       dnl right way to do it -- it's not portable to Macs, where .so's are called
       dnl .dylib's instead.
       dnl AC_CHECK_FILE($VTK_LIB/libvtkIO.so, [enablevtk=yes], [enablevtk=no])
       dnl AC_CHECK_FILE($VTK_LIB/libvtkCommon.so, [enablevtk=yes], [enablevtk=no])
       
       dnl AC_HAVE_LIBRARY (library, [action-if-found], [action-if-not-found], [other-libraries])
       dnl Note: Basically tries to compile a function which calls main().  

       dnl Save original value of LIBS, then append $VTK_LIB
       old_LIBS="$LIBS"
       LIBS="$old_LIBS -L$VTK_LIB"

       dnl Try to compile test prog to check for existence of VTK libraries
       dnl AC_HAVE_LIBRARY uses the LIBS variable.
       AC_HAVE_LIBRARY([vtkIO], [enablevtk=yes], [enablevtk=no])
       
       if (test $enablevtk = yes); then
         AC_HAVE_LIBRARY([vtkCommon], [enablevtk=yes], [enablevtk=no])
       fi

       dnl As of VTK 5.4 it seems we also need vtkFiltering
       if (test $enablevtk = yes); then
         AC_HAVE_LIBRARY([vtkFiltering], [enablevtk=yes], [enablevtk=no])
       fi
       
       dnl Reset $LIBS
       LIBS="$old_LIBS"
     fi
     
     dnl If both the header file and the required libs were found, continue.
     if (test $enablevtk = yes); then
       VTK_INCLUDE="-I$VTK_INC"
       VTK_LIBRARY="\$(libmesh_RPATHFLAG)$VTK_LIB -L$VTK_LIB -lvtkIO -lvtkCommon -lvtkFiltering"
       AC_DEFINE(HAVE_VTK, 1, [Flag indicating whether the library will be compiled with VTK support])
       AC_MSG_RESULT(<<< Configuring library with VTK support >>>)
     fi
  fi

  dnl Substitute the substitution variables
  AC_SUBST(VTK_INCLUDE)
  AC_SUBST(VTK_LIBRARY)	
  AC_SUBST(enablevtk)
])
