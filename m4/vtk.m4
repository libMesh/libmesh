dnl Declare compilation test function preamble to
dnl be used later with AC_LANG_PROGRAM
m4_define([_AX_CXX_COMPILE_VTK_preamble],
          [[
             @%:@include "vtkSmartPointer.h"
             @%:@include "vtkCellArray.h"
             @%:@include "vtkUnstructuredGrid.h"
             @%:@include "vtkPoints.h"
             @%:@include "vtkDoubleArray.h"
             @%:@include "vtkXMLPUnstructuredGridWriter.h"
             @%:@include "vtkImageThreshold.h"
          ]])

dnl Declare compilation test function body to
dnl be used later with AC_LANG_PROGRAM
m4_define([_AX_CXX_COMPILE_VTK_body],
          [[
             vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
             vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
             vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
             vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
             vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
             vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
          ]])

dnl ----------------------------------------------------------------
dnl VTK Mesh I/O API (by Wout Ruijter) requires VTK headers and lib
dnl ----------------------------------------------------------------
AC_DEFUN([CONFIGURE_VTK],
[
  AC_ARG_VAR([VTK_INCLUDE], [path to VTK header files])
  AC_ARG_VAR([VTK_DIR],     [path to VTK installation])

  AC_ARG_ENABLE(vtk,
                AS_HELP_STRING([--disable-vtk],
                               [build without VTK file I/O support]),
               [case "${enableval}" in
                 yes)  enablevtk=yes ;;
                  no)  enablevtk=no ;;
                   *)  AC_MSG_ERROR(bad value ${enableval} for --enable-vtk) ;;
                esac],
                [enablevtk=$enableoptional])


  if (test $enablevtk = yes); then

    dnl Honor VTK_DIR if it is set
    if test "x$VTK_DIR" = x; then
       VTK_DIR=/usr
    fi

    dnl Look for VTK location in the environment, then default paths
    VTK_LS_CHECK=$(dirname $(ls -d $VTK_DIR/include/vtk*/vtkConfigure.h 2>/dev/null | tail -n 1) 2>/dev/null)
    if test "x$VTK_INC" = x; then
      if test "x$VTK_INCLUDE" != x; then
        if test -d $VTK_INCLUDE ; then
          VTK_INC=$VTK_INCLUDE
        fi
      elif (test -d $VTK_LS_CHECK); then
        VTK_INC=$VTK_LS_CHECK
      elif test "x$VTK_DIR" != x; then
        VTK_INC=$VTK_DIR/include
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
                AS_HELP_STRING([--with-vtk-include=PATH],[Specify the path for VTK header files]),
                withvtkinc=$withval,
                withvtkinc=no)

    dnl User-specific library path
    AC_ARG_WITH(vtk-lib,
                AS_HELP_STRING([--with-vtk-lib=PATH],[Specify the path for VTK libs]),
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
       vtkincFound=no;
       AC_CHECK_HEADERS($VTK_INC/vtkConfigure.h, vtkincFound=yes)

       if (test $vtkincFound = no); then
         AC_MSG_RESULT(VTK header files not found!)
         enablevtk=no;
       fi

       dnl Discover the major, minor, and build versions of VTK by looking in
       dnl vtkConfigure.h and vtkVersionMacros.h
       if (test x$enablevtk = xyes); then

         dnl If we have the vtkVersionMacros.h (VTK 6.x) find the version there
         if (test -r $VTK_INC/vtkVersionMacros.h) ; then
           vtkmajor=`grep "define VTK_MAJOR_VERSION" $VTK_INC/vtkVersionMacros.h | sed -e "s/.*#define VTK_MAJOR_VERSION[ ]*//g"`
           vtkminor=`grep "define VTK_MINOR_VERSION" $VTK_INC/vtkVersionMacros.h | sed -e "s/.*#define VTK_MINOR_VERSION[ ]*//g"`
           vtkbuild=`grep "define VTK_BUILD_VERSION" $VTK_INC/vtkVersionMacros.h | sed -e "s/.*#define VTK_BUILD_VERSION[ ]*//g"`
         dnl Otherwise (VTK 5.x) find the version numbers in vtkConfigure.h
         elif (test -r $VTK_INC/vtkConfigure.h); then
           vtkmajor=`grep "define VTK_MAJOR_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/.*#define VTK_MAJOR_VERSION[ ]*//g"`
           vtkminor=`grep "define VTK_MINOR_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/.*#define VTK_MINOR_VERSION[ ]*//g"`
           vtkbuild=`grep "define VTK_BUILD_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/.*#define VTK_BUILD_VERSION[ ]*//g"`
         fi

         vtkversion=$vtkmajor.$vtkminor.$vtkbuild
         vtkmajorminor=$vtkmajor.$vtkminor

         AC_SUBST(vtkversion)
         AC_SUBST(vtkmajor)
         AC_SUBST(vtkbuild)

         AC_DEFINE_UNQUOTED(DETECTED_VTK_VERSION_MAJOR, [$vtkmajor],
           [VTK's major version number, as detected by vtk.m4])

         AC_DEFINE_UNQUOTED(DETECTED_VTK_VERSION_MINOR, [$vtkminor],
           [VTK's minor version number, as detected by vtk.m4])

         AC_DEFINE_UNQUOTED(DETECTED_VTK_VERSION_SUBMINOR, [$vtkbuild],
           [VTK's subminor version number, as detected by vtk.m4])
       fi

       if (test x$enablevtk = xyes); then
         dnl Save original value of LIBS, then append $VTK_LIB
         old_LIBS="$LIBS"
         old_CPPFLAGS="$CPPFLAGS"

         # If this compiler supports -rpath commands, create a
         # variable for them now that can be used in $LIBS below.  We
         # ran across an issue where GCC's linker actually needed
         # -rpath flags in order to *link* a test program.  From the
         # man page for GNU ld:
         #   The -rpath option is also used when locating shared objects
         #   which are needed by shared objects explicitly included in
         #   the link; see the description of the -rpath-link option.
         if (test "x$RPATHFLAG" != "x" -a -d $VTK_LIB); then
           VTK_RPATH_FLAGS="${RPATHFLAG}${VTK_LIB}"
         fi

         dnl VTK 5.x
         if (test $vtkmajor -eq 5); then
           VTK_LIBRARY="-L$VTK_LIB -lvtkIO -lvtkCommon -lvtkFiltering -lvtkImaging"

         dnl VTK 6.1.x
         elif (test $vtkmajor -eq 6 -a $vtkminor -eq 1); then
           VTK_LIBRARY_WITH_VERSION="-L$VTK_LIB -lvtkIOCore-$vtkmajorminor -lvtkCommonCore-$vtkmajorminor -lvtkCommonDataModel-$vtkmajorminor \
                                     -lvtkFiltersCore-$vtkmajorminor -lvtkIOXML-$vtkmajorminor -lvtkImagingCore-$vtkmajorminor \
                                     -lvtkIOImage-$vtkmajorminor -lvtkImagingMath-$vtkmajorminor"

           # Some Linux distributions (Arch) install VTK without the
           # "libfoo-6.x.so" naming scheme, so we try to handle that
           # situation as well.
           VTK_LIBRARY_NO_VERSION="-L$VTK_LIB -lvtkIOCore -lvtkCommonCore -lvtkCommonDataModel \
                                   -lvtkFiltersCore -lvtkIOXML -lvtkImagingCore \
                                   -lvtkIOImage -lvtkImagingMath"

         dnl VTK 6.2.x and above
         else # elif (test $vtkmajor -eq 6 -a $vtkminor -eq 2); then
           VTK_LIBRARY_WITH_VERSION="-L$VTK_LIB -lvtkIOCore-$vtkmajorminor -lvtkCommonCore-$vtkmajorminor -lvtkCommonDataModel-$vtkmajorminor \
                                     -lvtkFiltersCore-$vtkmajorminor -lvtkIOXML-$vtkmajorminor -lvtkImagingCore-$vtkmajorminor \
                                     -lvtkIOImage-$vtkmajorminor -lvtkImagingMath-$vtkmajorminor -lvtkIOParallelXML-$vtkmajorminor"

           # Some Linux distributions (Arch) install VTK without the
           # "libfoo-6.x.so" naming scheme, so we try to handle that
           # situation as well.
           VTK_LIBRARY_NO_VERSION="-L$VTK_LIB -lvtkIOCore -lvtkCommonCore -lvtkCommonDataModel \
                                   -lvtkFiltersCore -lvtkIOXML -lvtkImagingCore \
                                   -lvtkIOImage -lvtkImagingMath -lvtkIOParallelXML"
         fi

         dnl Try to compile test prog to check for existence of VTK libraries.
         dnl AC_LINK_IFELSE uses the LIBS variable.
         if (test $vtkmajor -gt 5); then
           CPPFLAGS="$CPPFLAGS -I$VTK_INC"

           dnl 1. First try linking the NO_VERSION library names
           LIBS="$old_LIBS $VTK_RPATH_FLAGS $VTK_LIBRARY_NO_VERSION"
           AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_VTK_preamble],[_AX_CXX_COMPILE_VTK_body])],
                          [enablevtk=yes
                           VTK_LIBRARY=$VTK_LIBRARY_NO_VERSION],
                          [enablevtk=no])

           dnl 2. If that didn't work, try linking the library names with version numbers.
           if (test x$enablevtk = xno); then
             LIBS="$old_LIBS $VTK_RPATH_FLAGS $VTK_LIBRARY_WITH_VERSION"
             AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_VTK_preamble],[_AX_CXX_COMPILE_VTK_body])],
                            [enablevtk=yes
                             VTK_LIBRARY=$VTK_LIBRARY_WITH_VERSION],
                            [enablevtk=no])
           fi

           dnl 3. If that also failed, print a message that we failed.
           if (test x$enablevtk = xno); then
             AC_MSG_RESULT(<<< Linking a test program against the VTK libraries failed >>>)
           fi

           dnl Reset $LIBS, $CPPFLAGS
           LIBS="$old_LIBS"
           CPPFLAGS="$old_CPPFLAGS"

         dnl Check for VTK 5.x libraries
         else
           LIBS="$old_LIBS $VTK_RPATH_FLAGS $VTK_LIBRARY"
           CPPFLAGS="$CPPFLAGS -I$VTK_INC"

           dnl AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
           AC_CHECK_LIB([vtkIO], main, [enablevtk=yes], [enablevtk=no], [-lvtkCommon -lvtkFiltering -lvtkImaging])

           if (test $enablevtk = yes); then
             AC_CHECK_LIB([vtkCommon], main, [enablevtk=yes], [enablevtk=no])
           fi

           dnl As of VTK 5.4 it seems we also need vtkFiltering
           if (test $enablevtk = yes); then
             AC_CHECK_LIB([vtkFiltering], main, [enablevtk=yes], [enablevtk=no])
           fi

           dnl Reset $LIBS, $CPPFLAGS
           LIBS="$old_LIBS"
           CPPFLAGS="$old_CPPFLAGS"
         fi
       fi

       dnl If both the header file and the required libs were found, continue.
       if (test x$enablevtk = xyes); then
         VTK_INCLUDE="-I$VTK_INC"

         if (test "x$RPATHFLAG" != "x" -a -d $VTK_LIB); then # add the VTK_LIB to the linker run path, if it is a directory
           VTK_LIBRARY="${RPATHFLAG}${VTK_LIB} $VTK_LIBRARY"
         fi

         AC_DEFINE(HAVE_VTK, 1, [Flag indicating whether the library will be compiled with VTK support])
         AC_MSG_RESULT(<<< Configuring library with VTK version $vtkversion support >>>)
       else
         VTK_LIBRARY=''
         AC_MSG_RESULT(<<< Configuring library without VTK support >>>)
       fi
    fi
  fi

  dnl Substitute the substitution variables
  AC_SUBST(VTK_INCLUDE)
  AC_SUBST(VTK_LIBRARY)
])
