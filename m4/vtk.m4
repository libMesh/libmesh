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
             @%:@include "vtkMPIController.h"
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
             vtkSmartPointer<vtkMPIController> controller = vtkSmartPointer<vtkMPIController>::New();
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
                [AS_CASE("${enableval}",
                         [yes], [enablevtk=yes],
                         [no],  [enablevtk=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-vtk)])],
                [enablevtk=$enableoptional])

  dnl Setting --enable-vtk-required causes an error to be emitted during
  dnl configure if VTK is not successfully detected. This is useful for app
  dnl codes which require VTK (like MOOSE-based apps), since it prevents
  dnl situations where libmesh is accidentally built without VTK support
  dnl (which may take a very long time), and then the app fails to compile,
  dnl requiring you to redo everything.
  AC_ARG_ENABLE(vtk-required,
                AC_HELP_STRING([--enable-vtk-required],
                               [Error if VTK is not detected by configure]),
                [AS_CASE("${enableval}",
                         [yes], [vtkrequired=yes],
                         [no],  [vtkrequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-vtk-required)])],
                     [vtkrequired=no])

  AS_IF([test "$enablevtk" = "yes"],
        [
    dnl Honor VTK_DIR if it is set
    AS_IF([test "x$VTK_DIR" = "x"], [VTK_DIR=/usr])

    dnl Look for VTK location in the environment, then default paths
    VTK_LS_CHECK=$(dirname $(ls -d $VTK_DIR/include/vtk*/vtkConfigure.h 2>/dev/null | tail -n 1) 2>/dev/null)
    AS_IF([test "x$VTK_INC" = "x"],
          [
            AS_IF([test "x$VTK_INCLUDE" != "x"], [AS_IF([test -d $VTK_INCLUDE], [VTK_INC=$VTK_INCLUDE])],
                  [test -d $VTK_LS_CHECK],       [VTK_INC=$VTK_LS_CHECK],
                  [test "x$VTK_DIR" != "x"],     [VTK_INC=$VTK_DIR/include],
                  [VTK_INC="/usr/include/vtk"])
          ])

    AS_IF([test "x$VTK_LIB" = "x"],
          [
            AS_IF([test "x$VTK_DIR" != x], [VTK_LIB=$VTK_DIR/lib],
                  [VTK_LIB="/usr/lib"])
          ])

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

    AS_IF([test "$withvtkinc" != "no"], [VTK_INC="$withvtkinc"])
    AS_IF([test "$withvtklib" != "no"], [VTK_LIB="$withvtklib"])

    dnl Initialize Makefile/config.h substitution variables
    VTK_INCLUDE=""
    VTK_LIBRARY=""

    dnl Properly let the substitution variables
    AS_IF([test "$enablevtk" = "yes"],
          [
       dnl Check for existence of a header file in the specified location
       vtkincFound=no;
       AC_CHECK_HEADERS($VTK_INC/vtkConfigure.h, vtkincFound=yes)

       AS_IF([test "$vtkincFound" = "no"],
             [
               AC_MSG_RESULT(VTK header files not found!)
               enablevtk=no;
             ])

       dnl Discover the major, minor, and build versions of VTK by looking in
       dnl vtkConfigure.h and vtkVersionMacros.h
       AS_IF([test "x$enablevtk" = "xyes"],
             [
               dnl If we have the vtkVersionMacros.h (VTK 6.x) find the version there
               AS_IF([test -r $VTK_INC/vtkVersionMacros.h],
                     [
                       vtkmajor=`grep "define VTK_MAJOR_VERSION" $VTK_INC/vtkVersionMacros.h | sed -e "s/.*#define VTK_MAJOR_VERSION[ ]*//g"`
                       vtkminor=`grep "define VTK_MINOR_VERSION" $VTK_INC/vtkVersionMacros.h | sed -e "s/.*#define VTK_MINOR_VERSION[ ]*//g"`
                       vtkbuild=`grep "define VTK_BUILD_VERSION" $VTK_INC/vtkVersionMacros.h | sed -e "s/.*#define VTK_BUILD_VERSION[ ]*//g"`
                     ],
                     dnl Otherwise (VTK 5.x) find the version numbers in vtkConfigure.h
                     [test -r $VTK_INC/vtkConfigure.h],
                     [
                       vtkmajor=`grep "define VTK_MAJOR_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/.*#define VTK_MAJOR_VERSION[ ]*//g"`
                       vtkminor=`grep "define VTK_MINOR_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/.*#define VTK_MINOR_VERSION[ ]*//g"`
                       vtkbuild=`grep "define VTK_BUILD_VERSION" $VTK_INC/vtkConfigure.h | sed -e "s/.*#define VTK_BUILD_VERSION[ ]*//g"`
                     ])
               vtkversion=$vtkmajor.$vtkminor.$vtkbuild
               vtkmajorminor=$vtkmajor.$vtkminor
             ])

       AS_IF([test "x$enablevtk" = "xyes"],
             [
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
         AS_IF([test "x$RPATHFLAG" != "x" && test -d $VTK_LIB],
               [VTK_RPATH_FLAGS="${RPATHFLAG}${VTK_LIB}"])

         dnl VTK 5.x
         AS_IF([test $vtkmajor -eq 5],
               [VTK_LIBRARY="-L$VTK_LIB -lvtkIO -lvtkCommon -lvtkFiltering -lvtkImaging -lvtkParallel"],

               dnl VTK 6.1.x
               dnl Not sure if -lvtkParallelMPI and -lvtkParallelCore existed in VTK-6.1.x
               [test $vtkmajor -eq 6 && test $vtkminor -le 1],
               [
                 VTK_LIBRARY_WITH_VERSION="-L$VTK_LIB -lvtkIOCore-$vtkmajorminor -lvtkCommonCore-$vtkmajorminor -lvtkCommonDataModel-$vtkmajorminor \
                                           -lvtkFiltersCore-$vtkmajorminor -lvtkIOXML-$vtkmajorminor -lvtkImagingCore-$vtkmajorminor \
                                           -lvtkIOImage-$vtkmajorminor -lvtkImagingMath-$vtkmajorminor \
                                           -lvtkParallelMPI-$vtkmajorminor -lvtkParallelCore-$vtkmajorminor"

                 dnl Some Linux distributions (Arch) install VTK without the
                 dnl "libfoo-6.x.so" naming scheme, so we try to handle that
                 dnl situation as well.
                 VTK_LIBRARY_NO_VERSION="-L$VTK_LIB -lvtkIOCore -lvtkCommonCore -lvtkCommonDataModel \
                                         -lvtkFiltersCore -lvtkIOXML -lvtkImagingCore \
                                         -lvtkIOImage -lvtkImagingMath \
                                         -lvtkParallelMPI -lvtkParallelCore"
               ],
               [
                 dnl VTK 6.2.x and above
                 dnl Not sure if -lvtkParallelMPI and -lvtkParallelCore existed in VTK-6.2.x, but it does in VTK-7.x.
                 VTK_LIBRARY_WITH_VERSION="-L$VTK_LIB -lvtkIOCore-$vtkmajorminor -lvtkCommonCore-$vtkmajorminor -lvtkCommonDataModel-$vtkmajorminor \
                                           -lvtkFiltersCore-$vtkmajorminor -lvtkIOXML-$vtkmajorminor -lvtkImagingCore-$vtkmajorminor \
                                           -lvtkIOImage-$vtkmajorminor -lvtkImagingMath-$vtkmajorminor -lvtkIOParallelXML-$vtkmajorminor \
                                           -lvtkParallelMPI-$vtkmajorminor -lvtkParallelCore-$vtkmajorminor"

                 dnl Some Linux distributions (Arch) install VTK without the
                 dnl "libfoo-6.x.so" naming scheme, so we try to handle that
                 dnl situation as well.
                 VTK_LIBRARY_NO_VERSION="-L$VTK_LIB -lvtkIOCore -lvtkCommonCore -lvtkCommonDataModel \
                                         -lvtkFiltersCore -lvtkIOXML -lvtkImagingCore \
                                         -lvtkIOImage -lvtkImagingMath -lvtkIOParallelXML \
                                         -lvtkParallelMPI -lvtkParallelCore"
               ])

         dnl Try to compile test prog to check for existence of VTK libraries.
         dnl AC_LINK_IFELSE uses the LIBS variable.
         AS_IF([test $vtkmajor -gt 5],
               [
                 CPPFLAGS="$CPPFLAGS -I$VTK_INC"

                 dnl 1. First try linking the NO_VERSION library names
                 LIBS="$old_LIBS $VTK_RPATH_FLAGS $VTK_LIBRARY_NO_VERSION"
                 AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_VTK_preamble],[_AX_CXX_COMPILE_VTK_body])],
                                [enablevtk=yes
                                 VTK_LIBRARY=$VTK_LIBRARY_NO_VERSION],
                                [enablevtk=no])

                 dnl 2. If that didn't work, try linking the library names with version numbers.
                 AS_IF([test "x$enablevtk" = "xno"],
                       [
                         LIBS="$old_LIBS $VTK_RPATH_FLAGS $VTK_LIBRARY_WITH_VERSION"
                         AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_VTK_preamble],[_AX_CXX_COMPILE_VTK_body])],
                                        [enablevtk=yes
                                         VTK_LIBRARY=$VTK_LIBRARY_WITH_VERSION],
                                        [enablevtk=no])
                       ])

                 dnl 3. If that also failed, print a message that we failed.
                 AS_IF([test "x$enablevtk" = "xno"],
                       [
                         AC_MSG_RESULT(<<< Linking a test program against the VTK libraries failed >>>)
                         AC_MSG_RESULT([<<< libMesh requires VTK to be configured with -DVTK_Group_MPI:BOOL=ON >>>])
                       ])

                 dnl Reset $LIBS, $CPPFLAGS
                 LIBS="$old_LIBS"
                 CPPFLAGS="$old_CPPFLAGS"
               ],
               [
                 dnl Check for VTK 5.x headers and libraries
                 LIBS="$old_LIBS $VTK_RPATH_FLAGS $VTK_LIBRARY"
                 CPPFLAGS="$CPPFLAGS -I$VTK_INC"

                 dnl test compiling and linking a test program.
                 dnl AC_LINK_IFELSE (input, [action-if-true], [action-if-false])
                 AC_LINK_IFELSE([AC_LANG_PROGRAM([_AX_CXX_COMPILE_VTK_preamble],[_AX_CXX_COMPILE_VTK_body])],
                                [enablevtk=yes],
                                [enablevtk=no])

                 dnl If linking failed, print a status message
                 AS_IF([test "x$enablevtk" = "xno"],
                       [
                         AC_MSG_RESULT(<<< Linking a test program against the VTK libraries failed >>>)
                         AC_MSG_RESULT([<<< libMesh requires VTK to be configured with -DVTK_USE_PARALLEL:BOOL=ON -DVTK_USE_MPI:BOOL=ON >>>])
                       ])

                 dnl Reset $LIBS, $CPPFLAGS
                 LIBS="$old_LIBS"
                 CPPFLAGS="$old_CPPFLAGS"
               ])
       ])

       dnl If both the header file and the required libs were found, continue.
       AS_IF([test "x$enablevtk" = "xyes"],
             [
               VTK_INCLUDE="-I$VTK_INC"

               dnl add the VTK_LIB to the linker run path, if it is a directory
               AS_IF([test "x$RPATHFLAG" != "x" && test -d $VTK_LIB],
                      [VTK_LIBRARY="${RPATHFLAG}${VTK_LIB} $VTK_LIBRARY"])

               AC_SUBST(vtkversion)
               AC_SUBST(vtkmajor)
               AC_SUBST(vtkbuild)

               AC_DEFINE_UNQUOTED(DETECTED_VTK_VERSION_MAJOR, [$vtkmajor],
                 [VTK's major version number, as detected by vtk.m4])

               AC_DEFINE_UNQUOTED(DETECTED_VTK_VERSION_MINOR, [$vtkminor],
                 [VTK's minor version number, as detected by vtk.m4])

               AC_DEFINE_UNQUOTED(DETECTED_VTK_VERSION_SUBMINOR, [$vtkbuild],
                 [VTK's subminor version number, as detected by vtk.m4])

               AC_DEFINE(HAVE_VTK, 1, [Flag indicating whether the library will be compiled with VTK support])
               AC_MSG_RESULT(<<< Configuring library with VTK version $vtkversion support >>>)
             ],
             [
               VTK_LIBRARY=''
               AC_MSG_RESULT(<<< Configuring library without VTK support >>>)
             ])
    ])
  ])

  dnl Substitute the substitution variables
  AC_SUBST(VTK_INCLUDE)
  AC_SUBST(VTK_LIBRARY)

  dnl If VTK is not enabled, but it *was* required, error out now
  dnl instead of compiling libmesh in an invalid configuration.
  AS_IF([test "$enablevtk" = "no" && test "$vtkrequired" = "yes"],
        dnl We return error code 4 here, since 0 means success and 1 is
        dnl indistinguishable from other errors.  Ideally, all of the
        dnl AC_MSG_ERROR calls in our m4 files would return a different
        dnl error code, but currently this is not enforced.
        [AC_MSG_ERROR([*** VTK was not found, but --enable-vtk-required was specified.], 4)])
])
