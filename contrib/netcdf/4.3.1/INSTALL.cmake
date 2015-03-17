Documentation for building netCDF with CMake
********************************************************

This document describes building the netCDF 4.3 C libraries
using KitWare's CMake utility.  By integrating CMake, we are 
able to  provide cross-platform (Windows as well as Linux/Unix)
configuration and build support.

For instructions on getting the latest version of netCDF,
as well as how to get and build the netCDF dependencies,
please refer to the 'INSTALL' file found in the root 
netCDF source directory.


Getting CMake
*************

CMake, a software configuration and testing tool, is maintained by
Kitware.  CMake is available in many linux package management
systems, as well as the 'macports' package management system for 
OSX. 

CMake may also be downloaded for these platforms, as well as Windows,
from the CMake website at http://www.cmake.org.

Building NetCDF with CMake
**************************

The netCDF-C library and utilities requires 3rd party libraries for
full functionality.

  *  Building with NetCDF-4 and the Remote Data Client
  *  Building NetCDF with Classic Library Only
  *  Building with HDF4 Support
  *  Building with Parallel I/O Support

Note that CMake encourages 'out-of-source-tree' builds, i.e. the
directory used to build netCDF is not the root of the netCDF
file structure.  For example, it is fairly common practice to
create a 'build' directory inside the source directory. The examples
in this file will use syntax which assumes the user is currently
located in c:\netcdf-src-dir\build\.

Building NetCDF with Classic Library Only
*****************************************

It is possible to build the netCDF C libraries and utilities so that
only the netCDF classic and 64-bit offset formats are supported, or
the remote data access client is not built.  (See

  http://www.unidata.ucar.edu/netcdf/docs/netcdf_format.html

for more information about the netCDF format variants.  See the
netCDF-DAP site

  http://opendap.org/netCDF-DAP 

for more information about remote client access to data on OPeNDAP
servers.)

To build without support for the netCDF-4 formats or the additional
netCDF-4 functions, but with remote access, use:

Windows:
	> cmake -DCMAKE_INSTALL_PREFIX=C:\Home\Ed\Local -DENABLE_NETCDF_4=OFF
	> cmake --build .
	> cmake --build . --target RUN_TESTS
	> cmake --build . --target INSTALL

Linux/Unix:	
  	> cmake -DCMAKE_INSTALL_PREFIX=/home/ed/local -DENABLE_NETCDF_4=OFF
   	> make 
	> make test install

(Replace ``/home/ed/local'' with the name of the directory where
netCDF is to be installed.)

Starting with version 4.1.1 the netCDF C libraries and utilities have
supported remote data access, using the OPeNDAP protocols.  To build 
with full support for netCDF-4 APIs and format but without remote
client access, use:

Windows:
	> cmake -DCMAKE_INSTALL_PREFIX=C:\Home\Ed\Local -DENABLE_DAP=OFF
	> cmake --build .
	> cmake --build . --target RUN_TESTS
	> cmake --build . --target INSTALL

Linux/Unix
        > cmake -DCMAKE_INSTALL_PREFIX=/home/ed/local -DENABLE_DAP=OFF
	> make
	> make test install
	

If you get the message that netCDF installed correctly, then you are
done!



Building with HDF4 Support
**************************

The netCDF-4 library can (since version 4.1) read HDF4 data files, if
they were created with the SD (Scientific Data) API. To enable this
feature, use the -DENABLE_HDF4=ON option. The location for the HDF4
header files and library must be set in the CPPFLAGS and LDFLAGS
options.

For HDF4 access to work, the library must be build with netCDF-4
features.

Building with Parallel I/O Support
**********************************

For parallel I/O to work, HDF5 must be installed with
-DENABLE_PARALLEL=ON, and an MPI library (and related libraries) must be
made available to the HDF5 configure. This can be accomplished with
the mpicc wrapper script, in the case of MPICH2 (assuming you are building
within the `netcdf/build` directory).

  CC=mpicc cmake .. -DENABLE_PARALLEL=ON -DCMAKE_INSTALL_PREFIX=/shecky/local_par 
  make check install

If the HDF5 used by netCDF has been built with parallel I/O, then
netCDF will also be built with support for parallel I/O. This allows
parallel I/O access to netCDF-4/HDF5 files.  (See

  http://www.unidata.ucar.edu/netcdf/docs/netcdf_format.html

for more information about the netCDF format variants.)

If parallel I/O access to netCDF classic and 64-bit offset files is
also needed, the parallel-netcdf library should also be installed,
(and the replacement pnetcdf.h at

  ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/pnetcdf.h

must be used). Then configure netCDF with the --enable-pnetcdf flag.

Linking to NetCDF
*****************

For static build, to use netCDF-4 you must link to all the libraries,
netCDF, HDF5, zlib, szip (if used with HDF5 build), and curl (if the
remote access client has not been disabled). This will mean -L options
to your build for the locations of the libraries, and -l (lower-case
L) for the names of the libraries.

For shared builds, only -lnetcdf is needed. All other libraries will
be found automatically.

On Windows, when using Visual Studio and compiling a shared library (.dll),
netcdf will be built with an 'import' library, named netcdf.lib. Visual
Studio projects should link against this import library, instead of linking
against the netcdf.dll file directly.

To specify static libraries with CMake, you must use the 'BUILD_SHARED_LIBS=OFF'
flag when invoking CMake.

 	netcdf/build> cmake .. -D"BUILD_SHARED_LIBS=OFF" 

