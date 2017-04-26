Getting and Building netCDF {#getting_and_building_netcdf}
=============================

[TOC]

This document is for getting and building the netCDF C library and utilities for the most recent released version.  Other libraries that depend on the netCDF C library, such as the Fortran, Python, Java, and C++ libraries, are available as separate distributions that can be optionally built and installed after the C library is successfully installed.  The netCDF-Java library is independent of the netCDF C library unless writing netCDF-4 files from Java is required.

Getting netCDF-C {#getting}
=========================

* For information regarding the netCDF-Fortran libraries, see \subpage building_netcdf_fortran.
* Functionality to make it easier to build netcdf-fortran as part of   the netcdf-c build for *non-MSVC* builds may be enabled at configure time by using the following **Highly Experimental** options:

 * Autotools: `--enable-remote-fortran-bootstrap`
 * CMake:  `-DENABLE_REMOTE_FORTRAN_BOOTSTRAP=ON`

For more details, see the <a href="https://github.com/Unidata/netcdf-c/blob/master/RELEASE_NOTES.md">draft instructions</a> in the Release Notes under the `4.3.3-rc3` section.

Getting pre-built netCDF-C libraries. {#sec_get_pre_built}
-------------------------------------

The easiest way to get netCDF is through a package management program,
such as rpm, yum, homebrew, macports, adept, and others. NetCDF is
available from many different repositories, including the default Red
Hat and Ubuntu repositories.

When getting netCDF from a software repository, you should get a
development version that includes the netcdf.h header file. A
development version will typically have a name such as "netcdf-devel"
or "libnetcdf-dev".

Instructions for installing and using pre-built libraries for Windows may be found here: \ref winbin.

Getting the latest netCDF-C Source Code {#sec_get_source}
----------------------------------------

The netCDF-C source code is hosted from the <a href="http://github.com/Unidata/netcdf-c" >Unidata GitHub repository</a>.


Two options are available for building from source:

* The latest release.
* The developer snapshot.

### The latest release {#sec_latest_release}

The latest full release may be <a href="http://github.com/Unidata/netcdf-c/releases" >downloaded from GitHub</a>.

Source files are available in `.tar.gz` and `.zip` formats.

### The developer snapshot {#sec_dev_snapshot}

The developer snapshot may be cloned from GitHub directly by using the `git` command.

> $ git clone http://github.com/Unidata/netcdf-c netcdf-c

**Note:**

*The developer snapshot release contains bug-fixes and new features added since the last full release, but may also contain new bugs, as it is not tested as extensively as the full release.*

Building netCDF-C {#building}
===========================

The netCDF-C library and utilities require third-party libraries for
full functionality. (See \ref architecture).
* \ref build_default
* \ref build_classic
* \ref build_hdf4
* \ref build_parallel
* \ref building_netcdf_fortran
* \ref configure_options

Requirements {#netcdf_requirements}
----------------------------------

* HDF5 1.8.9 or later (for netCDF-4 support)
* zlib 1.2.5 or later (for netCDF-4 compression)
* curl 7.18.0 or later (for DAP remote access client support)


CMake and Windows support {#sub}
--------------------------------

* \ref netCDF-CMake
* \subpage winbin

Building with netCDF-4 and the Remote Data Client {#build_default}
--------------------------------

The usual way of building netCDF requires the HDF5, zlib, and curl libraries. Versions required are at least HDF5 1.8.9, zlib 1.2.5, and curl 7.18.0 or later.

HDF5 and zlib packages are available from the <a href="http://www.hdfgroup.org/downloads/">HDF5 downloads site</a> and the <a href="http://www.zlib.net/">zlib home site</a>. If you wish to use the remote data client code, then you will also need libcurl, which can be obtained from the <a href="http://curl.haxx.se/download.html">curl website</a>.

> Note that for building netCDF, it is not necessary to build the HDF5 Fortran, C++, or Java API's. Only the HDF5 C library is used, even for netCDF Fortran or C++ libraries.

### Optional: szip support {#op_szip_support}

*Optionally*, you can also build netCDF-4 with the szip library (a.k.a. szlib). If building with szlib, get szip 2.0 or later. NetCDF cannot create szipped data files, but can read HDF5 data files that have used szip. To determine whether license restrictions on the use of szip apply to your situation, see the <a href="http://www.hdfgroup.org/doc_resource/SZIP/">web page on szip compression in HDF products</a>.

If `make check` fails for either `zlib` or `HDF5`, the problem must be resolved before the netCDF-4 installation can continue. For HDF5 problems, see the <a href="http://www.hdfgroup.org/services/support.html">HDF5 help services</a>.

### Building zlib from source {#build_zlib_from_source}

To build zlib from source, specify where you want to install zlib in a shell variable you will also use later (ZDIR, for example), and build it like this from the top-level zlib source directory

~~~~{.py}
    $ # Build and install zlib
    $ ZDIR=/usr/local
    $ ./configure --prefix=${ZDIR}
    $ make check
    $ make install   # or sudo make install, if root permissions required
~~~~

### Building hdf5 from source {#build_hdf5_from_source}

Next, specify where you want to install HDF5 in another shell variable, for example H5DIR, and build it from the HDF5 top-level source directory:

~~~~{.py}
    $ # Build and install HDF5
    $ H5DIR=/usr/local
    $ ./configure --with-zlib=${ZDIR} --prefix=${H5DIR}
    $ make check
    $ make install   # or sudo make install, if root permissions required
~~~~

If you are building HDF5 with the optional szip library, include the `--with-szlib=` option to specify where it was installed.

In all cases, the installation location specified with the `--prefix` option *must be different* from the source directory where the software is being built.

### Building netCDF-4 and the Remote Data Client from source {#build_nc4_dap_from_source}

Before building netCDF, you may need to add `${H5DIR}/lib` to the LD_LIBRARY_PATH environment variable if that lib directory is not searched by default. See the <a href="http://www.unidata.ucar.edu/netcdf/docs/faq.html#Shared%20Libraries">netCDF FAQ</a> for more details on using shared libraries.

Indicate where you want to install netCDF in another shell variable, for example NCDIR. Then run the netCDF configure script, specifying where HDF5 was installed using the CPPFLAGS and LDFLAGS environment variables. For example, from the top-level netCDF source directory:

~~~~{.py}
    $ # Build and install netCDF-4
    $ NCDIR=/usr/local
    $ CPPFLAGS=-I${H5DIR}/include LDFLAGS=-L${H5DIR}/lib ./configure --prefix=${NCDIR}
    $ make check
    $ make install  # or sudo make install
~~~~

If you don't provide a `--prefix` option, installation will be in `/usr/local/`, in subdirectories lib/, include/, and bin/.  The installation location specified with the `--prefix` option must be different from the source directory where the software is being built.

> WARNING: you should be able to use parallel 'make all'. But 'make check' will probably fail if you use parallel make. This is because historically, there are inter-dependencies between test programs. It is unlikely that this will be fixed any time soon, if ever.

Building netCDF with Classic Library Only {#build_classic}
---------------------------------------

It is possible to build the netCDF C libraries and utilities so that
only the netCDF classic and 64-bit offset formats are supported, or
the remote data access client is not built. (See \ref netcdf_format
for more information about the netCDF format variants.  See the <a
href="http://www.opendap.org/documentation">DAP documentation and
support site</a> for more information about remote client access to
data on OPeNDAP servers.)

If necessary, set the NCDIR shell variable to indicate where netCDF should be
installed. Then to build a netCDF-3 library without support for the
netCDF-4 formats or functions, but with remote client access, use:

~~~~{.py}
    $ # Build and install netCDF-3 from netCDF-4 source
    $ ./configure --prefix=${NCDIR} --disable-netcdf-4
    $ make check install
~~~~

To build with full support for netCDF-4 API's and format but without
remote client access, use:

~~~~{.py}
    $ # Build and install netCDF-4 without DAP client support
    $ ./configure --prefix=${NCDIR} --disable-dap
    $ make check install
~~~~

To build without netCDF-4 support or remote client access, use:

~~~~{.py}
    $ # Build and install minimal netCDF-3 with no DAP client support
    $ ./configure --prefix=${NCDIR} --disable-netcdf-4 --disable-dap
    $ make check install
~~~~

If you get the message that netCDF installed correctly, then you are
done!

Building with HDF4 Support {#build_hdf4}
---------------------

The netCDF-4 library can read HDF4 data files, if they were created
with the SD (Scientific Data) API.

For this to work, you must build the HDF4 library with the
configure option `--disable-netcdf`
to prevent it from building an HDF4 version of the netCDF-2 library
that conflicts with the netCDF-2 functions that are built into the Unidata
netCDF library.

Then, when building netCDF-4, use the `--enable-hdf4`.
option to configure. The location for the HDF4 header files and
library must be specified in the CPPFLAGS and LDFLAGS environment variables
or configure options.

For HDF4 access to work, the library must be built with netCDF-4
features.

Here's an example, assuming the HDF5 library has been built and
installed in H5DIR and you will build and install the HDF4 library in
H4DIR (which could be the same as H5DIR). From the top-level HDF4
source directory:

~~~~{.py}
    $ # Build and install HDF4
    $ ./configure --enable-shared --disable-netcdf --disable-fortran --prefix=${H4DIR}
    $ make check
    $ make install
~~~~

Then from the top-level netCDF directory:

~~~~{.py}
    $ # Build and install netCDF-4 with HDF4 access enabled
    $ CPPFLAGS="-I${H5DIR}/include -I${H4DIR}/include" \
      LDFLAGS="-L${H5DIR}/lib -L${H4DIR}/lib" \
      ./configure --enable-hdf4 --enable-hdf4-file-tests
    $ make check
    $ make install
~~~~

Building with Parallel I/O Support {#build_parallel}
--------------

For parallel I/O to work, HDF5 must be installed with
–enable-parallel, and an MPI library (and related libraries) must be
made available to the HDF5 configure. This can be accomplished with
an mpicc wrapper script.

The following works from the top-level HDF5 source directory to build
HDF5 with parallel I/O:

~~~~{.py}
    $ # Build and install HDF5 with parallel support
    $ CC=mpicc ./configure --enable-parallel --prefix=${H5DIR}
    $ make check
    $ make install
~~~~

If the HDF5 used by netCDF has been built with parallel I/O, then netCDF will also be built with inherited support for parallel I/O. This allows parallel I/O access to netCDF-4/HDF5 files.  (See /ref netcdf_formats for more information about the netCDF format variants.)

From the top-level netCDF-4 source directory, the following builds netCDF-4 with parallel I/O, assuming H5DIR specifies where parallel HDF5 was installed:

~~~~{.py}
    $ # Build, test, and install netCDF-4 with HDF5 parallel support
    $ CC=mpicc CPPFLAGS=-I${H5DIR}/include LDFLAGS=-L${H5DIR}/lib \
      ./configure --disable-shared --enable-parallel-tests --prefix=${NCDIR}
    $ make check
    $ make install
~~~~

If parallel I/O access to netCDF classic and 64-bit offset files is
needed, an alternate
[parallel-netcdf library](https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/WikiStart),
referred to as "PnetCDF", must also be installed. Assume it was
installed in the directory named by the PNDIR shell variable.
Then, from the top-level netCDF-4 source directory, configure netCDF
with the "--enable-pnetcdf" option:

~~~~{.py}
    $ # Build, test, and install netCDF-4 with pnetcdf support
    $ CC=mpicc CPPFLAGS="-I${H5DIR}/include -I${PNDIR}/include" \
      LDFLAGS="-L${H5DIR}/lib -L${PNDIR}/lib" ./configure \
	  --disable-shared --enable-pnetcdf  --enable-parallel-tests \
	  --prefix=${NCDIR}
    $ make check
    $ make install
~~~~

Linking to netCDF-C {#linking}
-------------------

For static builds of applications that use netCDF-4 you must link to all the libraries,
netCDF, HDF5, zlib, szip (if used with HDF5 build), and curl (if the
remote access client has not been disabled). This will require -L options
to your build for the locations of the libraries, and -l (lower-case
L) for the names of the libraries.

For example, you might build other applications with netCDF-4 by
setting the LIBS environment variable, assuming NCDIR, H5DIR, and ZDIR
indicate where netCDF, HDF5, and zlib are installed:

~~~~{.py}
    LIBS="-L${NCDIR}/lib -lnetcdf -L${H5DIR}/lib -lhdf5_hl -lhdf5 -L${ZDIR}/lib -lz -lm"
~~~~

For shared builds, only `-L${NCDIR}/lib -lnetcdf` is
needed. All other libraries will be found automatically.

The `pkg-config` or `nc-config` utilities can be
used to specify build options for software that uses netCDF.  For
example, to compile and link an application named myapp.c with a
netCDF-C libraries, whether shared or static, you can use

~~~~{.py}
    $ cc -o myapp myapp.c `nc-config --cflags --libs`
~~~~

or

~~~~{.py}
    $ PKG_CONFIG_PATH=${NCDIR}/lib/pkgconfig:$PKG_CONFIG_PATH
    $ export PKG_CONFIG_PATH
    $ cc -o myapp myapp.c `pkg-config --cflags --libs netcdf`
~~~~

configure options {#configure_options}
-----------------------------

These options are used for `autotools`-based builds.yup

Note: `--disable` prefix indicates that the option is normally enabled.
<table>
<tr><th>Option<th>Description<th>Dependencies
<tr><td>--disable-doxygen<td>Disable generation of documentation.<td>doxygen
<tr><td>--disable-fsync<td>disable fsync support<td>kernel fsync support

<tr><td>--disable-netcdf-4<td>build netcdf-3 without HDF5 and zlib<td>
<tr><td>--disable-netcdf4<td>synonym for disable-netcdf-4
<tr><td>--enable-hdf4<td>build netcdf-4 with HDF4 read capability<td>HDF4, HDF5 and zlib
<tr><td>--enable-hdf4-file-tests<td>test ability to read HDF4 files<td>selected HDF4 files from Unidata ftp site
<tr><td>--enable-pnetcdf<td>build netcdf-4 with parallel I/O for classic and
                          64-bit offset files using parallel-netcdf
<tr><td>--enable-extra-example-tests<td>Run extra example tests<td>--enable-netcdf-4,GNU sed
<tr><td>--enable-parallel-tests <td>run extra parallel IO tests<td>--enable-netcdf-4, parallel IO support
<tr><td>--enable-logging<td>enable logging capability<td>--enable-netcdf-4
<tr><td>--disable-dap<td>build without DAP client support.<td>libcurl
<tr><td>--disable-dap-remote-tests<td>disable dap remote tests<td>--enable-dap
<tr><td>--enable-dap-long-tests<td>enable dap long tests<td>
<tr><td>--enable-extra-tests<td>run some extra tests that may not pass because of known issues<td>
<tr><td>--enable-ffio<td>use ffio instead of posixio (ex. on the Cray)<td>
<tr><td>--disable-examples<td>don't build the netCDF examples during make check
                          (examples are treated as extra tests by netCDF)<td>
<tr><td>--disable-v2<td>turn off the netCDF version 2 API<td>
<tr><td>--disable-utilities<td>don't build netCDF utilities ncgen, ncdump, and nccopy<td>
<tr><td>--disable-testsets<td>don't build or run netCDF tests<td>
<tr><td>--enable-large-file-tests <td>Run tests which create very large data
		          files<td>~13 GB disk space required, but recovered when
                          tests are complete). See option --with-temp-large to
                          specify temporary directory
<tr><td>--enable-benchmarks<td>Run benchmarks. This is an experimental feature.
			  The benchmarks are extra tests, used to check netCDF performance.
    <td>sample data files from the Unidata ftp site
<tr><td>--disable-extreme-numbers
<td>don't use extreme numbers during testing, such as MAX_INT - 1<td>
<tr><td>--disable-shared<td>don't build shared libraries<td>
<tr><td>--disable-static<td>don't build static libraries<td>
<tr><td>--disable-largefile<td>omit support for files larger than 2GB<td>
<tr><td>--enable-mmap<td>Use mmap to implement NC_DISKLESS<td>System-provided `mmap` or `mremap` functions
<tr><td>--enable-valgrind-tests <td>build with valgrind-tests; static builds only<td>valgrind
</table>

Build Instructions for netCDF-C using CMake {#netCDF-CMake}
===========================================

## Overview {#cmake_overview}

Starting with netCDF-C 4.3.0, we are happy to announce the inclusion of CMake support.  CMake will allow for building netCDF on a wider range of platforms, include Microsoft Windows with Visual Studio.  CMake support also provides robust unit and regression testing tools.  We will also maintain the standard autotools-based build system in parallel.

In addition to providing new build options for netCDF-C, we will also provide pre-built binary downloads for the shared versions of netCDF for use with Visual Studio.


##  Requirements {#cmake_requirements}
The following packages are required to build netCDF-C using CMake.

* netCDF-C Source Code
* CMake version 2.8.12 or greater.
* Optional Requirements:
	* HDF5 Libraries for netCDF4/HDF5 support.
	* libcurl for DAP support.

<center>
<img src="deptree.jpg" height="250px" />
</center>

## The CMake Build Process {#cmake_build}

There are four steps in the Build Process when using CMake

1. Configuration: Before compiling, the software is configured based on the desired options.
2. Building: Once configuration is complete, the libraries are compiled.
3. Testing: Post-build, it is possible to run tests to ensure the functionality of the netCDF-C libraries.
4. Installation: If all tests pass, the libraries can be installed in the location specified during configuration.

For users who prefer pre-built binaries, installation packages are available at \ref winbin

### Configuration {#cmake_configuration}

The output of the configuration step is a project file based on the appropriate configurator specified.  Common configurators include:

* Unix Makefiles
* Visual Studio
* CodeBlocks
* ... and others

### Common CMake Options {#cmake_common_options}

| **Option** | **Autotools** | **CMake** |
| :------- | :---- | :----- |
Specify Install Location | --prefix=PREFIX | -D"CMAKE\_INSTALL\_PREFIX=PREFIX"
Enable/Disable netCDF-4 | --enable-netcdf-4<br>--disable-netcdf-4 | -D"ENABLE\_NETCDF\_4=ON" <br> -D"ENABLE\_NETCDF\_4=OFF"
Enable/Disable DAP | --enable-dap <br> --disable-dap | -D"ENABLE\_DAP=ON" <br> -D"ENABLE\_DAP=OFF"
Enable/Disable Utilities | --enable-utilities <br> --disable-utilities | -D"BUILD\_UTILITIES=ON" <br> -D"BUILD\_UTILITIES=OFF"
Specify shared/Static Libraries | --enable-shared <br> --enable-static | -D"BUILD\_SHARED\_LIBS=ON" <br> -D"BUILD\_SHARED\_LIBS=OFF"
Enable/Disable Tests | --enable-testsets <br> --disable-testsets | -D"ENABLE\_TESTS=ON" <br> -D"ENABLE\_TESTS=OFF"
Specify a custom library location | Use *CFLAGS* and *LDFLAGS* | -D"CMAKE\_PREFIX\_PATH=/usr/custom_libs/"

A full list of *basic* options can be found by invoking `cmake [Source Directory] -L`. To enable a list of *basic* and *advanced* options, one would invoke `cmake [Source Directory] -LA`.

### Configuring your build from the command line. {#cmake_command_line}

The easiest configuration case would be one in which all of the dependent libraries are installed on the system path (in either Unix/Linux or Windows) and all the default options are desired. From the build directory (often, but not required to be located within the source directory):

> $ cmake [Source Directory]

If you have libraries installed in a custom directory, you may need to specify the **CMAKE\_PREFIX_PATH** variable to tell cmake where the libraries are installed. For example:

> $ cmake [Source Directory] -DCMAKE\_PREFIX\_PATH=/usr/custom_libraries/

## Building {#cmake_building}

The compiler can be executed directly with 'make' or the appropriate command for the configurator which was used.

> $ make

Building can also be executed indirectly via cmake:

> $ cmake --build [Build Directory]

## Testing {#cmake_testing}

Testing can be executed several different ways:

> $ make test

or

> $ ctest

or

> $ cmake --build [Build Directory] --target test

### Installation {#cmake_installation}

Once netCDF has been built and tested, it may be installed using the following commands:

> $ make install

or

> $ cmake --build [Build Directory] --target install

## See Also {#cmake_see_also}

For further information regarding netCDF and CMake, see \ref cmake_faq
