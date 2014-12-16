Build Instructions for NetCDF-C using CMake {#netCDF-CMake}
===========================================

[TOC]

# Overview {#cmake_overview}

Starting with netCDF-C 4.3.0, we are happy to announce the inclusion of CMake support.  CMake will allow for building netCDF on a wider range of platforms, include Microsoft Windows with Visual Studio.  CMake support also provides robust unit and regression testing tools.  We will also maintain the standard autotools-based build system in parallel.

In addition to providing new build options for netCDF-C, we will also provide pre-built binary downloads for the shared versions of netCDF for use with Visual Studio.  

		
# Requirements {#cmake_requirements}
The following packages are required to build netCDF-C using CMake.

* netCDF-C Source Code
* CMake version 2.8.9 or greater.
* Optional Requirements:
	* HDF5 Libraries for netCDF4/HDF5 support.
	* libcurl for DAP support.

<center>
<img src="deptree.jpg" height="250px" />
</center>

# The CMake Build Process {#cmake_build}

There are four steps in the Build Process when using CMake

1. Configuration: Before compiling, the software is configured based on the desired options.
2. Building: Once configuration is complete, the libraries are compiled.
3. Testing: Post-build, it is possible to run tests to ensure the functionality of the netCDF-C libraries.
4. Installation: If all tests pass, the libraries can be installed in the location specified during configuration.

For users who prefer pre-built binaries, installation packages are available at \ref winbin

## Configuration {#cmake_configuration}

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

## Installation {#cmake_installation}

Once netCDF has been built and tested, it may be installed using the following commands:

> $ make install

or 

> $ cmake --build [Build Directory] --target install

# See Also {#cmake_see_also}

For further information regarding NetCDF and CMake, see \ref cmake_faq
