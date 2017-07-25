Installing and Using netCDF-C Libraries in Windows {#winbin}
==================================================

\brief NetCDF-C Libraries in a Windows Environment may be used under multiple sets of circumstances.

[TOC]

There are several development environments available for programmers who develop on Windows.

* `Microsoft Visual Studio `
* `MSYS/MinGW`
* `Cygwin`

For developers using `Microsoft Visual Studio`, you may download using the Windows build instructions, or you may download the pre-built netCDF-C libraries from this page.

For developers using `MSYS/MinGW` or `Cygwin`, you may build netCDF-C using the Linux/Unix build instructions.

> For complex builds that include netCDF-4 and/or DAP support this may prove tricky, as it is time consuming to collect all of the dependencies.  In these cases it may be easier to use the pre-built `netcdf` packages provided by the `MSYS` and `Cygwin` environments.

Users who prefer to build the netCDF-C libraries from source in a Windows environment using Microsoft Visual Studio are referred to \ref netCDF-CMake

# Getting pre-built netCDF-C Libraries for Visual Studio {#msvc-prebuilt}

These libraries have been built using Visual Studio 2012.  The downloads are installer packages which contain the netCDF-C libraries and utilities (ncgen, ncgen3, ncdump and nccopy), as well as the associated dependencies.


## Included Dependencies {#msvc-inc-deps}

The included dependencies and versions are as follows:

* `libhdf5`: 1.8.17
* `libcurl`: 7.35.0
* `zlib`:    1.2.8

## Latest Release (netCDF-C 4.4.1.1) {#msvc-latest-release}

Configuration		| 32-bit 						| 64-bit |
:-------------------|:--------							|:-------|
netCDF 3		| [netCDF4.4.1.1-NC3-32.exe][r1]		| [netCDF4.4.1.1-NC3-64.exe][r6]
netCDF3+DAP		| [netCDF4.4.1.1-NC3-DAP-32.exe][r2]	| [netCDF4.4.1.1-NC3-DAP-64.exe][r6]
netCDF4			| [netCDF4.4.1.1-NC4-32.exe][r3]		| [netCDF4.4.1.1-NC4-64.exe][r7]
netCDF4+DAP		| [netCDF4.4.1.1-NC4-DAP-32.exe][r4]	| [netCDF4.4.1.1-NC4-DAP-64.exe][r8]

# Using the netCDF-C Libraries with Visual Studio {#msvc-using}

In order to use the netcdf libraries, you must ensure that the .dll files (along with any dependencies from deps/shared/bin) are on the system path. In order to compile a program using these libraries, you must first link your program against the appropriate 'import' (.lib) libraries.

## Install Hierarchy {#msvc-install-hierarchy}

When installed, the netCDF libraries are placed in the specified locations, along with the netCDF-C utilities and dependencies.

<center>
<IMG SRC="InstallTreeWindows.png" width="400"/>
</center>

# Notes {#msvc-notes}

*The following points should be considered when using the netCDF-C libraries on Windows.*

1. When building the netCDF-C libraries with netCDF4 support, using the `Debug` libraries may cause extraneous warnings. These warnings are related to cross-dll memory management, and appear to be harmless. You can safely ignore them by using the `Release` libraries. [NCF-220]


[r1]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC3-32.exe
[r2]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC3-DAP-32.exe
[r3]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC4-32.exe
[r4]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC4-DAP-32.exe
[r6]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC3-64.exe
[r6]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC3-DAP-64.exe
[r7]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC4-64.exe
[r8]: http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.4.1.1-NC4-DAP-64.exe
