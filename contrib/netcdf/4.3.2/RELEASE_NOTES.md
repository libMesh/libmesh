Release Notes {#release_notes}
===============================

\brief Release notes file for the netcdf-c package.

This file contains a high-level description of this package's evolution. Releases are in reverse chronological order (most recent first). Note that, as of netcdf 4.2, the netcdf-c++ and netcdf-fortran libraries have been separated into their own libraries.

## 4.3.2 Released 2014-04-23

* As part of an ongoing project, the Doxygen-generated netcdf documentation has been reorganized.  The goal is to make the documentation easier to parse, and to eliminate redundant material.  This project is ongoing.

* The oc .dodsrc reader was improperly handling the user name and password entries. [NCF-299](https://bugtracking.unidata.ucar.edu/browse/NCF-299)

* CTestConfig.cmake has been made into a template so that users may easily specify the location of an alternative CDash-based Dashboard using the following two options:

	* `NC_TEST_DROP_SITE` - Specify an alternative Dashboard by URL or IP address.

	* `NC_CTEST_DROP_LOC_PREFIX` - Specify a prefix on the remote webserver relative to the root directory. This lets CTest accommodate dashboards that do not live at the top level of the web server.
	
* Return an error code on open instead of an assertion violation for truncated file.

### 4.3.2-rc2 Released 2014-04-15

* Cleaned up a number of CMake inconsistencies related to CMake usage, parallel builds.
* Now passing -Wl,--no-undefined to linker when appropriate.
* Corrected an issue preventing large file tests from running correctly under Windows.
* Misc Bug Fixes detected by static analysis.

### 4.3.2-rc1 Released 2014-03-20

* Pre-built Windows downloads will now be bundled with the latest (as of the time of this writing) versions of the various dependencies:
	* `hdf5: 1.8.12`
	* `zlib: 1.2.8`
	* `libcurl: 7.35.0`
	
* Added a separate flag to enable DAP AUTH tests. These tests are disabled by default.  The flags for autotools and CMAKE-based builds are (respectively):
	* --enable-dap-auth-tests
	* -DENABLE\_DAP\_AUTH\_TESTS

* Fixed small default chunk size for 1-dimensional record variables.  [NCF-211](https://bugtracking.unidata.ucar.edu/browse/NCF-211)

* Cleaned up type handling in netCDF-4 to fix bugs with fill-values.

* Corrected "BAIL" macros to avoid infinite loop when logging is disabled and an error occurs.

* Refactored how types are used for attributes, variables, and committed types, clarifying and categorizing fields in structs, and eliminating duplicated type information between variables and types they use.

* Made type structure sharable by committed datatypes and variables that use it.

* Handled string datatypes correctly, particularly for fill value attributes. Expanded testing for string fill values.

* Simplified iteration of objects in the file when it's opened, tracking fewer objects and using less memory.

* Enabled netCDF-4 bit-for-bit reproducibility for nccopy and other applications (thanks to Rimvydas Jasinskas and Quincey Koziol) by turning off HDF5 object creation, access, and modification time tracking.  [NCF-290](https://bugtracking.unidata.ucar.edu/browse/NCF-290)

* Addressed an issue where `cmake`-based builds would not properly create a `pkg-config` file. This file is now created properly by `cmake`.  [NCF-288](https://bugtracking.unidata.ucar.edu/browse/NCF-288)

* Addressed an issue related to old DAP servers. [NCF-287](https://bugtracking.unidata.ucar.edu/browse/NCF-287)

## 4.3.1.1 Released 2014-02-05

This is a bug-fix-only release for version 4.3.1.

* Corrected a DAP issue reported by Jeff Whitaker related to non-conforming servers.

* Corrected an issue with DAP tests failing in a 64-bit Cygwin environment. [NCF-286](https://bugtracking.unidata.ucar.edu/browse/NCF-286)

## 4.3.1 Released 2014-01-16

* Add an extended format inquiry method to the netCDF API: nc\_inq\_format\_extended. NC\_HAVE\_INQ\_FORMAT\_EXTENDED is defined in netcdf.h [NCF-273]

[NCF-273]:https://bugtracking.unidata.ucar.edu/browse/NCF-273


### 4.3.1-rc6 Released 2013-12-19

* Fixed fill value handling for string types in nc4\_get\_vara().

* Corrected behavior of nc\_inq\_unlimdim and nv\_inq\_unlimdims to report dimids
  in same order as nc\_inq\_dimids.

* Addressed an issue reported by Jeff Whitaker regarding `nc_inq_nvars` returning an incorrect number of dimensions (this issue was introduced in 4.3.1-rc5).  Integrated a test contributed by Jeff Whitaker.

* A number of previously-disabled unit tests were reviewed and made active. 


### 4.3.1-rc5 Released 2013-12-06

* When opening a netCDF-4 file, streamline the iteration over objects in the underlying HDF5 file.

* Fixed netCDF-4 failure when renaming a dimension and renaming a variable using that dimension, in either order. [NCF-177]

[NCF-177]:https://bugtracking.unidata.ucar.edu/browse/NCF-177

* When compiling with `hdf4` support, both autotools and cmake-based builds now properly look for the `libjpeg` dependency and will link against it when found (or complain if it's not).  Also added `ENABLE_HDF4_FILE_TESTS` option to CMake-based builds.

* Fixed bug in ncgen; it was not properly filling empty string constants ("") to be the proper length. [NCF-279]

[NCF-279]:https://bugtracking.unidata.ucar.edu/browse/NCF-279

* Fixed bug in ncgen where it was interpreting int64 constants
  as uint64 constants. [NCF-278]

[NCF-278]:https://bugtracking.unidata.ucar.edu/browse/NCF-278

* Fixed bug in handling Http Basic Authorization. The code was actually there but was not being executed. [NCF-277]

[NCF-277]:https://bugtracking.unidata.ucar.edu/browse/NCF-277

* Added hack to the DAP code to address a problem with the Columbia.edu server. That server does not serve up proper DAP2 DDS replies. The Dataset {...} name changes depending on if the request has certain kinds of constraints. [NCF-276]

[NCF-276]:https://bugtracking.unidata.ucar.edu/browse/NCF-276

* Fixed bugs with ncdump annotation of values, using -b or -f
  options. [NCF-275]

[NCF-275]:https://bugtracking.unidata.ucar.edu/browse/NCF-275


### 4.3.1-rc4 Released 2013-11-06

* Addressed an issue on Windows where `fstat` would report an incorrect file size on files > 4GB.  [NCF-219]


* Added better documentation about accessing ESG datasets.
  See http://www.unidata.ucar.edu/software/netcdf/docs/esg.html.

* Corrected an issue with CMake-based builds enabling HDF4 support where the HDF4 libraries were in a non-standard location.

* Fix bug introduced by [NCF-267] where octal constants above
'\177' were not recognized as proper octal constants. [NCF-271]

[NCF-271]:https://bugtracking.unidata.ucar.edu/browse/NCF-271

* Fixed an issue where the `netcdf.3` man page was not being installed by CMake-based builds. [Github](https://github.com/Unidata/netcdf-c/issues/3)



### 4.3.1-rc3 Released 2013-09-24

* Modify ncgen to support NUL characters in character array
  constants. [NCF-267]

[NCF-267]:https://bugtracking.unidata.ucar.edu/browse/NCF-267

* Modify ncgen to support disambiguating references to
  an enum constant in a data list. [NCF-265]
  
[NCF-265]:https://bugtracking.unidata.ucar.edu/browse/NCF-265

* Corrected bug in netCDF-4 dimension ID ordering assumptions, resulting in access that works locally but fails through DAP server. [NCF-166]

[NCF-166]:https://bugtracking.unidata.ucar.edu/browse/NCF-166

* Added a new configuration flag, `NC_USE_STATIC_CRT` for CMake-based Windows builds.  The default value is 'OFF'.  This will allow the user to define whether to use the shared CRT libraries (\\MD) or static CRT libraries (\\MT) in Visual Studio builds.

* Ensure netCDF-4 compiles with OpenMPI as an alternative to MPICH2. [NCF-160]

[NCF-160]:https://bugtracking.unidata.ucar.edu/browse/NCF-160

* Addressed issue with hanging Parallel netCDF-4 using HDF5 1.8.10. [NCF-240]

[NCF-240]:https://bugtracking.unidata.ucar.edu/browse/NCF-240

* Addressed issue with Large File Support on Windows, using both 32 and 64-bit builds. [NCF-219]

[NCF-219]:https://bugtracking.unidata.ucar.edu/browse/NCF-219

* Removed deprecated directories:
	* librpc/
	* udunits/
	* libcf/
	* libcdmr/

### 4.3.1-rc2 Released 2013-08-19

* Added `configure` and accompanying configuration files/templates to release repository.  **These will only be added to tagged releases on GitHub**.

* Integrated a fix by Quincey Koziol which addressed a variation of [NCF-250], *Fix issue of netCDF-4 parallel independent access with unlimited dimension hanging*.

[NCF-250]:https://www.unidata.ucar.edu/jira/browse/NCF-250

* Integrated change contributed by Orion Poplawski which integrated GNUInstallDirs into the netCDF-C CMake system; this will permit systems that install into lib64 (such as Fedora) to `make install` without problem.

* Corrected an error with the CMake config files that resulted in the `netcdf.3` manpage not being built or installed.

### 4.3.1-rc1 Released 2013-08-09

* Migrated from the netCDF-C `subversion` repository to a publically available GitHub repository available at https://github.com/Unidata/netCDF-C.  This repository may be checked out (cloned) with the following command:

	$ git clone https://github.com/Unidata/netCDF-C.git

* Note: in this release, it is necessary to generate the `configure` script and makefile templates using `autoreconf` in the root netCDF-C directory.:
	
	$ autoreconf -i -f 

* Added `nc_rename_grp` to allow for group renaming in netCDF-4 files. [NCF-204]

[NCF-204]: https://bugtracking.unidata.ucar.edu/browse/NCF-204

* Added a `NC_HAVE_RENAME_GRP` macro to netcdf.h, [as per a request by Charlie Zender][cz1]. This will allow software compiling against netcdf to easily query whether or not nc\_rename\_grp() is available.

[cz1]: https://www.unidata.ucar.edu/esupport/staff/index.php?_m=tickets&_a=viewticket&ticketid=22442

* Added Greg Sjaardema's contributed optimization for the nc4\_find\_dim\_len function in libsrc4/nc4internal.c. The patch eliminates several malloc/free calls that exist in the original coding.

* Added support for dynamic loading, to compliment the dynamic loading support introduced in hdf 1.8.11.  Dynamic loading support depends on libdl, and is enabled as follows: [NCF-258]
	* autotools-based builds: --enable-dynamic-loading
	* cmake-based builds: -DENABLE\_DYNAMIC\_LOADING=ON
	
[NCF-258]: https://www.unidata.ucar.edu/jira/browse/NCF-258

* Fix issue of netCDF-4 parallel independent access with unlimited dimension hanging.  Extending the size of an unlimited dimension in HDF5 must be a collective operation, so now an error is returned if trying to extend in independent access mode. [NCF-250]

[NCF-250]: https://bugtracking.unidata.ucar.edu/browse/NCF-250

* Fixed bug with netCDF-4's inability to read HDF5 scalar numeric attributes. Also allow, in addition to zero length strings, a new NULL pointer as a string value. to improve interoperability with HDF5. This required a new CDL constant, 'NIL', that can be output from ncdump for such a string value in an HDF5 or netCDF-4 file. The ncgen utility was also modified to properly handle such NIL values for strings. [NCF-56]

[NCF-56]: https://bugtracking.unidata.ucar.edu/browse/NCF-56

* Parallel-build portability fixes, particularly for OpenMPI and gcc/gfortran-4.8.x on OSX.

* Fix contributed by Nath Gopalaswamy to large file problem reading netCDF classic or 64-bit offset files that have a UINT32_MAX flag for large last record size of a variable that has values larger than 1 byte.  This problem had previously been fixed for *writing* such data, but was only tested with an ncbyte variable.

* Fixed various minor documentation problems.

## 4.3.0 Released 2013-04-29

* fsync: Changed default in autotools config file; fsync must now be
explicitely enabled instead of explicitely disabled. [NCF-239]

[NCF-239]: https://www.unidata.ucar.edu/jira/browse/NCF-239

* Fixed netCDF-4 bug where odometer code for libdap2 mishandled stride \> 1. Bug reported by Ansley Manke. [NCF-249]

[NCF-249]: https://www.unidata.ucar.edu/jira/browse/NCF-249

* Fixed netCDF-4 bug so netCDF just ignores objects of HDF5 reference type in
the file, instead of rejecting the file. [NCF-29]

[NCF-29]: https://www.unidata.ucar.edu/jira/browse/NCF-29

* Fixed netCDF-4 bug with particular order of creation of dimensions,
coordinate variables, and subgroups resulting in two dimensions with the
same dimension ID. [NCF-244]

[NCF-244]: https://www.unidata.ucar.edu/jira/browse/NCF-244

* Fixed netCDF-4 bug with a multidimensional coordinate variable in a
subgroup getting the wrong dimension IDs for its dimensions. [NCF-247]

[NCF-247]: https://www.unidata.ucar.edu/jira/browse/NCF-247

* Fixed bug with incorrect fixed-size variable offsets in header getting
written when schema changed for files created by parallel-netcdf. Thanks
to Wei-keng Liao for developing and contributing the fix. [NCF-234]

[NCF-234]: https://www.unidata.ucar.edu/jira/browse/NCF-234

* Fixed bug in handling old servers that do not do proper Grid to
Structure conversions. [NCF-232]

[NCF-232]: https://www.unidata.ucar.edu/jira/browse/NCF-232

* Replaced the oc library with oc2.0

* Fix bug with nc\_get\_var1\_uint() not accepting unsigned ints larger
than 2\*\*31. [NCF-226]

[NCF-226]: https://www.unidata.ucar.edu/jira/browse/NCF-226

* Fix to convert occurrences of '/' in DAP names to %2f. [NCF-223]

[NCF-223]: https://www.unidata.ucar.edu/jira/browse/NCF-223

* Fix bug in netCDF-4 with scalar non-coordinate variables with same name
as dimensions. [NCF-222]

[NCF-222]: https://www.unidata.ucar.edu/jira/browse/NCF-222

* Fix bug in which calling netCDF-4 functions in which behavior that
should not depend on order of calls sometimes produces the wrong
results. [NCF-217]

[NCF-217]: https://www.unidata.ucar.edu/jira/browse/NCF-217

* Merged in nccopy additions from Martin van Driel to support -g and -v
options for specifying which groups or variables are to be copied.
[NCF-216]

[NCF-216]: https://www.unidata.ucar.edu/jira/browse/NCF-216

* Merged in parallel-netcdf bugs fixes from Greg Sjaardema. [NCF-214]

[NCF-214]: https://www.unidata.ucar.edu/jira/browse/NCF-214

* Modify ncgen so that if the incoming file has a special attribute, then
it is used to establish the special property of the netcdf file, but the
attribute is not included as a real attribute in the file. [NCF-213].

[NCF-213]: https://www.unidata.ucar.edu/jira/browse/NCF-213

* Added library version info to the user-agent string so that the server
logs will be more informative. [NCF-210]

[NCF-210]: https://www.unidata.ucar.edu/jira/browse/NCF-210

* Added work around for bad servers that sometimes sends DAP dataset with
duplicate field names. [NCF-208]

[NCF-208]: https://www.unidata.ucar.edu/jira/browse/NCF-208

* Fixed bug with strided access for NC\_STRING type. [NCF-206]

[NCF-206]: https://www.unidata.ucar.edu/jira/browse/NCF-206

* Prevented adding an invalid \_FillValue attribute to a variable (with
nonmatching type or multiple values), to avoid later error when any
record variable is extended. [NCF-190]

[NCF-190]: https://www.unidata.ucar.edu/jira/browse/NCF-190

* Fix bug in which some uses of vlen within compounds causes HDF5 errors.
[NCF-155]

[NCF-155]: https://www.unidata.ucar.edu/jira/browse/NCF-155

* Fixed ncdump bug in display of data values of variables that use
multiple unlimited dimensions. [NCF-144]

[NCF-144]: https://www.unidata.ucar.edu/jira/browse/NCF-144

* Fix bug in which interspersing def\_var calls with put\_var calls can
lead to corrupt metadata in a netCDF file with groups and inherited
dimensions. [NCF-134]

[NCF-134]: https://www.unidata.ucar.edu/jira/browse/NCF-134

* Building shared libraries works with DAP and netCDF4 functionality.
[NCF-205] [NCF-57]

[NCF-205]: https://www.unidata.ucar.edu/jira/browse/NCF-205
[NCF-57]: https://www.unidata.ucar.edu/jira/browse/NCF-57

* 32-and-64-bit builds are working under MinGW on Windows. [NCF-112]

[NCF-112]: https://www.unidata.ucar.edu/jira/browse/NCF-112

* Config.h for Windows compiles are included in the build. [NCF-98]

[NCF-98]: https://www.unidata.ucar.edu/jira/browse/NCF-98

* NetCDF-4 dependency on NC\_MAX\_DIMS has been removed. [NCF-71]

[NCF-71]: https://www.unidata.ucar.edu/jira/browse/NCF-71

* 64-bit DLL's are produced on Windows. [NCF-65]

[NCF-65]: https://www.unidata.ucar.edu/jira/browse/NCF-65

* DLL Packaging issues are resolved. [NCF-54]

[NCF-54]: https://www.unidata.ucar.edu/jira/browse/NCF-54

* The CMake build system (with related ctest and cdash systems for
testing) has been integrated into netCDF-C. This allows for Visual
Studio-based builds in addition to gcc-based builds. This requires at
least CMake version 2.8.8. This replaces/supplements the cross-compiled
set of Visual-Studio compatible netCDF libraries introduced in netCDF
4.2.1-rc1.

## 4.2.1.1 Released 2012-08-03

* Patched libdap2/ncdap3.c to fix DAP performance bug remotely accessing large files (> 2GiB).

* Patched ncdump/dumplib.c to properly escape special characters in CDL output from ncdump for netCDF-4 string data.


### 4.2.1 Released 2012-07-18

* Added a specific NC\_MMAP mode flag to modify behavior of NC\_DISKLESS.

* Changed the file protections for NC\_DISKLESS created files to 0666
[NCF-182]

* Fixed ncdump to report error when an unsupported option is specified.
[NCF-180]

* Fixed documentation of CDL char constants in Users Guide and ncgen man
page.

* Fixed memory leak detected by valgrind in one of the HDF5 tests.

* Fixed problem with \#elif directives in posixio.c revealed by PGI
compilers.

### 4.2.1-rc1 Released 2012-06-18

* Ported static and shared libraries (DLL's) for both 32- and 64-bit
Windows, including support for DAP remote access, with netCDF-3 and
netCDF-4/HDF5 support enabled. The environment for this build is
MSYS/MinGW/MinGW64, but the resulting DLLs may be used with Visual
Studio. [NCF-112] [NCF-54] [NCF-57] [NCF-65]

* Implemented diskless files for all netCDF formats. For nc\_create(),
diskless operation performs all operations in memory and then optionally
persists the results to a file on close. For nc\_open(), but only for
netcdf classic files, diskless operation caches the file in-memory,
performs all operations on the memory resident version and then writes
all changes back to the original file on close.
[NCF-110][NCF-109][NCF-5]

* Added MMAP support. If diskless file support is enabled, then it is
possible to enable implementation of diskless files using the operating
system's MMAP facility (if available). The enabling flag is
"--enable-mmap". This is most useful when using nc\_open() and when only
parts of files, a single variable say, need to be read.

* Added configure flag for --disable-diskless.

* Added nccopy command-line options to exploit diskless files, resulting
in large speedups for some operations, for example converting unlimited
dimension to fixed size or rechunking files for faster access. Upgraded
doxygen and man-page documentation for ncdump and nccopy utilities,
including new -w option for diskless nccopy, with an example.

* Modified Makefile to allow for concurrent builds and to support builds
outside the source tree, e.g. 'mkdir build; cd build;
SOURCE-DIR/configure' where SOURCE-DIR is the top-level source
directory.

* Fixed some netCDF-4 bugs with handling strings in non-netCDF-4 HDF5
files. [NCF-150]

* Fixed bug using nccopy to compress with shuffling that doesn't compress
output variables unless they were already compressed in the input file.
[NCF-162]

* Fixed bug in 64-bit offset files with large records, when last record
variable requires more than 2\*\*32 bytes per record. [NCF-164]

* Fix bug in which passing a NULL path to nc\_open causes failure.
[NCF-173]

* Fixed ncgen bugs in parsing and handling opaque data.

* Fixed ncdump bug, not escaping characters special to CDL in enumeration
labels. [NCF-169]

* Fixed bug reading netCDF int into a C longlong or writing from longlong
to external int on 32-bit platforms with classic format files. The upper
32 bits of the longlong were not cleared on read or used on write.
[NCF-171]

* Resolved some erroneous returns of BADTYPE errors and RANGE errors due
to conflating C memory types with external netCDF types when accessing
classic or 64-bit offset files. [NCF-172]

* Fixed bug with ncdump -t interpreting unit attribute without base time
as a time unit. [NCF-175]

* Changed port for testing remote access test server to increase
reliability of tests.

* Modified ncio mechanism to support multiple ncio packages, so that it is
possible to have e.g. posixio and memio operating at the same time.

* Generation of documentation is disabled by default. Use --enable-doxygen
to generate. [NCF-168]

* Added description of configure flags to installation guide.

* Clarified documentation of arguments to nc**open() and nc**create() and
their default values.

* Fixed doxygen installation guide source file to preserve line breaks in
code and scripts. [NCF-174]

* Cleaned up a bunch of lint issues (unused variables, etc.) and some
similar problems reported by clang static analysis.

* Updated and fixed pkg-config source file netcdf.pc.in to work with
separated netCDF language-specific packages. Also fixed nc-config to
call nf-config, ncxx-config, and ncxx4-config for for backward
compatibility with use of nc-config in current Makefiles. [NCF-165]
[NCF-179]

* 4.2 Released 2012-03-19 (Note: Jira entries include reference to
'[NCF-XX]')

* Completely rebuilt the DAP constraint handling. This primarily affects
users who specify a DAP constraint as part of their URL. [NCF-120]

* Fixed cause of slow nccopy performance on file systems with many records
and large disk block size or many record variables, by accessing data a
record at a time instead of a variable at a time. [NCF-142]

* Performance improvement to DAP code to support fetching partial
variables into the cache; especially important when using nc\_get\_var()
API. A partial variable is one that has ranges attached to the
projection variables (e.g. x[1:10][20:21]) [NCF-157]

* Separate the Fortran and C++ libraries and release the C library and
ncdump/ncgen/nccopy without Fortran or C++. [NCF-24]

* Documentation mostly migrated to Doxygen, from Texinfo. [NCF-26]

* Properly convert vara start/count parameters to DAP [NCF-105][NCF-106]

* Fixed major wasted space from previous default variable chunk sizes
algorithm. [NCF-81]

* Fixed bug in nccopy, in which compression and chunking options were
ignored for netCDF-4 input files. [NCF-79]

* Fixed bug in ncgen in which large variables (more than 2**18 elements)
duplicates the first 2**18 values into subsequent chunks of data
[NCF-154].

* Applied Greg Sjaardema's nccopy bug fix, not compressing output
variables f they were not already using compression on the input file
when shuffle specified. [NCF-162]

* Fixed problem when a URL is provided that contains only a host name.
[NCF-103]

* Fixed behavior of ncgen flags so that -o =\> -lb and, in the absence of
any other markers, make the default be -k1 [NCF-158]

* Created a text INSTALL file for netCDF-4.2 release. [NCF-161]

* Fixed bug in ncgen for vlen arrays as fields of compound types where
datalists for those types was improperly interpreted [NCF-145] (but see
NCF-155).

* Improve use of chunk cache in nccopy utility, making it practical for
rechunking large files. [NCF-85]

* Fixed nccopy bug copying a netCDF-4 file with a chunksize for an
unlimited dimension that is larger than the associated dimension size.
[NCF-139]

* Fixed nccopy bug when rechunking a netCDF-4 file with a chunkspec option
that doesn't explicitly specify all dimensions. [NCF-140]

* Fixed bug in netCDF-4 files with non-coordinate variable with the same
name as a dimension. [NCF-141]

* Incorporated Wei Huang's fix for bug where netCDF-4 sometimes skips over
too many values before adding fill values to an in-memory buffer.
[NCF-143]

* Fixed ncgen bug with netCDF-4 variable-length constants (H/T to Lynton
Appel). [NCF-145]

* Incorporated Peter Cao's performance fixes using HDF5 link iterator for
any group with many variables or types. [NCF-148]

* Incorporated Constantine Khroulev's bug fix for invalid usage of
MPI\_Comm\_f2c in nc\_create\_par. [NCF-135]

* Fixed turning off fill values in HDF5 layers when NOFILL mode is set in
netCDF-4 API (thanks to Karen Schuchardt). [NCF-151]

* Fixed bug with scalar coordinate variables in netCDF-4 files, causing
failure with --enable-extra-tests [NCF-149]

* Cleaned up the definition and use of nulldup. [NCF-92][NCF-93][NCF-94]

* Fixed various '\#include' bugs. [NCF-91][NCF-96][NCF-127]

* v2 API functions modified to properly call the external API instead of
directly calling the netcdf-3 functions. [NCF-100]

* Fixed problem with 64-bit offset format where writing more than 2\*\*31
records resulted in erroneous NC\_EINVALCOORDS error. [NCF-101]

* Restored original functionality of ncgen so that a call with no flags,
only does the syntax check. [NCF-104]

* Corrected misc. test bugs [NCF-107]

* Modified ncdump to properly output various new types (ubyte, ushort,
uint, int64, and uint64). [NCF-111]

* Fixed incorrect link flag for szip in configure.ac [NCF-116]

* ncdump -t now properly parses ISO "T" separator in date-time strings.
[NCF-16]

* ncdump -t "human time" functionality now available for attributes and
bounds variables [NCF-70]

* In ncdump, add -g option to support selection of groups for which data
is displayed. [NCF-11]

* Now supports bluefire platform [NCF-52]

* ncdump now properly displays values of attributes of type NC\_USHORT as
signed shorts [NCF-82]

* Rename some code files so that there are no duplicate filenames.
[NCF-99]

* Demonstration of netCDF-4 Performance Improvement with KNMI Data
[NCF-113]

* Dimension size in classic model netCDF-4 files now allows larger sizes
than allowed for 64-bit offset classic files. [NCF-117]

* ncdump now reports correct error message when "-x" option specifying
NcML output is used on netCDF-4 enhanced model input. [NCF-129]

* Fixed bug causing infinite loop in ncdump -c of netCDF-4 file with
subgroup with variables using inherited dimensions. [NCF-136]

## 4.1.3 2011-06-17

* Replace use of --with-hdf5= and other such configure options that
violate conventions and causes build problems. Set environment variables
CPPFLAGS, LDFLAGS, and LD\_LIBRARY\_PATH instead, before running
configure script. [NCF-20]

* Detect from configure script when szlib is needed [NCF-21]

* Fix bug that can silently zero out portions of a file when writing data
in nofill mode beyond the end of a file, crossing disk-block boundaries
with region to be written while in-memory buffer is in a specific state.
This bug was observed disabling fill mode using Lustre (or other large
blksize file system) and writing data slices in reverse order on disk.
[NCF-22]

* Fix bug that prevents netCDF-4/HDF5 files created with netCDF-4.1.2 from
being read by earlier versions of netCDF or HDF5 versions before 1.8.7.
[NCF-23]

* Fix bug in configure that did not make the search for the xdr library
depend on --enable-dap. [NCF-41]

* Fix ncgen bug that did not use the value of a \_Format attribute in the
input CDL file to determine the kind of output file created, when not
specified by the -k command-line flag. [NCF-42]

* Fix ncgen bug, not properly handling unsigned longlong parsing. [NCF-43]

* Fix DAP client code to suppress variables with names such as "x.y",
which DAP protocol interprets as variable "y" inside container "x". Such
variables will be invisible when accessed through DAP client. [NCF-47]

* Define uint type for unsigned integer, if not otherwise available.
Symptom was compile error involving uint in putget.c. [NCF-49]

* Fix username+password handling in the DAP client code. [NCF-50]

* Add test for handling parallel I/O problem from f77 when user forgets to
turn on one of the two MPI flags. [NCF-60]

* Resolved "make check" problems when ifort compiler. Some "make install"
problems remain when using MPI and shared libraries. [NCF-61]

* Fix problem with f90\_def\_var not always handle deflate setting when
compiler was ifort. [NCF-67]

* Check that either MPIIO or MPIPOSIX flag is set when parallel create or
open is called. Also fix examples that didn't set at least one of these
flags. [NCF-68]

* Improve documentation on handling client-side certificates [NCF-48]

* Document that array arguments, except in varm functions, must point to
contiguous blocks of memory. [NCF-69]

* Get netCDF-4 tests working for DLLs generated with mingw. [NCF-6]

* Make changes necessary for upgrading to HDF5 1.8.7 [NCF-66]

### 4.1.3-rc1 2011-05-06 

* Stop looking for xdr if --disable-dap is used.

* Don't try to run (some) fortran configure tests on machines with no
fortran.

* Allow nccopy to rechunk with chunksizes larger than current dimension
lengths.

* Initial implementation of CDMREMOTE is complete; needs comprehensive
testing.

### 4.1.3-beta1 2011-04-29

* Fixed szlib not linking bug.

* Fixed dreaded "nofill bug", lurking in netCDF classic since at least
1999. Writing more than a disk block's worth of data that crossed disk
block boundaries more than a disk block beyond the end of file in nofill
mode could zero out recently written earlier data that hadn't yet been
flushed to disk.

* Changed setting for H5Pset\_libver\_bounds to ensure that all netCDF-4
files can be read by HDF5 1.8.x.

* Merged libncdap3 and libncdap4 into new libdap2 library. The suffix dap2
now refers to the dap protocol. This is in prep for adding dap4 protocol
support.

* Took out --with-hdf5 and related options due to high cost of maintaining
this non-standard way of finding libraries.

## 4.1.2 2011-03-29

* Changes in build system to support building dlls on cygwin/mingw32.

* Changes to fix portability problems and get things running on all test
platforms.

* Some minor documentation fixes.

* Fixed opendap performance bug for nc\_get\_vars; required adding
nc\_get\_var{s,m} to the dispatch table.

* Now check for libz in configure.ac.

* Fixed some bugs and some performance problems with default chunksizes.

### 4.1.2-beta2 2011-01-11

* Add "-c" option to nccopy to specify chunk sizes used in output in terms
of list of dimension names.

* Rewrite netCDF-4 attribute put code for a large speedup when writing
lots of attributes.

* Fix nc-config --libs when static dependent libraries are not installed
in the same directory as netCDF libraries (thanks to Jeff Whitaker).

* Build shared libraries by default, requiring separate Fortran library.
Static libraries now built only with --disable-shared.

* Refactor of HDF5 file metadata scan for large speedup in opening files,
especially large files.

* Complete rewrite of the handling of character datalist constants. The
heuristics are documented in ncgen.1.

* Eliminate use of NC\_MAX\_DIMS and NC\_MAX\_VARS in ncdump and nccopy,
allocating memory as needed and reducing their memory footprint.

* Add documentation for new nc\_inq\_path() function.

* Use hashing to speedup lookups by name for files with lots of dimensions
and variables (thanks to Greg Sjaardema).

* Add options to nccopy to support uniform compression of variables in
output, shuffling, and fixing unlimited dimensions. Documented in
nccopy.1 man page and User's Guide.

### 4.1.2-beta1 2010-07-09

* Fix "ncdump -c" bug identifying coordinate variables in groups.

* Fix bug in libsrc/posixio.c when providing sizehint larger than default,
which then doesn't get used (thanks to Harald Anlauf).

* Fix netCDF-4 bug caused when doing enddef/redef and then defining
coordinate variable out of order.

* Fixed bug in man4 directory automake file which caused documentation to
be rebuilt after make clean.

* Turned off HDF5 caching when parallel I/O is in use because of its
memory use.

* Refactoring of netCDF code with dispatch layer to decide whether to call
netCDF classic, netCDF-4, or opendap version of a function.

* Refactoring of netCDF-4 memory internals to reduce memory use and end
dependence on NC\_MAX\_DIMS and NC\_MAX\_NAME.

* Modified constraint parser to be more compatible with a java version of
the parser.

* Modified ncgen to utilize iterators internally; should be no user
visible effect.

* Fixed two large-file bugs with using classic format or 64-bit offset
format and accessing multidimensional variables with more than 2\*\*32
values.

## 4.1.1 2010-04-01

* Fixed various build issues.

* Fixed various memory bugs.

* Fixed bug for netCDF-4 files with dimensions and coord vars written in
different orders, with data writes interspersed.

* Added test for HDF5-1.8.4 bug.

* Added new C++ API from Lynton Appel.

## 4.1 2010-01-30

* Much better memory leak checking with valgrind.

* Added per-variable chunk cache control for better performance. Use
nc\_set\_var\_chunk\_cache / nf\_set\_var\_chunk\_cache /
nf90\_set\_var\_chunk\_cache to set the per-variable cache.

* Automatically set per-variable chunk cache when opening a file, or
creating a variable, so that the cache is big enough for more than one
chunk. (Can be overridden by user). Settings may be changed with
configure options --max-default-chunk-size and
--default-chunks-in-cache.

* Better default chunks size. Now chunks are sized to fit inside the
DEFAULT\_CHUNK\_SIZE (settable at configure time with
--with-default-chunk-size= option.)

* Added nccopy utility for converting among netCDF format variants or to
copy data from DAP servers to netCDF files.

* The oc library has been modified to allow the occurrence of alias
definitions in the DAS, but they will be ignored.

* The old ncgen has been moved to ncgen3 and ncgen is now the new ncgen4.

* Modified --enable-remote-tests to be on by default.

* Fixed the nc\_get\_varm code as applied to DAP data sources.

* Added tests for nc-config.

* Many documentation fixes.

* Added capability to use the parallel-netcdf (a.k.a. pnetcdf) library to
perform parallel I/O on classic and 32-bit offset files. Use the
NC\_PNETCDF mode flag to get parallel I/O for non-netcdf-4 files.

* Added libcf library to netCDF distribution. Turn it on with configure
option --with-libcf.

* Added capability to read HDF4 files created with the SD (Scientific
Data) API.

* The DAP support was revised to closely mimic the original libnc-dap
support.

* Significantly revised the data handling mechanism in ncgen4 to more
closely mimic the output from the original ncgen.

* Added prototype NcML output capability to ncgen4. It is specified by the
-lcml flag.

* Added capability to read HDF5 files without dimension scales. This will
allow most existing HDF5 datasets to be read by netCDF-4.

* Fixed bug with endianness of default fill values for integer types when
variables are created with a non-native endiannesss and use the default
fill value.

* Significant refactoring of HDF5 type handling to improve performance and
handle complicated nesting of types in cross-platform cases.

* Added UDUNITS2 to the distribution. Use --with-udunits to build udunits
along with netcdf.

* Made changes suggested by HDF5 team to relax creation-order requirement
(for read-only cases) which allows HDF5 1.6.x files to be retrofitted
with dimension scales, and be readable to netCDF-4.

* Handle duplicate type names within different groups in ncdump. Fix group
path handling in absolute and relative variable names for "-v" option.

* Added nc-config shell script to help users build netCDF programs without
having to figure out all the compiler options they will need.

* Fixed ncdump -s bug with displaying special attributes for classic and
64-bit offset files.

* For writers, nc\_sync() now calls fsync() to flush data to disk sooner.

* The nc\_inq\_type() function now works for primitive types.

## 4.0.1 2009-03-26

* Added optional arguments to F90 API to nf90\_open/create,
nf90\_create\_var, and nf90\_inquire\_variable so that all netCDF-4
settings may be accomplished with optional arguments, instead of
separate function calls.

* Added control of HDF5 chunk cache to allow for user performance tuning.

* Added parallel example program in F90.

* Changed default chunking to better handle very large variables.

* Made contiguous the default for fixed size data sets with no filters.

* Fixed bug in nc\_inq\_ncid; now it returns NC\_ENOGRP if the named group
is not found.

* Fixed man pages for C and F77 so that netCDF-4 builds will result in man
pages that document new netCDF-4 functions.

* Added OPeNDAP support based on a new C-only implementation. This is
enabled using --enable-dap option and requires libcurl. The configure
script will attempt to locate libcurl, but if it fails, then its
location must be specified by the --with-curl option.

### 4.0.1-beta2 2008-12-26

* Changed chunksizes to size\_t from int.

* Fixed fill value problem from F77 API.

* Fixed problems in netcdf-4 files with multi-dimensional coordinate
variables.

* Fixed ncgen to properly handle CDL input that uses Windows line endings
("\r\n"), instead of getting a syntax error.

* Added "-s" option to ncdump to display performance characterisitics of
netCDF-4 files as special virtual attributes, such as \_Chunking,
\_DeflateLevel, \_Format, and \_Endianness.

* Added "-t" option to ncdump to display times in human readable form as
strings. Added code to interpret "calendar" attribute according to CF
conventions, if present, in displaying human-readable times.

* Added experimental version of ncgen4 capable of generating netcdf-4 data
files and C code for creating them. In addition, it supports the special
attributes \_Format, etc.

* 4.0.1-beta1 2008-10-16

* Fixed Fortran 90 int64 problems.

* Rewrote HDF5 read/write code in accordance with performance advice from
Kent.

* Fixed memory leaks in gets/puts of HDF5 data.

* Fixed some broken tests for parallel I/O (i.e. MPI) builds.

* Fixed some cross-compile problems.

* Rewrote code which placed bogus errors on the HDF5 error stack, trying
to open non-existant attributes and variables. Now no HDF5 errors are
seen.

* Removed man subdirectory. Now man4 subdirectory is used for all builds.

* Changed build so that users with access to parallel make can use it.

* Added experimental support for accessing data through OPeNDAP servers
using the DAP protocol (use --enable-opendap to build it).

* Fixed ncdump bugs with array field members of compound type variables.
Fixed ncdump bug of assuming default fill value for data of type
unsigned byte.

## 4.0 2008-05-31

* Introduced the use of HDF5 as a storage layer, which allows use of
groups, user-defined types, multiple unlimited dimensions, compression,
data chunking, parallel I/O, and other features. See the netCDF Users
Guide for more information.

## 3.6.3 2008-05-31

* In ncdump and ncgen, added CDL support for UTF-8 encoding of characters
in names and for escaped special chars in names. Made sure UTF-8 names
are normalized using NFC rules before storing or comparing.

* Handle IEEE NaNs and infinities in a platform-independent way in ncdump
output.

* Added support for ARM representation of doubles, (thanks to Warren
Turkal).

* Fixed bug in C++ API creating 64-bit offset files. (See
http://www.unidata.ucar.edu/software/netcdf/docs/known\_problems.html\#cxx\_64-bit.)

* Fixed bug for variables larger than 4 GB. (See
http://www.unidata.ucar.edu/software/netcdf/docs/known\_problems.html\#large\_vars\_362.)

* Changed the configure.ac to build either 3.6.x or 4.x build from the
same configure.ac.

* Build now checks gfortran version and handles it cleanly, also Portland
Group in Intel fortran, with various configurations.

* A Fortran netcdf.inc file is now created at build time, based on the
setting of --disable-v2.

* Documentation has been fixed in several places.

* Upgraded to automake 1.10, autoconf 2.62, and libtool 2.2.2.

* Includes missing Windows Visual Studio build files.

* Fixed missing include of config.h in a C++ test program.

* Fixed maintainer-clean in man directory.

* Fixed --enable-c-only and make check.

* Fixed behavior when opening a zero-length file.

* Many portability enhancements to build cleanly on various platforms.

* Turned on some old test programs which were not being used in the build.

## 3.6.2 2007-03-05

* Released.

### 3.6.2 beta6 2007-01-20

* Fine tuning of build system to properly handle cygwin, Mingw, and
strange configuration issues.

* Automake 1.10 has a problem with running our tests on MinGW, so I'm
switching back to automake 1.9.6 for this release.

### 3.6.2 beta5 2006-12-30

* Now netCDF configuration uses autoconf 2.61, and automake 1.10. (Thanks
to Ralf Wildenhues for the patches, and all the autotools help in
general!)

* Final major revision of netCDF tutorial before the 3.6.2 release.

* Now netCDF builds under MinGW, producing a windows DLL with the C and
F77 APIs. Use the --enable-shared --enable-dll --disable-cxx
--disable-f90 flags to configure. (C++ and F90 have never been built as
windows DLLs, but might be in a future release if there is user
interest). This has all been documented in the netCDF Porting and
Installation Guide.

* Now extreme numbers (i.e. those close to the limits of their type) can
be turned off in nc\_test/nf\_test, with --disable-extreme-numbers. It
is turned off automatically for Solaris i386 systems.

* Added --enable-c-only option to configure. This causes only the core
netCDF-3 C library to be built. It's the same as --disable-f77
--disable-cxx --disable-v2 --disable-utilities.

* Added --disable-utilities to turn off building and testing of
ncgen/ncdump.

* Fix a long-standing bug in nf90\_get\_att\_text() pointed out by Ryo
Furue, to make sure resulting string is blank-padded on return. This is
fixed in the Fortran-90 interface, but is impractical to fix in the
Fortran-77 interface implemented via cfortran.h.

* Now large file tests are run if --enable-large-file-tests is used in the
configure.

* For Cray users, the ffio module is used if the --enable-ffio option is
passed to configure.

* Unrolled loops in byte-swapping code used on little-endian platforms to
reduce loop overhead. This optimization resulted in a 22% speedup for
some applications accessing floats or ints (e.g. NCO utilities ncap and
ncbo) and a smaller speedup for shorts or doubles.

* Added "-k" option to ncdump and ncgen, for identifying and specifying
the kind of netCDF file, one of "classic", "64-bit-offset", "hdf5", or
"hdf5-nc3". Removed output of kind of netCDF file in CDL comment
produced by ncdump.

* Fixed bug of ncdump seg-faulting if invoked incorrectly with option like
"-c" or "-h" but no file name.

### 3.6.2 beta4 2006-08-15

* Changed F77/F90 man pages from netcdf.3f and netcdf.3f90 to
netcdf\_f77.3 and netcdf\_f90.3. Also fixed broken install of man pages.

* Changed configure script so that "-g -O2" is no longer set as CFLAGS,
CXXFLAGS, and FFLAGS by default if a GNU compiler is being used. Now
nothing is set.

* Changed configure script so that fortran flag is set in config.h.

* Updated Installation and Porting Guide, C++ Interface Guide, F77 and F90
Interface Guides.

* Build with static libraries by default.

* Added configure option --enable-separate-fortran, which causes the
fortran library to be built separately. This is turned on automatically
for shared libraries.

* Improved clarity of error messages.

* Changed configuration to get cygwin DLL and mingw DLL builds working,
for the C library only (i.e. no F77, F90, or C++ APIs).

* Changed type of ncbyte in C++ interface from unsigned char to signed
char, for consistency with C interface. The C++ documentation warned
this change would eventually occur.

* Changed the C++ interface to use only the netCDF-3 C interface instead
of the older netCDF-2 C interface. This has the added benefit that
on-the-fly numeric conversions are now supported using get methods, for
example you can get data of any type as double. When using --disable-v2
flag to configure, the C++ interface can now be built and installed.

### 3.6.2 beta3 2006-05-24

* Changed to use default prefix of /usr/local instead of package-based
prefix of previous releases of netCDF. Use the --prefix argument to the
configure script to override the default.

* Made separate fortran library file, instead of appending fortran library
functions to the C library file, if --enable-separate-fortran is used
during configure (it's turned on automatically if --enable-shared is
used). If uses, the fortran API users must link to *both* the C library
and the new fortran library, like this: -lnetcdff -lnetcdf

* Added netCDF examples in C, C++, F77, F90, and CDL. See the examples
subdirectory.

* Added the NetCDF Tutorial.

* Minor fixes to some of the netCDF documentation.

* Made it possible to build without V2 API using --disable-v2 from
configure.

* Switched to new build system, with automake and libtool. Now shared
libraries are built (as well as static ones) on platforms which support
it. For more information about shared libraries, see
http://www.unidata.ucar.edu/software/netcdf/docs/faq.html\#shared\_intro

* Fixed ncdump crash that happened when no arguments were used.

* Fixed for building with gfortran 4.1.0.

* Important fix for machines whose SIZEOF\_SIZE\_T != SIZEOF\_LONG, such
as NEC-SX, thanks to Stephen Leak.

* Fixed C++ on AIX platform.

* Fixed 64-bit builds on AIX platform.

* Removed bad assertion that could be triggered in rare cases when reading
a small file.

* Added comments in v1hpg.c to clarify purpose of each internal function.

* Make sure filesize is determined in nc\_close() *after* buffers get
flushed.

* Fix long-standing problem resulting in files up to 3 bytes longer than
necessary if there is exactly one record variable of type byte, char, or
short and if the number of values per record for that variable is not
divisible by 4 (or 2 in the case of short). Now the filesize determined
from header info by NC\_calcsize should be correct in all cases.

## 3.6.1 2006-01-31

* Updated installation manual for 3.6.1.

* Changed installation to try to provide correct compiler flags for
compiling in 64-bit mode on Sun, Irix, AIX, and HPUX. (HPUX doesn't work
for me, however). Now run configure with --enable-64bit to get a 64 bit
compile.

* Fixed long-standing bug that would cause small netCDF files to be padded
on the end with zero bytes to 4096 bytes when they were opened and
changed. Now small files should stay small after you change a value.

* Fixed bug in assertions in putget.c that would only be noticed if you
change the manifest constant NC\_MAX\_DIMS in netcdf.h to be different
from NC\_MAX\_VAR\_DIMS.

* Moved test ftest.F from fortran to nf\_test directory, and fixed bug in
ftest.F which caused it to return 0 even if tests failed (no tests were
failing, however). Also renamed some test output files to make things a
little clearer.

* If open for writing, pad with up to 3 extra zero bytes before close to
the correct canonical length, calculated from the header. Previously
files could be short due to not padding when writing in NOFILL mode.

* Doubled arbitrary limits on number of dimensions, variables, attributes,
and length of names.

* Change name of nc\_get\_format() to nc\_inq\_format(). Add analogous
interfaces for nf\_inq\_format(), nf90\_inquire(), and
NcFile::get\_format() to f77, f90, and C++ interfaces. Document new
function in texinfo files. Add minimal test to nc\_test, nf\_test.

### 3.6.1-beta3 2005-02-17

* Added function nc\_get\_format(int ncid, int\* formatp) that returns
either NC\_FORMAT\_CLASSIC or NC\_FORMAT\_64BIT for a CDF1 or CDF2 file,
respectively.

* Added test to nc\_test that detects whether format version was changed
after a file is reopened and define mode is entered.

* Correctly configure for Intel ifort Fortran compiler on Linux.

### 3.6.0-p1 2005-02-18

* Fixed bug that changes CDF2 files to CDF1 files if CDF2 file is reopened
for write access and either an attribute is changed or define mode is
entered.

### 3.6.1-beta2 2005-1-6

* Fixed absoft compile problem. Maybe.

### 3.6.1-beta1 2005-1-3

* Fixed Cygwin C++ problem.

* Fixed large file problem in MS Visual C++.NET environment.

* More information in installation and porting guide.

## 3.6.0 2004-12-16

* Added texinfo source for the documentation.

* Added large file tests to Windows directory in distribution.

* Modified win32 visual studio project files so that m4 is no longer
required to build netcdf under visual studio.

* Modified rules.make to use install instead of cp, fixing install problem
for cygwin users.

* Modified configure/install stuff to support HP-UX.

* Modified configure/install stuff to support G95.

* In the f90 interface, applied Arnaud Desitter's fixes to correct
mismatches between scalar and array arguments, eliminating (legitimate)
complaints by the NAGWare f95 compiler. Also fixed bugs introduced in
3.6.0-beta5 in the mapped array interfaces.

### 3.6.0-beta6 2004-10-05

* Fixed AIX 64-bit/largefile install problems.

* Removed FAQ section from netcdf.texi User's Guide, in deference to
online version that can be kept up to date more easily.

### 3.6.0-beta5 2004-10-04

* Fixed assertion violation on 64-bit platforms when size of last fixed
size variable exceeds 2\^32 - 1.

* Removed another restriction on file size by making record size (derived
from other sizes, not part of the format) an off\_t instead of a
size\_t, when an off\_t is larger than a size\_t. This permits records
to be *much* larger in either classic format or 64-bit-offset format.

* Incorporated patch from Mathis Rosenhauer to improve performance of
Fortran 90 interface for calls to nf90\_put\_var\_TYPE(),
nf90\_get\_var\_TYPE(), nf90\_put\_vara\_TYPE(), and
nf90\_get\_vara\_TYPE() functions by not emulating them with the
corresponding nf90\_put\_varm\_TYPE() and nf90\_get\_varm\_TYPE() calls.

* Added tests for invalid offsets in classic format when defining multiple
large variables.

* Improved installation ease. Have configure script use Large File Support
as a default, if available.

* Add "extra\_test" as a target for testing Large File Support.

### 3.6.0-beta3 2004-08-24

* Upgraded to recent autoconf, changed configure to (hopefully) improve
installation. Also added macros to deal with large file systems.

* Added nf\_set\_default\_format to Fortran interface.

* Added testing to the set\_default\_format functions to nc\_test and
nf\_test.

* Added documentation to the man page for set\_default\_format functions.

* Added two new error return codes to C, f77, and f90 interfaces for
invalid dimension size and for bad variable size. Made test for max
dimension size depend on whether 64-bit offsets used. Fixed bug with
dimension sizes between 2\^31 and 2\^32 (for byte variables).

* Fixed ncdump to properly print dimensions larger than 2\^31.

* Fixed ncgen to properly handle dimensions between 2\^31 and 2\^32.

### 3.6.0-beta2 

* Added -v2 (version 2 format with 64-bit offsets) option to
ncgen, to specify that generated files or generated C/Fortran code
should create 64-bit offset files. Also added -x option to ncgen to
specify use of no-fill mode for fast creation of large files.

* Added function to set default create mode to C interface
(nc\_set\_default\_create).

* Added win32 directory, with NET subdirectory to hold .NET port of
netCDF. To use, open netcdf.sln with Visual Studio, and do a clean and
then a build of either the debug or release builds. Tests will be run as
part of the build process. VC++ with managed extensions is required
(i.e. VC++.NET).

* Added windows installer files to build windows binary installs.

### 3.6.0-beta1 

* By incorporating Greg Sjaardema's patch, added support for
64-bit offset files, which remove many of the restrictions relating to
very large files (i.e. larger than 2 GB.) This introduces a new data
format for the first time since the original netCDF format was
introduced. Files in this new 64-bit offset format can't be read by
earlier versions of netCDF. Users should continue to use the netCDF
classic format unless they need to create very large files.

* The test suite, nc\_test, will now be run twice, once for netCDF classic
format testing, and once for 64-bit offset format testing.

* The implementation of the Fortran-77 interface has been adapted to
version 4.3 of Burkhard Burow's "cfortran.h".

### 3.6.0-alpha 

* Added NEC SX specific optimization for NFILL tunable
parameter in libsrc/putget.c

Added support for the ifc Fortran-90 compiler creating files "netcdf.d"
and "typesizes.d" (instead of ".mod" files).

* Fixed access to iargc and getarg functions from Fortran-90 for NAG f90
compiler, contributed by Harald Anlauf.

## 3.5.1 2004-02-03

* Updated INSTALL.html for Mac OS X (Darwin).

* Made the installation of the netCDF Fortran-90 module file more robust
regarding the name of the file.

* Added support for eight-byte integers in Fortran90 interface.

* Increased advisory limits in C netcdf.h and Fortran netcdf.inc for
maximum number of dimensions, variables, and attributes.

* Changed C++ declarations "friend NcFile" to "friend class NcFile" in
cxx/netcdfcpp.h to conform to standard.

* Added Dan Schmitt's backward compatible extension to the C++ record
interface to work with arbitrary dimension slices.

* Added C++ documentation note that caller is responsible for deleting
pointer returned by Variable::values() method when no longer needed.

* Made C++ interface more standard; the result may not compile on some old
pre-standard C++ compilers.

* Fixed bug in ncgen when parsing values of a multidimensional char
variable that resulted in failure to pad a value with nulls on IRIX.

* Fixed ncdump bug adding extra quote to char variable data when using -fc
or -ff option.

* Fixed so compiling with -DNO\_NETCDF\_2 will work for building without
backward-compatibility netCDF-2 interfaces.

* Eliminated use of ftruncate(), because it fails on FAT32 file systems
under Linux.

* Initialized a pointer in putget.m4 (used to generate putget.c) that was
involved in uninitialized memory references when nc\_test is run under
Purify. Two users had reported seeing crashes resulting from this
problem in their applications.

* Reverted pointer initializations in putget.m4, after testing revealed
these caused a performance problem, resulting in many extra calls to
px\_pgin and px\_pgout when running nc\_test.

* Added checking of size of "dimids" vector in function
nf90\_inquire\_variable(...) and error-returning if it isn't
sufficiently capacious.

* Added variable index to ncvarget() and ncattinq() error messages and
attribute name to ncattinq() error message.

* Tweaked configure script to work with recent C++ compilers.

* Fixed a memory leak in C++ interface, making sure NcVar::cur\_rec[] gets
deleted in NcVar destructor.

* Reimplemented nc\_sync() fix of version 3.5.0 to eliminate performance
penalty when synchronization is unnecessary.

* Changed order of targets in Makefile to build Fortran interface last, as
a workaround for problem with make on AIX platforms.

## 3.5.0 2001-03-23

* Added Fortran 90 interface.

* Changed C macro TIMELEN in file cxx/nctst.cpp to TIMESTRINGLEN to avoid
clash with macro defined on AIX systems in /usr/include/time.h.

* Fixed miswriting of netCDF header when exiting define mode. Because the
header was always written correctly later, this was only a problem if
there was another reader of the netCDF file.

* Fixed explicit synchronizing between netCDF writer and readers via the
nc\_sync(), nf\_sync(), and ncsync() functions.

* Fixed a number of bugs related to attempts to support shrinking the
header in netCDF files when attributes are rewritten or deleted. Also
fixed the problem that nc\_\_endef() did not work as intended in
reserving extra space in the file header, since the extra space would be
compacted again on calling nc\_close().

* Fixed the "redef bug" that occurred when nc\_enddef() or nf\_enddef() is
called after nc\_redef() or nf\_redef(), the file is growing such that
the new beginning of a record variable is in the next "chunk", and the
size of at least one record variable exceeds the chunk size (see
netcdf.3 man page for a description of this tuning parameter and how to
set it). This bug resulted in corruption of some values in other
variables than the one being added.

* The "\*\*" tuning functions for the Fortran interface, nf\*\*create,
nf\*\*open, and nf\*\*enddef, are now documented in the Fortran interface
man pages.

* Add an 'uninstall' target to all the Makefiles. Dave Glowacki
<dglo@SSEC.WISC.EDU> 199810011851.MAA27335

* Added support for multiprocessing on Cray T3E. Hooks added by Glenn, but
the majority of the work was done at NERSC. Also includes changes to
ffio option specification. Patch rollup provided by R. K. Owen
<rkowen@Nersc.GOV>. The following functions are added to the public
interface. nc**create\_mp() nc**open\_mp() nc\_set\_base\_pe()
nc\_inq\_base\_pe()

* Fixed makefile URL for Win32 systems in INSTALL file.

* Made test for UNICOS system in the configure script case independent.

* Ported to the following systems: AIX 4.3 (both /bin/xlc and
/usr/vac/bin/xlc compilers) IRIX 6.5 IRIX64 6.5

* Changed the extension of C++ files from ".cc" to ".cpp". Renamed the C++
interface header file "netcdfcpp.h" instead of "netcdf.hh", changing
"netcdf.hh" to include "netcdfcpp.h" for backward compatibility.

* Treat "FreeBSD" systems the same as "BSD/OS" system w.r.t. Fortran and
"whatis" database.

* Corrected manual pages: corrected spelling of "enddef" (was "endef") and
ensured that the words "index" and "format" will be correctly printed.

* Updated support for Fortran-calling-C interface by updating
"fortran/cfortran.h" from version 3.9 to version 4.1. This new version
supports the Portland Group Fortran compiler (C macro "pgiFortran") and
the Absoft Pro Fortran compiler (C macro "AbsoftProFortran").

* Corrected use of non-integral-constant-expression in specifying size of
temporary arrays in file "libsrc/ncx\_cray.c".

* Added Compaq Alpha Linux workstation example to INSTALL file.

* Ported cfortran.h to Cygnus GNU Win32 C compiler (gcc for Windows).

* Fixed bug in ncdump using same CDL header name when called with multiple
files.

* Added new NULL data type NC\_NAT (Not A Type) to facilitate checking
whether a variable object has had its type defined yet, for example when
working with packed values.

* Fixed use of compile-time macro NO\_NETCDF\_2 so it really doesn't
include old netCDF-2 interfaces, as intended.

* Ported to MacOS X Public Beta (Darwin 1.2/PowerPC).

* Fixed C++ friend declarations to conform to C++ standard.

* Changed INSTALL file to INSTALL.html instead.

## 3.4 1998-03-09

* Fixed ncx\_cray.c to work on all CRAY systems, not just CRAY1. Reworked
USE\_IEG, which was incorrect. Reworked short support. Now USE\_IEG and
otherwise both pass t\_ncx.

* To better support parallel systems, static and malloc'ed scratch areas
which were shared in the library were eliminated. These were made
private and on the stack where possible. To support this, the macros
ALLOC\_ONSTACK and FREE\_ONSTACK are defined in onstack.h.

* The buffered i/o system implementation in posixio.c was reimplemented to
limit the number and size of read() or write() system calls and use
greater reliance on memory to memory copy. This saves a great deal of
wall clock time on slow (NFS) filesystems, especially during
nc\_endef().

* Added performance tuning "underbar underbar" interfaces nc**open(),
nc**create(), and nc\_\_enddef().

* The 'sizehint' contract between the higher layers and the ncio layer is
consistently enforced.

* The C++ interface has been updated so that the deprecated "nclong"
typedef should no longer be required, and casts to nclong no longer
necessary. Just use int or long as appropriate. nclong is still
supported for backwards compatibility.

* The ncdump utility now displays byte values as signed, even on platforms
where the type corresponding to a C char is unsigned (SGI, for example).
Also the ncdump and ncgen utilities have been updated to display and
accept byte attributes as signed numeric values (with a "b" suffix)
instead of using character constants.

* In libsrc/error.c:nc\_strerror(int), explain that NC\_EBADTYPE applies
to "\_FillValue type mismatch".

* Some changes to configure scripts (aclocal.m4), macros.make.in and
ncgen/Makefile to support NEC SUPER-UX 7.2.

* The "usage" messages of ncgen and ncdump include the string returned
from nc\_inq\_libvers().

* Corrected some casts in the library so that all phases of the arithmetic
computing file offsets occurs with "off\_t" type. This allows certain
larger netcdf files to be created and read on systems with larger
(64bit) off\_t.

* In ncgen, multidimensional character variables are now padded to the
length of last dimension, instead of just concatenating them. This
restores an undocumented but convenient feature of ncgen under netCDF-2.
Also, a syntax error is now reliably reported if the netcdf name is
omitted in CDL input.

* Fortran and C code generated by ncgen for netCDF components whose names
contain "-" characters will now compile and run correctly instead of
causing syntax errors.

* The library allows "." characters in names as well as "\_" and "-"
characters. A zero length name "" is explicitly not allowed. The ncgen
utility will now permit "." characters in CDL names as well.

* Memory leaks in the C++ interface NcVar::as\_\*() member functions and
NcFile::add\_var() member function are fixed. The documentation was
fixed where it indicated incorrectly that the library managed value
blocks that the user is actually responsible for deleting.

* he values of the version 2 Fortran error codes have been modified to
make the version 2 Fortran interface more backward compatible at the
source level.

* Added support for systems whose Fortran INTEGER*1 and INTEGER*2 types
are equivalent to the C "long" type but whose C "int" and "long" types
differ. An example of such a system is the NEC SX-4 with the "-ew"
option to the f90 compiler (sheesh, what a system!).

* Fixed Version 2 Fortran compatibility bug: NCVGTG, NCVGGC, NCVPTG, and
NCVPGC didn't work according to the Version 2 documentation if the
innermost mapping value (i.e. IMAP[1]) was zero (indicating that the
netCDF structure of the variable should be used).

## 3.3.1 1997-06-16

* One can now inquire about the number of attributes that a variable has
using the global variable ID.

* The FORTRAN interface should now work on more systems. In particular:

* It should now work with FORTRAN compilers whose "integer*1" datatype is
either a C "signed char", "short", or "int" and whose "integer*2"
datatype is either a C "short" or "int".

* It should now work with FORTRAN compilers that are extremely picky about
source code formatting (e.g. the NAG f90 compiler).

* The dependency on the non-POSIX utility m4(1) for generating the C and
FORTRAN manual pages has been eliminated.

* EXTERNAL statements have been added to the FORTRAN include-file
"netcdf.inc" to eliminate excessive warnings about "unused" variables
(which were actually functions) by some compilers (e.g. SunOS 4.1.3's
f77(1) version 1.x).

* Building the netCDF-3 package no longer requires the existence of the
Standard C macro RAND\_MAX.

* Fixed an ncdump bug resulting in ncdump reporting Attempt to convert
between text & numbers when \_FillValue attribute of a character
variable set to the empty string "".

* Made ncgen tests more stringent and fixed various bugs this uncovered.
These included bugs in handling byte attributes on platforms on which
char is unsigned, initializing scalar character variables in generated C
code under "-c" option, interspersing DATA statements with declaration
statements in generated Fortran code under "-f" option, handling empty
string as a value correctly in generated C and Fortran, and handling
escape characters in strings. The Fortran output under the "-f" option
was also made less obscure and more portable, using automatic conversion
with netCDF-3 interfaces instead of "BYTE", "INTEGER*1", or "INTEGER*2"
declarations.

* Fixed a C++ interface problem that prevented compiling the C++ library
with Digital's cxx compiler.

* Made ncgen "make test" report failure and stop if test resulted in a
failure of generated C or Fortran code.

* The file that you are now reading was created to contain a high-level
description of the evolution of the netCDF-3 package.

## 3.3 1997-05-15

* The production version of the netCDF-3 package was released.

* A comparison of the netCDF-2 and netCDF-3 releases can be found in the
file COMPATIBILITY.

*/
