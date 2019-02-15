FAQ                {#faq}
=======================

[TOC]

General            {#ncFAQGeneral}
=======

What Is netCDF?  {#What-Is-netCDF}
-----------------


NetCDF (network Common Data Form) is a set of interfaces for
array-oriented data access and a [freely](http://www.unidata.ucar.edu/software/netcdf/docs/copyright.html) distributed
collection of data access libraries for C, Fortran, C++, Java, and other
languages. The netCDF libraries support a machine-independent format for
representing scientific data. Together, the interfaces, libraries, and
format support the creation, access, and sharing of scientific data.

NetCDF data is:

-   *Self-Describing*. A netCDF file includes information about the data
    it contains.
-   *Portable*. A netCDF file can be accessed by computers with
    different ways of storing integers, characters, and floating-point
    numbers.
-   *Scalable*. A small subset of a large dataset may be accessed
    efficiently.
-   *Appendable*. Data may be appended to a properly structured netCDF
    file without copying the dataset or redefining its structure.
-   *Sharable*. One writer and multiple readers may simultaneously
    access the same netCDF file.
-   *Archivable*. Access to all earlier forms of netCDF data will be
    supported by current and future versions of the software.

The netCDF software was developed by Glenn Davis, Russ Rew, Ed Hartnett,
John Caron, Dennis Heimbigner, Steve Emmerson, Harvey Davies, and Ward
Fisher at the Unidata Program Center in Boulder, Colorado, with
[contributions](http://www.unidata.ucar.edu/software/netcdf/docs/credits.html) from many other netCDF users.

----------

How do I get the netCDF software package? {#HowdoIgetthenetCDFsoftwarepackage}
-----------------


The latest source distribution, which includes the C libraries and
utility programs, is available from [the NetCDF Downloads
page](/downloads/netcdf/index.jsp). Separate source distributions for
the Java library, Fortran libraries, and C++ libraries are also
available there. Installation instructions are available with the
distribution or [online](http://www.unidata.ucar.edu/software/netcdf/docs/building.html).

Binary distributions of netCDF are available for various platforms from
package management systems such as dpkg, RPM, fink, MacPorts, Homebrew,
OpenCSW, OpenPKG, and the FreeBSD Ports Collection.

----------

How do I convert netCDF data to ASCII or text? {#How-do-I-convert-netCDF-data-to-ASCII-or-text}
-----------------



One way to convert netCDF data to text is to use the **ncdump** tool
that is part of the netCDF software distribution. It is a command line
tool that provides a text representation of a netCDF file's data, just its
metadata, or just the data for specified
variables, depending on what arguments you use. For more information,
see the \ref ncdump_guide documentation.

Another good tool for conversion of netCDF data to text is the ["ncks" program](http://nco.sourceforge.net/nco.html#ncks-netCDF-Kitchen-Sink) that's one of the utility programs in the [NCO (NetCDF Operators)](software.html#NCO) package. Similar capabilities are available using programs from the [CDO (Climate Data Operators)](software.html#CDO) software, commands from [NCL (NCAR Command Language)](software.html#NCL), or various other packages such as [ANAX](http://science.arm.gov/~cflynn/ARM_Tested_Tools/), cdf2asc, and NOESYS, all "third party" netCDF utilities developed and supported by other organizations. You can find more information about these third-party packages on the [Software for Manipulating or Displaying NetCDF Data](software.html) page.

You can also get netCDF data in ASCII from an OPeNDAP server by using a
".ascii" extension with the URL that specifies the data. For details,
see the OPeNDAP page on [Using a Spreadsheet Application with DODS](http://www.opendap.org/useExcel).

Another freely available tool, [netcdf4excel](https://code.google.com/p/netcdf4excel/), has been developed as a netCDF add-in for MS Excel that can facilitate the conversion of netCDF data to and from text form.

Note that **ncdump** and similar tools can print metadata and data values
from netCDF files, but in general they don't understand coordinate
systems specified in the metadata, only variable arrays and their
indices. To interpret georeferencing metadata so you can print the data
within a latitude/longitude bounding box, for example, you need a higher
level tool that interprets conventions for specifying coordinates, such
as the CF conventions. Or you can write a small program using one of the
language APIs that provide netCDF support, for which [examples are available](http://www.unidata.ucar.edu/software/netcdf/examples/programs/).

----------

How do I convert ASCII or text data to netCDF? {#How-do-I-convert-ASCII-or-text-data-to-netCDF}
-----------------


One way to convert data in text form to netCDF is to use the **ncgen**
tool that is part of the netCDF software distribution. Using **ncgen** for
this purpose is a two-step process:

1.  Convert text data to a file in [CDL form](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#CDL-Syntax) using a text
    editor or text manipulation tools
2.  Convert the CDL representation to netCDF using the **ncgen** tool with
    the "-o" or "-b" option

For more information, see the [ncgen documentation](http://www.unidata.ucar.edu/software/netcdf/docs/ncgen-man-1.html).

If you have installed the NCAR Command Language
([NCL](http://www.ncl.ucar.edu/)) software, there are functions
available and described
[here](http://www.ncl.ucar.edu/Applications/list_io.shtml) and
[here](http://www.ncl.ucar.edu/Applications/read_ascii.shtml) for
reading ASCII and tables into NCL and writing the data out to netCDF
files.

With access to [MATLAB](http://www.mathworks.com/), you can create a
schema for the desired netCDF file using
[ncwriteschema](http://www.mathworks.com/help/techdoc/ref/ncwriteschema.html),
read the data using
[textscan](http://www.mathworks.com/help/techdoc/ref/textscan.html), and
write the data to a netCDF file using
[ncwrite](http://www.mathworks.com/help/techdoc/ref/ncwrite.html).

What's new in the latest netCDF release?


[Release notes](http://www.unidata.ucar.edu/software/netcdf/release-notes-latest.html) for the
latest netCDF release are available that describe new features and fixed
bugs since the previous release.

----------

What is the best way to represent [some particular data] using netCDF? {#What-is-the-best-way-to-represent-some-particular-data-using-netCDF}
-----------------

There are many ways to represent the same information in any
general-purpose data model. Choices left up to the user in the case of
netCDF include which information to represent as variables or as
variable attributes; what names to choose for variables, dimensions, and
attributes; what order to use for the dimensions of multidimensional
variables; what variables to include in the same netCDF file; and how to
use variable attributes to capture the structure and meaning of data. We
provide some guidelines in the NetCDF User's Guide (e.g., the section on
[Differences between Attributes and Variables](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Differences-between-Attributes-and-Variables.html#Differences%20between%20Attributes%20and%20Variables))
and in a new web document [Writing NetCDF Files: BestPractices](http://www.unidata.ucar.edu/software/netcdf/BestPractices.html), but we've found that
a little experience helps. Occasionally we have decided it was useful to
change the structure of netCDF files after experience with how the data
is used.

----------

What convention should be used for the names of netCDF files? {#What-convention-should-be-used-for-the-names-of-netCDF-files}
-----------------


NetCDF files should have the file name extension ".nc". The recommended
extension for netCDF files was changed from ".cdf" to ".nc" in 1994 in
order to avoid a clash with the NASA CDF file extension, and now it also
avoids confusion with "Channel Definition Format" files.



----------

Is there a mailing list for netCDF discussions and questions? {#Is-there-a-mailing-list-for-netCDF-discussions-and-questions}
-----------------

The netcdfgroup@unidata.ucar.edu mailing-list is intended for
discussions and announcements about netCDF interfaces, software, and
use. The volume of this list varies widely, from one message per month
to a dozen messages per day (especially after a new release). A message
posted to this mailing-list will be seen by several hundred people, so
it's usually not appropriate for asking simple questions about use. Such
questions should instead be sent to support-netcdf@unidata.ucar.edu.

If you would prefer to get only a single daily digest of the postings to
the netcdfgroup mailing-list, subscribe instead to the digest form of
the mailing-list, containing the same messages but appearing at most
once per day instead of whenever anyone sends a message to the group.

To subscribe or unsubscribe to either of these mailing lists, use one of
these mailing list actions:

* [subscribe: non-digest](mailto:netcdfgroup-join@unidata.ucar.edu) ]
* [subscribe: digest](mailto:netcdfgroup-request@unidata.ucar.edu?subject=subscribe%0A%20%20%20%20%20%20%20%20%20%20digest)
]
* [change subscription options](http://mailman.unidata.ucar.edu/mailman/options/netcdfgroup)
* [view posts](/mailing_lists/archives/netcdfgroup/)
* [search archives](/search.jsp).

----------

Where are some examples of netCDF datasets? {#Where-are-some-examples-of-netCDF-datasets}
-----------------

Here are some [example netCDF files](http://www.unidata.ucar.edu/software/netcdf/examples/files.html).

----------

What is the best way to handle time using netCDF? {#What-is-the-best-way-to-handle-time-using-netCDF}
-----------------


Discussions of conventions for representing time and handling
time-dependent data have been a past topic of discussion on the
netcdfgroup mailing list. When the subject comes up, interesting
discussions often result, so we've archived past discussions on this
subject at
[http://www.unidata.ucar.edu/software/netcdf/time/](http://www.unidata.ucar.edu/software/netcdf/time/).

A summary of Unidata's recommendations is available from
[http://www.unidata.ucar.edu/software/netcdf/time/recs.html](http://www.unidata.ucar.edu/software/netcdf/time/recs.html).
Briefly, we recommend use of the units conventions supported by the
[udunits library](/software/udunits/) for time and other units
attributes.

Other groups have established more specific conventions that include the
representation of time in netCDF files. For more information on such
conventions, see the NetCDF Conventions Page at
[http://www.unidata.ucar.edu/software/netcdf/conventions.html](http://www.unidata.ucar.edu/software/netcdf/conventions.html).

----------

Who else uses netCDF? {#Who-else-uses-netCDF}
-----------------

The netCDF mailing list has over 500 addresses (some of which are
aliases to more addresses) in thirty countries. Several groups have
[adopted netCDF as a standard](http://www.unidata.ucar.edu/software/netcdf/docs/standards.html) for
representing some forms of scientific data.

A somewhat dated description of some of the projects and groups that
have used netCDF is available from
[http://www.unidata.ucar.edu/software/netcdf/usage.html](http://www.unidata.ucar.edu/software/netcdf/usage.html).

----------

What are some references to netCDF? {#What-are-some-references-to-netCDF}
-----------------

A primary reference is the User's Guide:

Rew, R. K., G. P. Davis, S. Emmerson, and H. Davies, **NetCDF User's
Guide for C, An Interface for Data Access, Version 3**, April 1997.

Current online and downloadable documentation is available from the
[documentation directory](http://www.unidata.ucar.edu/software/netcdf/docs/).

Other references include:

Brown, S. A, M. Folk, G. Goucher, and R. Rew, "Software for Portable
Scientific Data Management," Computers in Physics, American Institute of
Physics, Vol. 7, No. 3, May/June 1993, pp. 304-308.

Fulker, D. W., "Unidata Strawman for Storing Earth-Referencing Data,"
Seventh International Conference on Interactive Information and
Processing Systems for Meteorology, Oceanography, and Hydrology, New
Orleans, La., American Meteorology Society, January 1991.

Jenter, H. L. and R. P. Signell, 1992. "[NetCDF: A Freely-Available Software-Solution to Data-Access Problems for Numerical Modelers](http://www.unidata.ucar.edu/software/netcdf/papers/jenter_signell_92.pdf)". Proceedings
of the American Society of Civil Engineers Conference on Estuarine and
Coastal Modeling. Tampa, Florida.

Kuehn, J.A., "Faster Libraries for Creating Network-Portable
Self-Describing Datasets", Proceedings of the 37th Cray User Group
Meeting, (Barcelona, Spain, March 1996), Cray User Group, Inc.

Rew, R. K. and G. P. Davis, "NetCDF: An Interface for Scientific Data
Access," IEEE Computer Graphics and Applications, Vol. 10, No. 4, pp.
76-82, July 1990.

Rew, R. K. and G. P. Davis, "The Unidata netCDF: Software for Scientific
Data Access," Sixth International Conference on Interactive Information
and Processing Systems for Meteorology, Oceanography, and Hydrology,
Anaheim, California, American Meteorology Society, pp. 33-40, February
1990.

Rew, R. K. and G. P. Davis, " [Unidata's netCDF Interface for Data Access: Status and Plans](/netcdf/ams97.html)," Thirteenth International Conference on Interactive Information and Processing Systems for Meteorology, Oceanography, and Hydrology, Anaheim, California, American Meteorology Society, February 1997.

----------

I'm submitting a paper for publication and want to include a citation for use of netCDF software. What reference should I use? {#How-should-I-cite-use-of-netCDF-software}
-----------------

The registered Digital Object Identifier for all versions of netCDF software is `http://doi.org/10.5065/D6H70CW6`.

The following can be used as a citation:

Unidata, (_year_): Network Common Data Form (netCDF) version _nc_version_ [software]. Boulder, CO: UCAR/Unidata. (http://doi.org/10.5065/D6H70CW6)

where _year_ is the year in which the work being described was done and _nc_version_ is the version of netCDF used. For example:

Unidata, (2015): Network Common Data Form (netCDF) version 4.3.3.1 [software]. Boulder, CO: UCAR/Unidata. (http://doi.org/10.5065/D6H70CW6)

----------

Is there a document describing the actual physical format for a Unidata netCDF file? {#Is-there-a-document-describing-the-actual-physical-format-for-a-Unidata-netCDF-file}
-----------------

A short document that specifies the [format of netCDF classic and 64-bit offset files](http://earthdata.nasa.gov/sites/default/files/esdswg/spg/rfc/esds-rfc-011/ESDS-RFC-011v2.00.pdf) has been approved as a standard by the NASA ESDS Software Process Group.

In addition, the NetCDF User's Guide contains an
[appendix](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#File-Format) with the same format specification.

The ["NetCDF File Structure and Performance"](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#Structure) chapter provides a less formal explanation of the format of netCDF data to help clarify the performance implications of different data organizations.

If users only access netCDF data through the documented interfaces, future changes to the format will be transparent.

----------

Installation and Porting                      {#Installation-and-Porting}
================

What does netCDF run on? {#What-does-netCDF-run-on}
-----------------

We test releases on the following operating systems with various compilers:

-   AIX
-   HPUX
-   IRIX, IRIX64
-   Linux
-   MacOS X
-   Solaris
-   Windows (some versions, see below)

The [NetCDF Installation and Porting Guide](@ref getting_and_building_netcdf) explains how to build netCDF from source on various platforms. Often, it's as easy as running

~~~~ {.boldcode}
  ./configure
  make check install
~~~~

----------


How can I use current versions of netCDF-4 with Windows? {#HowcanIusecu}
------------------


See [http://www.unidata.ucar.edu/software/netcdf/docs/winbin.html](http://www.unidata.ucar.edu/software/netcdf/win_netcdf).

How can I use netCDF-4.1 with Windows? {#HowcanIusenetCDF41withWindows}
-----------------


We recently (Summer of 2010) refactored the core building of the netCDF
library. Unfortunately this hopelessly broke the existing port to
Microsoft Visual Studio. Resources permitting, the development of a new
Visual Studio port will be undertaken in the second half of 2010 at
Unidata. Until then, no Visual Studio port of the latest version of the
library is available.

Users are advised that the netCDF build is known to work with Cygwin,
the free POSIX layer for Windows. Building netCDF with Cygwin, and
including the netCDF, HDF5, zlib, and Cygwin DLLs, will allow you to
access the netCDF C library on Windows, even from Visual Studio builds.

We understand that Windows users are most comfortable with a Visual
Studio build, and we intend to provide one.

The Visual Studio port is complicated by the following factors:

-   No configure script support on windows - the Unix build system uses
    a configure script to determine details of the build platform and
    allow the user to specify settings. Windows has no mechanism for
    this other than statically set properties. A Windows-only config.h
    file needs to be created for windows using Cygwin, then included
    with the distribution. Since this contains the version string, it
    must be updated "by hand" before each release.
-   No m4 on windows - the Unix build uses the macro language m4 to
    generate some of the C code in the netCDF library (for example,
    libsrc/putget.c). M4 must be run under Cygwin to generate these
    files, and then they must be statically added to the windows
    distribution. Each new version of netCDF these files should be
    checked for changes. We are restricting new use of m4 for netCDF
    compiles, but that doesn't help with the existing files.
-   No user options on Windows - since Windows does not support a
    configure step, all user options must be pre-set in the Visual
    Studio property lists. As a simplification, many options available
    to Unix users will be unavailable to builders on Windows, such as
    --disable-dap, --disable-netcdf-4, and --disable-shared.
-   Large files (\> 2 GB) have proved to be a problem area in past
    Windows builds.
-   Previous Windows ports have not had to deal with the new OPeNDAP
    client.

Unidata is a community supported organization, and we welcome
collaboration with users who would like to assist with the windows port.
Users should be sure to start with the netCDF daily snapshot, not a
previous release of netCDF.

NOTE: [Paratools](http://www.paratools.com/) has contributed
[instructions for how to build netCDF-4.1.3](http://www.paratools.com/Azure/NetCDF) as a Windows DLL using the MinGW cross compiler.

Nikolay Khabarov has contributed [documentation describing a netCDF-4.1.3 port](http://user.iiasa.ac.at/~khabarov/netcdf-win64-and-win32-mingw/) using MinGW to build native Windows 64-bit and 32-bit DLLs. Current limitations include leaving out support for Fortran and C++ interfaces, NetCDF-4, HDF5, the old version 2 API, and DAP access. The netCDF classic format and 64-bit offset format are fully supported. Links are provided to compiled 32-bit and 64-bit DLLs and static libraries.

A developer on the GMT Wiki has posted [detailed instructions for using CMake](http://gmtrac.soest.hawaii.edu/projects/gmt/wiki/BuildingNetCDF) and MS Visual C++ on Windows to build netCDF-4.1.3, including OPeNDAP support.

Another developer has contributed an unsupported native Windows build of
netCDF-4.1.3 with 32- and 64-bit versions, Fortran bindings, and OPeNDAP
support. The announcement of the availability of that port is
[here](http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2011/msg00363.html).

User Veit Eitner has contributed a port of 4.1.1 to Visual Studio,
including an F90 port to Intel Fortran. Download [source (ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/netcdf-4.1.1-win32-src.zip)](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/netcdf-4.1.1-win32-src.zip) or [binary](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/netcdf-4.1.1-win32-bin.zip) versions. This port was done before the code was refactored in 4.1.2.

How can I use netCDF-4 with Windows? {#How-can-I-use-netCDF-4-with-Windows}
-----------------


Note that we have not ported the F90 or C++ APIs to the Windows
platform, only the C and F77 APIs. User contributions of ports to F90
windows compilers are very welcome (send them to
support-netcdf@unidata.ucar.edu).

On windows, NetCDF consists of a DLL and the ncgen/ncdump executables.
The easiest course is to download one of the pre-built DLLs and
utilities and just install them on your system.

Unlike Unix builds, the Visual Studio build **always** requires HDF5,
zlib, and szlib in all cases. All Windows DLL users must also have the
HDF5, zlib, and szlib DLLs. These are now available from the Unidata FTP
site:

-   [zlib DLLs for 32-bit Windows](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/zlib123-vs2005.zip)
-   [szlib DLLs for 32-bit Windows](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/szip21-vs6-enc.zip)
-   [HDF5 DLLs for 32-bit Windows](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/5-181-win-vs2005.zip)

Two versions of the netCDF DLLs are available, for different Fortran
compilers:

-   [NetCDF for Intel and Portland Group Fortran compilers.](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/win32_vs_PGI_dll_4.0.1.zip)
-   [NetCDF for other Fortran compilers.](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/win32_vs_f2c_dll_4.0.1.zip)

To use netCDF, install the DLLs in /system/win32 and the .h files in a
directory known to your compiler, and define the DLL\_NETCDF
preprocessor macro before including netcdf.h.

The netCDF-4 library can also be built using Visual Studio 2008. Open
the solution file win32/NET/netcdf.sln.

If you install the header files in \\include directory, the netCDF
solution file will work without modifications. Otherwise the properties
of the netcdf project must be changed to include the proper header
directory.

Both the debug and release builds work. The release build links to
different system libraries on Windows, and will not allow debuggers to
step into netCDF library code. This is the build most users will be
interested in. The debug build is probably of interest only to netCDF
library developers.

As of version 4.0.1 (March 2009), the DLL build does not yet include any
testing of the extended netCDF-4 data model. The netCDF4/HDF5 format is
extensively tested in the classic model, but tests for groups,
user-defined types, and other features of the expanded netCDF-4 data
model have not yet been ported to Windows.

The [NetCDF Installation and Porting Guide](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-install/index.html) documents how to
use netCDF with Windows.

Some users have built and released netCDF with Intel Fortran on Windows.
See the [ifort entry in other builds document](http://www.unidata.ucar.edu/software/netcdf/docs/other-builds.html#ifort-361-windows).

Windows is a complicated platform to build on. Some useful explanations
of the oddities of Windows can be found here:

-   Cygwin documentation for [Building and Using DLLs](http://cygwin.com/cygwin-ug-net/dll.html)
-   [OpenLDAP FAQ answer: MinGW Support in Cygwin](http://www.openldap.org/faq/data/cache/301.html), by Jon
    Leichter.
-   [cygwin mailing list explanation of Windows DL requirements.](http://cygwin.com/ml/cygwin/2000-06/msg00688.html)
-   [-mno-cygwin - Building Mingw executables using Cygwin](http://www.delorie.com/howto/cygwin/mno-cygwin-howto.html)

Once you have the netCDF DLL, you may wish to call it from Visual Basic.
The [netCDF VB wrapper](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/netcdf_vb_net_wrapper.zip) will help you do this.

The SDS ([Scientific DataSet](http://research.microsoft.com/en-us/projects/sds/)) library and tools provide .Net developers a way to read, write and share scalars, vectors, and multidimensional grids using CSV, netCDF, and other file formats. It currently uses netCDF version 4.0.1. In addition to .Net libraries, SDS provides a set of utilities and packages: an sds command line utility, a DataSet Viewer application and an add-in for Microsoft Excel 2007 (and later versions).

----------

How do I build and install netCDF for a specific development environment? {#How-do-I-build-and-install-netCDF-for-a-specific-development-environment}
-----------------

You have to build and install the netCDF C library first, before you build and install other language libraries that depend on it, such as Fortran, C++, or Python netCDF libraries. The netCDF Java library is mostly independent of the netCDF C library, unless you need to write netCDF-4 files from Java, in which case you will also need an installed netCDF C library.

For more details, see
[NetCDF Installation and Porting Guide](@ref getting_and_building_netcdf).


----------

How can I tell if I successfully built and installed netCDF? {#How-can-I-tell-if-I-successfully-built-and-installed-netCDF}
-----------------


We make build output from various platforms [available](../builds) for
comparison with your output. In general, you can ignore compiler
warnings if the "make test" step is successful. Lines that begin with
"\*\*\*" in the "make test" output indicate results from tests. The C
and Fortran-77 interfaces are tested extensively, but only rudimentary
tests are currently used for the C++ and Fortran-90 interfaces.

How can I tell what version I'm using? {#How-can-I-tell-what-version-Im-using}
-----------------


If you invoke

~~~~ {.boldcode}
  ncdump --version
~~~~

the last line of the resulting output will identify the version
associated with the **ncdump** utility. You can also call one of the
functions `nc_inq_libvers()`, `nf_inq_libvers()`, or
`nf90_inq_libvers()` from C, Fortran-77, or Fortran-90 programs to get a
version string.

----------

Where does netCDF get installed? {#Where-does-netCDF-get-installed}
-----------------


The netCDF installation directory can be set at the time configure is
run using the --prefix argument. If it is not specified, /usr/local is
used as the default prefix.

For more information see the [NetCDF Installation and Porting Guide](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-install).

Formats, Data Models, and Software Releases {#formatsdatamodelssoftwarereleases}
===========================================

In different contexts, "netCDF" may refer to a data model, a software
implementation with associated application program interfaces (APIs), or
a data format. Confusion may arise in discussions of different versions
of the data models, software, and formats. For example, compatibility
commitments require that new versions of the software support all
previous versions of the format and data model. This section of FAQs is
intended to clarify netCDF versions and help users determine what
version to build and install.

How many netCDF formats are there, and what are the differences among them? {#How-many-netCDF-formats-are-there-and-what-are-the-differences-among-them}
-----------------


There are four netCDF format variants:

-   the classic format
-   the 64-bit offset format
-   the 64-bit data format
-   the netCDF-4 format
-   the netCDF-4 classic model format

(In addition, there are two textual representations for netCDF data,
though these are not usually thought of as formats: CDL and NcML.)

The **classic format** was the only format for netCDF data created
between 1989 and 2004 by the reference software from Unidata. It is
still the default format for new netCDF data files, and the form in
which most netCDF data is stored. This format is also referred as CDF-1 format.

In 2004, the **64-bit offset format** variant was added. Nearly
identical to netCDF classic format, it allows users to create and access
far larger datasets than were possible with the original format. (A
64-bit platform is not required to write or read 64-bit offset netCDF
files.) This format is also referred as CDF-2 format.

In 2008, the **netCDF-4 format** was added to support per-variable
compression, multiple unlimited dimensions, more complex data types, and
better performance, by layering an enhanced netCDF access interface on
top of the HDF5 format.

At the same time, a fourth format variant, **netCDF-4 classic model
format**, was added for users who needed the performance benefits of the
new format (such as compression) without the complexity of a new
programming interface or enhanced data model.

In 2016, the **64-bit data format** variant was added. To support large
variables with more than 4-billion array elements, it replaces most of the
32-bit integers used in the format specification with 64-bit integers. It also
adds support for several new data types including unsigned byte, unsigned
short, unsigned int, signed 64-bit int and unsigned 64-bit int. A 64-bit
platform is required to write or read 64-bit data netCDF files. This format is
also referred as CDF-5 format.

With each additional format variant, the C-based reference software from
Unidata has continued to support access to data stored in previous
formats transparently, and to also support programs written using
previous programming interfaces.

Although strictly speaking, there is no single "netCDF-3 format", that
phrase is sometimes used instead of the more cumbersome but correct
"netCDF classic CDF-1, 64-bit offset CDF-2, or 64-bit data CDF-5 format" to
describe files created by the netCDF-3 (or netCDF-1 or netCDF-2) libraries.
Similarly "netCDF-4 format" is sometimes used informally to mean "either the
general netCDF-4 format or the restricted netCDF-4 classic model format". We
will use these shorter phrases in FAQs below when no confusion is likely.

A more extensive description of the netCDF formats and a formal specification
of the classic and 64-bit formats is available as a [NASA ESDS community
standard](https://earthdata.nasa.gov/sites/default/files/esdswg/spg/rfc/esds-rfc-011/ESDS-RFC-011v2.00.pdf).

The 64-bit data CDF-5 format specification is available in
http://cucis.ece.northwestern.edu/projects/PnetCDF/CDF-5.html.

How can I tell which format a netCDF file uses? {#How-can-I-tell-which-format-a-netCDF-file-uses}
-----------------


The short answer is that under most circumstances, you should not care,
if you use version 4.0 or later of the netCDF library to access data in
the file. But the difference is indicated in the first four bytes of the
file, which are 'C', 'D', 'F', '\\001' for the classic netCDF CDF-1 format;
'C', 'D', 'F', '\\002' for the 64-bit offset CDF-2 format;
'C', 'D', 'F', '\\005' for the 64-bit data CDF-5 format; or '\\211', 'H',
'D', 'F' for an HDF5 file, which could be either a netCDF-4 file or a
netCDF-4 classic model file. (HDF5 files may also begin with a
user-block of 512, 1024, 2048, ... bytes before what is actually an
8-byte signature beginning with the 4 bytes above.)

With netCDF version 4.0 or later, there is an easy way that will
distinguish between netCDF-4 and netCDF-4 classic model files, using the
"-k" option to **ncdump** to determine the kind of file, for example:

~~~~~ {.boldcode}
  ncdump -k foo.nc
  classic
~~~~~


In a program, you can call the function
[nc_inq_format](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c.html#nc_005finq-Family)(or [nf90_inq_format](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html#Compiling-and-Linking-with-the-NetCDF-Library) for the Fortran-90 interface) to determine the format variant of an open netCDF file.

Finally, on a Unix system, one way to display the first four bytes of a
file, say foo.nc, is to run the following command:

~~~~ {.boldcode}
      od -An -c -N4 foo.nc
~~~~

which will output

~~~~ {.boldcode}
      C   D   F 001
~~~~

~~~~ {.boldcode}
      C   D   F 002
~~~~

~~~~ {.boldcode}
      C   D   F 005
~~~~

~~~~ {.boldcode}
    211   H   D   F
~~~~

depending on whether foo.nc is a classic CDF-1, 64-bit offset CDF-2, 64-bit
data CDF-5, or netCDF-4 file, respectively. This method cannot be used to
distinguish between netCDF-4 and netCDF-4 classic model variants, or between a
netCDF-4 file and a different kind of HDF5 file.

----------

How many netCDF data models are there? {#How-many-netCDF-data-models-are-there}
-----------------

There are only two netCDF data models, the [classic
model](/netcdf/workshops/2008/datamodel/NcClassicModel.html) and the [enhanced
model](/netcdf/workshops/2008/netcdf4/Nc4DataModel.html) (also called the
netCDF-4 data model). The classic model is the simpler of the two, and is used
for all data stored in classic CDF-1 format, 64-bit offset CDF-2 format, 64-bit
data CDF-5 format, or netCDF-4 classic model format. The enhanced model
(sometimes also referred to as the netCDF-4 data model) is an extension of the
classic model that adds more powerful forms of data representation and data
types at the expense of some additional complexity. Although data represented
with the classic model can also be represented using the enhanced model,
datasets that use enhanced model features, such as user-defined data types,
cannot be represented with the classic model. Use of the enhanced model
requires storage in the netCDF-4 format.

How many releases of the C-based netCDF software are supported? {#How-many-releases-of-the-C-based-netCDF-software-are-supported}
-----------------


When netCDF version 4.0 was released in June 2008, version 3.6.3 was
released simultaneously, and both releases were supported by Unidata.
Version 3.6.3 supported only the classic and 64-bit offset formats.
Version 4.0 supported both of those format variants by default, and also
the netCDF-4 and netCDF-4 classic model formats, if built using a
previously installed HDF5 library and using the "--enable-netcdf-4"
configure option. Software built from the netCDF-4.0 release without
specifying "--enable-netcdf-4" (the default) was identical to software
built with netCDF-3.6.3. Starting from version 4.4.0, netCDF added support
for CDF-5 format.

Both netCDF-3 and netCDF-4 C libraries are part of a single software
release. The netCDF software may be built to support just the classic
CDF-1 and 64-bit offset CDF-2 formats (the default), 64-bit data CDF-5 format,
or to also support the netCDF-4 and netCDF-4 classic model formats, if the
HDF5-1.8.x library is installed. Unidata no longer supports a separate
netCDF-3-only version of the software, but instead supports both the classic
and enhanced data models and all four format variants in a single source
distribution.

This does not indicate any plan to drop support for netCDF-3 or the
formats associated with netCDF-3. Support for earlier formats and APIs
will continue with all future versions of netCDF software from Unidata.

Should I get netCDF-3 or netCDF-4? {#Should-I-get-netCDF-3-or-netCDF-4}
-----------------

By downloading a current version of netCDF-4, you have the choice to
build either

-   the default netCDF-3 libraries, which support classic CDF-1, 2, and 5
    formats, and the classic data model; or
-   the netCDF-4 libraries, which support netCDF-4 and netCDF-4 classic
    model formats, as well as classic formats, and the
    enhanced data model.

Which version to build depends on how you will use the software.

Installing the simpler netCDF-3 version of the software is recommended
if the following situations apply:

-   all the data you need to access is available in netCDF classic
    formats
-   you are installing netCDF in order to support another software
    package that uses only netCDF-3 features
-   you plan to only write data in a form that netCDF-3 software and
    applications can access
-   you want to delay upgrading to support netCDF-4 until netCDF-4
    formats are more widely used
-   you cannot install the prerequisite HDF5 1.8 software required to
    build and install netCDF-4

Installing the netCDF-4 version of the software is required for any of
the following situations:

-   you need to access netCDF data that makes use of netCDF-4
    compression or chunking
-   you need to access data in all netCDF formats including netCDF-4 or
    netCDF-4 classic model formats
-   you need to write non-record variables larger than 4GiB or record variables with more than 4GiB per record (see ["Have all netCDF size limits been eliminated?"](http://www.unidata.ucar.edu/software/netcdf/docs/faq.html#Large%20File%20Support10))
-   you are installing netCDF to support other software packages that
    require netCDF-4 features
-   you want to write data that takes advantage of compression,
    chunking, or other netCDF-4 features
-   you want to be able to read netCDF-4 classic model data with no
    changes to your current software except relinking with the new
    library
-   you want to benchmark your current applications with the new
    libraries to determine whether the benefits are significant enough
    to justify the upgrade
-   you need to use parallel I/O with netCDF-4 or netCDF-4 classic files

What is the "enhanced data model" of netCDF-4, and how does it differ from the netCDF-3 classic data model? {#whatisenhanceddatamodel}
-------------


The enhanced model (sometimes referred to as the netCDF-4 data model) is
an extension to the [classic model](/netcdf/workshops/2008/datamodel/NcClassicModel.html) that adds more powerful forms of data representation and data types at the expense of some additional complexity. Specifically, it adds six new primitive data types, four kinds of user-defined data types, multiple unlimited
dimensions, and groups to organize data hierarchically and provide
scopes for names. A [picture](/netcdf/workshops/2008/netcdf4/Nc4DataModel.html) of the enhanced data model, with the extensions to the classic model
highlighted in red, is available from the online netCDF workshop.

Although data represented with the classic model can also be represented
using the enhanced model, datasets that use features of the enhanced
model, such as user-defined data types, cannot be represented with the
classic model. Use of added features of the enhanced model requires that
data be stored in the netCDF-4 format.

Why doesn't the new netCDF-4 installation I built seem to support any of the new features? {#Whydoesnt-the-new-netCDF-4-installation-I-built-seem-to-support-any-of-the-new-features}
-----------------


If you built the software from source without access to an HDF5 library,
then only the netCDF-3 library was built and installed. The current
release will build full netCDF-4 support if the HDF5 1.8.x library is
already installed where it can be found by the configure script or
cmake.

Will Unidata continue to support netCDF-3? {#Will-Unidata-continue-to-support-netCDF-3}
-----------------


Yes, Unidata has a commitment to preserving backward compatibility.

Because preserving access to archived data for future generations is
very important:

-   New netCDF software will provide read and write access to *all*
    earlier forms of netCDF data.
-   C and Fortran programs using documented netCDF APIs from previous
    releases will be supported by new netCDF software (after recompiling
    and relinking, if needed).
-   Future releases of netCDF software will continue to support data
    access and API compatibility.

To read compressed data, what changes do I need to make to my netCDF-3 program? {#To-read-compressed-data-what-changes-do-I-need-to-make-to-my-netCDF-3-program}
-----------------


None. No changes to the program source are needed, because the library
handles decompressing data as it is accessed. All you need to do is
relink your netCDF-3 program to the netCDF-4 library to recognize and
handle compressed data.

To write compressed data, what changes do I need to make to my netCDF-3 program? {#To-write-compressed-data-what-changes-do-I-need-to-make-to-my-netCDF-3-program}
-----------------


The **nccopy** utility in versions 4.1.2 and later supports a "-d *level*"
deflate option that copies a netCDF file, compressing all variables
using the specified level of deflation and default chunking parameters,
or you can specify chunking with the "-c" option.

To do this within a program, or if you want different variables to have
different levels of deflation, define compression properties when each
variable is defined. The function to call is
[nc_def_var_deflate](/netcdf-c.html#nc_005fdef_005fvar_005fdeflate)
for C programs, [nf90_def_var_deflate](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html#NF90_005fDEF_005fVAR_005fDEFLATE) for Fortran 90 programs, [NF_DEF_VAR_DEFLATE](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77.html#NF_005fDEF_005fVAR_005fDEFLATE) for Fortran 77. For C++ programs, the experimental cxx4 API may be used,
assuming you have configured with --enable-cxx-4.

Although default variable chunking parameters may be adequate,
compression can sometimes be improved by choosing good chunking
parameters when a variable is first defined. For example, if a 3D field
tends to vary a lot with vertical level, but not so much within a
horizontal slice corresponding to a single level, then defining chunks
to be all or part of a horizontal slice would typically produce better
compression than chunks that included multiple horizontal slices. There
are other factors in choosing chunk sizes, especially matching how the
data will be accessed most frequently. Chunking properties may only be
specified when a variable is first defined. The function to call is
[nc_def_var_chunking](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c.html#nc_005fdef_005fvar_005f)
for C programs,
[nf90_def_var_chunking](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html#NF90_005fDEF_005fVAR_005fCHUNKING)
for Fortran 90 programs, and
[NF_DEF_VAR_CHUNKING](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77.html#NF_005fDEF_005fVAR_005fCHUNKING)
for Fortran 77 programs. For C++ programs, the experimental cxx4 API may
be used, assuming you have configured with --enable-cxx-4.

If I create netCDF-4 classic model files, can they be read by IDL, MATLAB, R, Python and ArcGIS? {#If-I-create-netCDF-4-classic-model-files-can-they-be-read-by-IDL-MATLAB-R-Python-and-ArcGIS}
-----------------


IDL 8.0 ships with support for netCDF-4, including support for OPeNDAP
remote access.

MATLAB 2012a includes netCDF 4 support with OPeNDAP support turned on,
enabling remote access to many kinds of data, as well as use of groups,
compression, and chunking. An example is available demonstrating some of
the new functions. [NCTOOLBOX](http://nctoolbox.github.io/nctoolbox/),
uses netCDF-Java to provide read access to datasets in netCDF-4, GRIB,
GRIB2 and other formats through Unidata's Common Data Model.

R has the [ncdf4 package](http://cirrus.ucsd.edu/~pierce/ncdf/).

Python has the [netcdf4-python package](http://code.google.com/p/netcdf4-python/).

ArcGIS 10.0 can read netcdf4 using the Multidimensional Tools in
ArcToolbox, and in ArcGIS 10.1, the [Multidimensional Supplemental toolbox](http://esriurl.com/MultidimensionSupplementalTools) uses NetCDF4-Python to read OPeNDAP and netCDF4 files, taking advantage of CF conventions if they exist.

What applications are able to deal with *arbitrary* netCDF-4 files? {#What-applications-are-able-to-deal-with-arbitrary-netCDF-4-files}
-----------------

The netCDF utilities **ncdump**, **ncgen**, and **nccopy**, available in
the Unidata C-based netCDF-4 distribution, are able to deal with
arbitrary netCDF-4 files (as well as all other kinds of netCDF files).

How can I convert netCDF-3 files into netCDF-4 files? {#How-can-I-convert-netCDF-3-files-into-netCDF-4-files}
-----------------


Every netCDF-3 file can be read or written by a netCDF version 4 library, so in
that respect netCDF-3 files are already netCDF-4 files and need no conversion.
But if you want to convert a classic format file (CDF-1, 2, or 5) into a
netCDF-4 format or netCDF-4 classic model format file, the easiest way is to
use the **nccopy** utility. For example to convert a classic format file
foo3.nc to a netCDF-4 format file foo4.nc, use:

~~~~~~~~~~~~~~~~~~~~~~~~~ {.boldcode}
  nccopy -k netCDF-4 foo3.nc foo4.nc
~~~~~~~~~~~~~~~~~~~~~~~~~

To convert a classic format file foo3.nc to a netCDF-4 classic
model format file foo4c.nc, you could use:

~~~~~~~~~~~~~~~~~~~~~~~~~~ {.boldcode}
  nccopy -k netCDF-4-classic foo3.nc foo4c.nc
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have installed [NCO](http://www.unidata.ucar.edu/software/netcdf/software.html#NCO), the NCO
utility "ncks" can be used to accomplish the same task, as follows:

~~~~~~~~~~~~~~~~~~~~~~~~ {.boldcode}
  ncks -7 foo3.nc foo4c.nc
~~~~~~~~~~~~~~~~~~~~~~~~

Another method is available for relatively small files, using the **ncdump**
and **ncgen** utilities (built with a netCDF-4 library). Assuming "small3.nc"
is a small classic format netCDF file, you can create an equivalent netCDF-4
file named "small4.nc" as follows:

~~~~ {.boldcode}
  ncdump small3.nc > small.cdl
  ncgen -o small4.nc -k netCDF-4-classic small.cdl
~~~~

Why might someone want to convert netCDF-4 files into netCDF-3 files? {#Why-might-someone-want-to-convert-netCDF-4-files-into-netCDF-3-files}
-----------------


NetCDF-4 classic model files that use compression can be smaller than
the equivalent netCDF-3 files, so downloads are quicker. If they are
then unpacked and converted to the equivalent netCDF-3 files, they can
be accessed by applications that haven't yet upgraded to netCDF-4.

How can I convert netCDF-4 files into netCDF-3 files? {#How-can-I-convert-netCDF-4-files-into-netCDF-3-files}
-----------------


In general, you can't, because netCDF-4 files may have features of the
netCDF enhanced data model, such as groups, compound types,
variable-length types, or multiple unlimited dimensions, for which no
netCDF-3 representation is available. However, if you know that a
netCDF-4 file conforms to the classic model, either because it was
written as a netCDF-4 classic model file, because the program that wrote
it was a netCDF-3 program that was merely relinked to a netCDF-4
library, or because no features of the enhanced model were used in
writing the file, then there are several ways to convert it to a
netCDF-3 file.

You can use the **nccopy** utility. For
example to convert a netCDF-4 classic-model format file foo4c.nc to a
classic format file foo3.nc, use:

~~~~~~~~~~~~~~~~~~~~~~~~~ {.boldcode}
  nccopy -k classic foo4c.nc foo3.nc
~~~~~~~~~~~~~~~~~~~~~~~~~

If you have installed [NCO](http://www.unidata.ucar.edu/software/netcdf/docs/software.html#NCO), the NCO utility "ncks" can be used to accomplish the same task, as follows:

~~~~~~~~~~~~~~~~~~~~~~~~~ {.boldcode}
  ncks -3 foo4c.nc foo3.nc
~~~~~~~~~~~~~~~~~~~~~~~~~

For a relatively small netCDF-4 classic model file, "small4c.nc" for
example, you can also use the **ncdump** and **ncgen** utilities to create an
equivalent netCDF-3 classic format file named "small3.nc" as follows:

~~~~ {.boldcode}
  ncdump small4c.nc > small4.cdl
  ncgen -o small3.nc small4.cdl
~~~~

How can I convert HDF5 files into netCDF-4 files? {#How-can-I-convert-HDF5-files-into-netCDF-4-files}
-----------------


NetCDF-4 intentionally supports a simpler data model than HDF5, which
means there are HDF5 files that cannot be converted to netCDF-4,
including files that make use of features in the following list:

-   Multidimensional data that doesn't use shared dimensions implemented
    using HDF5 "dimension scales". (This restriction was eliminated in
    netCDF 4.1.1, permitting access to HDF5 datasets that don't use
    dimension scales.)
-   Non-hierarchical organizations of Groups, in which a Group may have
    multiple parents or may be both an ancestor and a descendant of
    another Group, creating cycles in the subgroup graph. In the
    netCDF-4 data model, Groups form a tree with no cycles, so each
    Group (except the top-level unnamed Group) has a unique parent.
-   HDF5 "references" which are like pointers to objects and data
    regions within a file. The netCDF-4 data model does not support
    references.
-   Additional primitive types not included in the netCDF-4 data model,
    including H5T\_TIME, H5T\_BITFIELD, and user-defined atomic types.
-   Multiple names for data objects such as variables and groups. The
    netCDF-4 data model requires that each variable and group have a
    single distinguished name.
-   Attributes attached to user-defined types.
-   Stored property lists
-   Object names that begin or end with a space

If you know that an HDF5 file conforms to the netCDF-4 enhanced data
model, either because it was written with netCDF function calls or
because it doesn't make use of HDF5 features in the list above, then it
can be accessed using netCDF-4, and analyzed, visualized, and
manipulated through other applications that can access netCDF-4 files.

The [ncks tool](http://nco.sourceforge.net/nco.html#ncks-netCDF-Kitchen-Sink) of the NCO collection of netCDF utilities can take simple HDF5 data as input and produce a netCDF file as output, so this may work:

~~~~ {.boldcode}
      ncks infile.hdf5 outfile.nc
~~~~

Another tool has been developed to convert HDF5-EOS Aura files to
netCDF-4 files, and it is currently undergoing testing and documentation
before release on the HDF5 web site.

How can I convert netCDF-4 files into HDF5 files? {#How-can-I-convert-netCDF-4-files-into-HDF5-files}
-----------------


Every netCDF-4 or netCDF-4 classic model file can be read or written by
the HDF5 library, version 1.8 or later, so in that respect netCDF-4
files are already HDF5 files and need no conversion.

The way netCDF-4 data objects are represented using HDF5 is described in
detail in the User Manual section ["C.3 The NetCDF-4 Format"](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#NetCDF_002d4-Format).

Why aren't different extensions used for the different formats, for example ".nc3" and ".nc4"? {#why-arent-different-extensions-used}
------------------

The file extension used for netCDF files is purely a convention. The
netCDF libraries don't use the file extension. A user can currently
create a netCDF file with any extension, even one not consistent with
the format of the file.

The **ncgen** utility uses ".nc" as a default extension for output, but this
can be overridden using the "-o" option to specify the name for the
output file. Recent versions of **ncgen** also have a "-k" option to specify
what kind of output file is desired, selecting any of the 4 format
variants, using either a numeric code or a text string. Most other
netCDF client software pays no attention to the file extension, so using
more explicit extensions by convention has no significant drawbacks,
except possibly causing confusion about format differences that may not
be important.

Why is the default of netCDF-4 to continue to create classic files, rather than netCDF-4 files? {#Why-is-the-default-of-netCDF-4-to-continue-to-create-classic-files-rather-than-netCDF-4-files}
-----------------


Until widely used netCDF client software has been adapted or upgraded to
read netCDF-4 data, classic file format is the default for
interoperability with most existing netCDF software.

Can netCDF-4 read arbitrary HDF5 files? {#Can-netCDF-4-read-arbitrary-HDF5-files}
-----------------


No, but it can read many HDF5 files, and more recent versions can access
more HDF5 data. If you want to access HDF5 data through netCDF
interfaces, avoid HDF5 features not included in the netCDF enhanced data
model. For more details see "[How can I convert HDF5 files into netCDF-4 files?](#fv15)", above.

I installed netCDF-3 with --enable-shared, but it looks like the libraries it installed were netCDF-4, with names like libnetcdf.4.dylib. What's going on? {#I-installed-netCDF-3-with---enable-shared-but-it-looks-like-the-libraries-it-installed-were-netCDF-4-with-names-like-libnetcdf4dylib-Whats-going-on}
-----------------


The number used for the shared library name is not related to the netCDF
library version number.

NetCDF-3.6.3 permits UTF-8 encoded Unicode names. Won't this break backward compatibility with previous software releases that didn't allow such names? {#NetCDF-363-permits-UTF-8-encoded-Unicode-names-Wont-this-break-backward-compatibility-with-previous-software-releases-that-didnt-allow-such-names}
-----------------


Earlier versions of the netCDF libraries have always been able to read
data with arbitrary characters in names. The restriction has been on
*creating* files with names that contained "invalid" special characters.
The check for characters used in names occurred when a program tried to
define a new variable, dimension, or attribute, and an error would be
returned if the characters in the names didn't follow the rules.
However, there has never been any such check on reading data, so
arbitrary characters have been permitted in names created through a
different implementation of the netCDF APIs, or through early versions
of netCDF software (before 2.4), which allowed arbitrary names.

In other words, the expansion to handle UTF-8 encoded Unicode characters
and special characters such as \`:' and \` ' still conforms with
Unidata's commitment to backwards compatibility. All old files are still
readable and writable by the new software, and programs that used to
work will still work when recompiled and relinked with the new
libraries. Files using new characters in names will still be readable
and writable by programs that used older versions of the libraries.
However, programs linked to older library versions will not be able to
create new data objects with the new less-restrictive names.

How difficult is it to convert my application to handle arbitrary netCDF-4 files? {#How-difficult-is-it-to-convert-my-application-to-handle-arbitrary-netCDF-4-files}
-----------------


Modifying an application to fully support the new enhanced data model
may be relatively easy or arbitrarily difficult :-), depending on what
your application does and how it is written. Use of recursion is the
easiest way to handle nested groups and nested user-defined types. An
object-oriented architecture is also helpful in dealing with
user-defined types.

We recommend proceeding incrementally, supporting features that are
easier to implement first. For example, handling the six new primitive
types is relatively straightforward. After that, using recursion (or the
group iterator interface used in **nccopy**) to support Groups is not too
difficult. Providing support for user-defined types is more of a
challenge, especially since they can be nested.

The utility program **nccopy**, provided in releases 4.1 and later, shows
how this can be done using the C interface. It copies an input netCDF
file in any of the format variants, handling nested groups, strings, and
any user-defined types, including arbitrarily nested compound types,
variable-length types, and data of any valid netCDF-4 type. It also
demonstrates how to handle variables that are too large to fit in memory
by using an iterator interface. Other generic utility programs can make
use of parts of **nccopy** for more complex operations on netCDF data.

----------

Shared Libraries {#Shared-Libraries}
================

What are shared libraries? {#What-are-shared-libraries}
-----------------


Shared libraries are libraries that can be shared by multiple running
applications at the same time. This **may** improve performance.

For example, if I have a library that provides function foo(), and I
have two applications that call foo(), then with a shared library, only
one copy of the foo() function will be loaded into memory, and both
programs will use it. With static libraries, each application would have
its own copy of the foo() function.

More information on shared libraries can be found at the following
external sites:

-   [The Program-Library HowTo](http://www.tldp.org/HOWTO/Program-Library-HOWTO/index.html),
    by David Wheeler.

-   [Wikipedia Library Entry](http://en.wikipedia.org/wiki/Library_(computer_science))

----------

Can I build netCDF with shared libraries? {#Can-I-build-netCDF-with-shared-libraries}
-----------------


Starting with version 3.6.2, netCDF can build shared libraries on
platforms that support them, but by default netCDF will build static
libraries only. To turn on shared libraries, use the --enable-shared
option to the [netCDF configure script](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-install/Running-the-configure-Script.html).

----------

How do I use netCDF shared libraries? {#How-do-I-use-netCDF-shared-libraries}
-----------------


With netCDF version 3.6.2, shared libraries can be built on platforms
that support them by using the --enable-shared argument to [netCDF configure script](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-install/Running-the-configure-Script.html).

Users of earlier versions of netCDF can build shared libraries by
setting flags correctly during builds.

When you use a static library, the code is copied from the library into
your program when the program is built. The library is only needed at
build time.

With a shared library the code in the library is not copied into your
executable, so the library is needed every time the program is run.

If you write a program that uses the netCDF shared library, the
operating system will have to find it every time your program is run. It
will look in these places:

1.  Directories you specified as shared library locations at **build
    time**. Unfortunately this is done differently with different
    compilers.

2.  Directories specified in the environment variable LD\_RUN\_PATH at
    **build time**.

3.  Directories specified in the OS-specific environment variable for
    this purpose at **run time**. (LD\_LIBRARY\_PATH on Linux and many
    other Unix variants, LOADLIBS on AIX systems, etc.)

4.  A default list of directories that includes /usr/lib (but don't
    install software there!), and may or may not contain places you
    might install netCDF, like /usr/local/lib.

5.  The directories specified in an OS file such as /etc/ld.conf.

By default the netCDF library will be installed in /usr/local/lib. (This
can be overridden with the --prefix option to the [netCDF configure script](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-install/Running-the-configure-Script.html)).

An external site by Arnaud Desitter has a [table of different tools and command line options relating to shared libraries](http://www.fortran-2000.com/ArnaudRecipes/sharedlib.html) on Linux, Solaris, HP-UX, Tru64, AIX, SGI, Win32, MacOS X, VMS (wow!), and OS/390.

For more information about how do to this in Linux users may find it
useful to read this external webpage, some documentation from Caldera, a
Linux distributor: [Specifying directories to be searched by the dynamic linker](http://osr507doc.sco.com/en/tools/ccs_linkedit_dynamic_dirsearch.html).

----------

Large File Support {#Large-File-Support}
================

Was it possible to create netCDF files larger than 2 GiBytes before version 3.6? {#Was-it-possible-to-create-netCDF-files-larger-than-2-GiBytes-before-version-36}
-----------------


Yes, but there are significant restrictions on the structure of large
netCDF files that result from the 32-bit relative offsets that are part
of the classic netCDF format. For details, see [NetCDF Classic Format Limitations](http://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf/NetCDF-Classic-Format-Limitations.html#NetCDF-Classic-Format-Limitations)
in the User's Guide.

----------

What is Large File Support? {#What-is-Large-File-Support}
-----------------


Large File Support (LFS) refers to operating system and C library
facilities to support files larger than 2 GiB. On a few 32-bit platforms
the default size of a file offset is still a 4-byte signed integer,
which limits the maximum size of a file to 2 GiB. Using LFS interfaces
and the 64-bit file offset type, the maximum size of a file may be as
large as 2^63^ bytes, or 8 EiB. For some current platforms, large file
macros or appropriate compiler flags have to be set to build a library
with support for large files. This is handled automatically in netCDF
3.6 and later versions.

More information about Large File Support is available from [Adding Large File Support to the Single UNIX Specification](http://www.unix.org/version2/whatsnew/lfs.html).

----------

What does Large File Support have to do with netCDF? {#What-does-Large-File-Support-have-to-do-with-netCDF}
-----------------


When the netCDF format was created in 1988, 4-byte fields were reserved
for file offsets, specifying where the data for each variable started
relative to the beginning of the file or the start of a record boundary.

This first netCDF format variant, the only format supported in versions
3.5.1 and earlier, is referred to as the netCDF *classic* format. The
32-bit file offset in the classic format limits the total sizes of all
but the last non-record variables in a file to less than 2 GiB, with a
similar limitation for the data within each record for record variables.
For more information see [Classic Format Limitations](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/NetCDF-Classic-Format-Limitations.html#NetCDF-Classic-Format-Limitations).

The netCDF classic format is also identified as *version 1* or *CDF1* in
reference to the format label at the start of a file.

With netCDF version 3.6 and later, a second variant of netCDF format is
supported in addition to the classic format. The new variant is referred
to as the *64-bit offset* format, *version 2*, or *CDF-2*. The primary
difference from the classic format is the use of 64-bit file offsets
instead of 32-bit offsets, but it also supports larger variable and
record sizes.

Starting from version 4.4.0, netCDF added support for CDF-5 format, which
allows multiple large variables with more than 4-billion array elements defined
in the file. This format is only supported on 64-bit machine platforms.

----------

Do I have to know which netCDF file format variant is used in order to access or modify a netCDF file? {#Do-I-have-to-know-which-netCDF-file-format-variant-is-used-in-order-to-access-or-modify-a-netCDF-file}
-----------------


No, version 3.6 and later versions of the netCDF C/Fortran library
detect which variant of the format is used for each file when it is
opened for reading or writing, so it is not necessary to know which
variant of the format is used. The version of the format will be
preserved by the library on writing. If you want to modify a classic
format file to use the CDF-2 or CDF-5 format so you can make it much
larger, you will have to create a new file and copy the data to it. The
**nccopy** utility available in version 4.1 can copy a classic file to a
CDF-2 or CDF-5 file.

----------

Will future versions of the netCDF library continue to support accessing files in the classic format? {#Will-future-versions-of-the-netCDF-library-continue-to-support-accessing-files-in-the-classic-format}
-----------------


Yes, the 3.6 library and all planned future versions of the library will
continue to support reading and writing files using the classic CDF-1 (32-bit
offset), 64-bit offset CDF-2, and 64-bit data CDF-5 format. There is no need to
convert existing archives from the classic to the 64-bit offset format.
Even netCDF-4, which introduces a third variant of the netCDF format
based on HDF5, continues to support accessing classic CDF-1, 2, and 5 format
files. NetCDF-4 HDF5 files have even fewer restrictions on size than CDF-1 and
CDF-2 files.

----------

Should I start using the new 64-bit offset format for all my netCDF files? {#Should-I-start-using-the-new-64-bit-offset-format-for-all-my-netCDF-files}
-----------------


No, we discourage users from making use of the 64-bit offset format
unless they need it for large files. It may be some time until
third-party software that uses the netCDF library is upgraded to 3.6 or
later versions that support the large file facilities, so we advise
continuing to use the classic netCDF format for data that doesn't
require file offsets larger than 32 bits. The library makes this
recommendation easy to follow, since the default for file creation is
the classic format.

----------

How can I tell if a netCDF file uses the classic format (CDF-1), 64-bit offset format (CDF-2) or 64-bit data format (CDF-5)? {#How-can-I-tell-if-a-netCDF-file-uses-the-classic-format-or-64-bit-offset-format}
-----------------


The short answer is that under most circumstances, you should not care,
if you use version 3.6.0 or later of the netCDF library. But the
difference is indicated in the first four bytes of the file, which are
'C', 'D', 'F', '\\001' for the classic CDF-1 format, 'C', 'D', 'F',
'\\002' for the 64-bit offset CDF-2 format, and 'C', 'D', 'F', '\\005' for the
64-bit data CDF-5 format. On a Unix system, one way to display the first four
bytes of a file, say foo.nc, is to run the following command:

~~~~ {.boldcode}
  od -An -c -N4 foo.nc
~~~~

which will output

~~~~ {.boldcode}
  C   D   F 001
~~~~

or

~~~~ {.boldcode}
  C   D   F 002
~~~~

or

~~~~ {.boldcode}
  C   D   F 005
~~~~

depending on whether foo.nc is a CDF-1, CDF-2, or CDF-5 netCDF file,
respectively.

With netCDF version 3.6.2 or later, there is an easier way, using the
"-k" option to **ncdump** to determine the kind of file, for example:

~~~~ {.boldcode}
  ncdump -k foo.nc
  classic
~~~~

----------

What happens if I create a 64-bit offset format netCDF file and try to open it with an older netCDF application that hasn't been linked with netCDF 3.6? {#What-happens-if-I-create-a-64-bit-offset-format-netCDF-file-and-try-to-open-it-with-an-older-netCDF-application-that-hasnt-been-linked-with-netCDF-36}
-----------------


The application will indicate an error trying to open the file and
present an error message equivalent to "not a netCDF file". This is why
it's a good idea not to create 64-bit offset netCDF files until you
actually need them.

----------

Can I create 64-bit offset files on 32-bit platforms? {#Can-I-create-64-bit-offset-files-on-32-bit-platforms}
-----------------


Yes, by specifying the appropriate file creation flag you can create
64-bit offset netCDF files the same way on 32-bit platforms as on 64-bit
platforms. You do not need to compile the C/Fortran libraries as 64-bit
to support access to 64-bit offset netCDF files.

----------

How do I create a 64-bit offset netCDF file from C, Fortran-77, Fortran-90, or C++? {#How-do-I-create-a-64-bit-offset-netCDF-file-from-C-Fortran-77-Fortran-90-or-Cpp}
-----------------


With netCDF version 3.6.0 or later, use the NC\_64BIT\_OFFSET flag when
you call nc\_create(), as in:

~~~~ {.boldcode}
  err = nc_create("foo.nc",
                  NC_NOCLOBBER | NC_64BIT_OFFSET,
                  &ncid);
~~~~

In Fortran-77, use the NF\_64BIT\_OFFSET flag when you call
nf\_create(), as in:

~~~~ {.boldcode}
  iret = nf_create('foo.nc',
                   IOR(NF_NOCLOBBER,NF_64BIT_OFFSET),
                   ncid)
~~~~

In Fortran-90, use the NF90\_64BIT\_OFFSET flag when you call
nf90\_create(), as in:

~~~~ {.boldcode}
  iret = nf90_create(path="foo.nc",
                     cmode=or(nf90_noclobber,nf90_64bit_offset),
                     ncid=ncFileID)
~~~~

In C++, use the Offset64Bits enum in the NcFile constructor, as in:

~~~~ {.boldcode}
  NcFile nc("foo.nc",
            FileMode=NcFile::New,
            FileFormat=NcFile::Offset64Bits);
~~~~

In Java, use the setLargeFile() method of the NetcdfFileWritable class.

----------

How do I create a 64-bit offset netCDF file using the ncgen utility? {#How-do-I-create-a-64-bit-offset-netCDF-file-using-the-ncgen-utility}
-----------------


A command-line option, '-k', specifies the kind of file format
variant. By default or if '-k classic' is specified, the generated
file will be in netCDF classic format.  If '-k 64-bit-offset' is
specified, the generated file will use the 64-bit offset format.

----------

Have all netCDF size limits been eliminated? {#Have-all-netCDF-size-limits-been-eliminated}
-----------------


The netCDF-4 HDF5-based format has no practical limits on the size of a
variable.

However, for the classic and 64-bit offset formats there are still
limits on sizes of netCDF objects. Each fixed-size variable (except the
last, when there are no record variables) and the data for one record's
worth of a single record variable (except the last) are limited in size
to a little less that 4 GiB, which is twice the size limit in versions
earlier than netCDF 3.6.

The maximum number of records remains 2^32^-1.

----------

Why are variables still limited in size? {#Why-are-variables-still-limited-in-size}
-----------------


While most platforms support a 64-bit file offset, many platforms only
support a 32-bit size for allocated memory blocks, array sizes, and
memory pointers. In C developer's jargon, these platforms have a 64-bit
`off_t` type for file offsets, but a 32-bit `size_t` type for size of
arrays. Changing netCDF to assume a 64-bit `size_t` would restrict
netCDF's use to 64-bit platforms.

----------

How can I write variables larger than 4 GiB? {#How-can-I-write-variables-larger-than-4-GiB}
-----------------


You can overcome the 4 GiB size barrier by using the netCDF-4 HDF5
format for your data. The only change required to the program that
writes the data is an extra flag to the file creation call, followed by
recompiling and relinking to the netCDF-4 library. Programs that access
the data would also need to be recompiled and relinked to the netCDF-4
library.

For classic and 64-bit offset netCDF formats, if you change the first
dimension of a variable from a fixed size to an unlimited size instead,
the variable can be much larger. Even though record variables are
restricted to 4 Gib per record, there may be 4 billion records. NetCDF
classic or 64-bit offset files can only have one unlimited dimension, so
this won't work if you are already using a record dimension for other
purposes.

It is also possible to overcome the 4 GiB variable restriction for a
single fixed size variable, when there are no record variables, by
making it the last variable, as explained in the example in [NetCDF Classic Format Limitations](http://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf/NetCDF-Classic-Format-Limitations.html#NetCDF-Classic-Format-Limitations).

----------

Why do I get an error message when I try to create a file larger than 2 GiB with the new library? {#Why-do-I-get-an-error-message-when-I-try-to-create-a-file-larger-than-2-GiB-with-the-new-library}
-----------------


There are several possible reasons why creating a large file can fail
that are not related to the netCDF library:

-   User quotas may prevent you from creating large files. On a Unix
    system, you can use the "ulimit" command to report limitations such
    as the file-size writing limit.

-   There is insufficient disk space for the file you are trying to
    write.

-   The file system in which you are writing may not be configured to
    allow large files. On a Unix system, you can test this with a
    commands such as

    ~~~~ {.boldcode}
      dd if=/dev/zero bs=1000000 count=3000 of=./largefile
      ls -l largefile
      rm largefile
    ~~~~

    which should write a 3 GByte file named "largefile" in the current
    directory, verify its size, and remove it.

If you get the netCDF library error "One or more variable sizes violate
format constraints", you are trying to define a variable larger than
permitted for the file format variant. This error typically occurs when
leaving "define mode" rather than when defining a variable. The error
status cannot be returned when a variable is first defined, because the
last fixed-size variable defined is permitted to be larger than other
fixed-size variables (when there are no record variables).

Similarly, the last record variable may be larger than other record
variables. This means that subsequently adding a small variable to an
existing file may be invalid, because it makes what was previously the
last variable now in violation of the format size constraints. For
details on the format size constraints, see the Users Guide sections
[NetCDF Classic Format Limitations](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#Classic-Limitations) and [NetCDF 64-bit Offset Format Limitations](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#64-bit-Offset-Limitations).

If you get the netCDF library error "Invalid dimension size" for a
non-negative size, you are exceeding the size limit of netCDF
dimensions, which must be less than 2,147,483,644 for classic files with
no large file support and otherwise less than 4,294,967,292.

----------

Do I need to use special compiler flags to compile and link my applications that use netCDF with Large File Support? {#Do-I-need-to-use-special-compiler-flags-to-compile-and-link-my-applications-that-use-netCDF-with-Large-File-Support}
-----------------


No, except that 32-bit applications should link with a 32-bit version of
the library and 64-bit applications should link with a 64-bit library,
similarly to use of other libraries that can support either a 32-bit or
64-bit model of computation. But note that a 32-bit version of the
netCDF library fully supports writing and reading 64-bit offset netCDF
files.

----------

Is it possible to create a "classic" format netCDF file with netCDF version 3.6.0 that cannot be accessed by applications compiled and linked against earlier versions of the library? {#isitpossibleclassic360}
----------------

No, classic files created with the new library should be compatible with
all older applications, both for reading and writing, with one minor
exception. The exception is due to a correction of a netCDF bug that
prevented creating records larger than 4 GiB in classic netCDF files
with software linked against versions 3.5.1 and earlier. This limitation
in total record size was not a limitation of the classic format, but an
unnecessary restriction due to the use of too small a type in an
internal data structure in the library.

If you want to always make sure your classic netCDF files are readable
by older applications, make sure you don't exceed 4 GiBytes for the
total size of a record's worth of data. (All records are the same size,
computed by adding the size for a record's worth of each record
variable, with suitable padding to make sure each record begins on a
byte boundary divisible by 4.)

----------

NetCDF and Other Software {#NetCDF-and-Other-Software}
================

What other software is available for accessing, displaying, and manipulating netCDF data? {#What-other-software-is-available-for-accessing-displaying-and-manipulating-netCDF-data}
-----------------


Utilities available in the current netCDF distribution from Unidata are
**ncdump**, for converting netCDF files to an ASCII human-readable form,
and **ncgen** for converting from the ASCII human-readable form back to
a binary netCDF file or a C or FORTRAN program for generating the netCDF
file. [Software for Manipulating or Displaying NetCDF Data](software.html) provides a list of other software useful for access, visualization, and analysis of netCDF data and data represented in other forms. Another useful [guide to netCDF utilities](http://nomads.gfdl.noaa.gov/sandbox/products/vis/data/netcdf/GFDL_VG_NetCDF_Utils.html) is available from NOAA's Geophysical Fluid Dynamics Laboratory.

----------

What other data access interfaces and formats are available for scientific data? {#What-other-data-access-interfaces-and-formats-are-available-for-scientific-data}
-----------------


The [Scientific Data Format Information FAQ](http://www.cv.nrao.edu/fits/traffic/scidataformats/faq.html) provides a somewhat dated description of other access interfaces and formats for scientific data, including [CDF](http://nssdc.gsfc.nasa.gov/cdf/cdf_home.html) and [HDF](http://hdf.ncsa.uiuc.edu/). A brief comparison of CDF, netCDF, and HDF is available in the [CDF FAQ](http://nssdc.gsfc.nasa.gov/cdf/html/FAQ.html). Another comparison is in Jan Heijmans' [An Introduction to Distributed Visualization](http://www.xi-advies.nl/downloads/AnIntroductionToDistributedVisualization.pdf). John May's book [*Parallel I/O for High Performance Computing*](http://www.llnl.gov/CASC/news/johnmay/John_May_book.html) includes a chapter on Scientific Data Libraries that describes netCDF and HDF5, with example source code for reading and writing files using both interfaces.

----------

What is the connection between netCDF and CDF? {#What-is-the-connection-between-netCDF-and-CDF}
-----------------


[CDF](http://cdf.gsfc.nasa.gov/) was developed at the NASA Space Science
Data Center at Goddard, and is freely available. It was originally a VMS
FORTRAN interface for scientific data access. Unidata reimplemented the
library from scratch to use [XDR](http://www.faqs.org/rfcs/rfc1832.html)
for a machine-independent representation, designed the
[CDL](http://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf/CDL-Syntax.htm) (network Common Data form Language) text
representation for netCDF data, and added aggregate data access, a
single-file implementation, named dimensions, and variable-specific
attributes.

NetCDF and CDF have evolved independently. CDF now supports many of the
same features as netCDF (aggregate data access, XDR representation,
single-file representation, variable-specific attributes), but some
differences remain (netCDF doesn't support native-mode representation,
CDF doesn't support named dimensions). There is no compatibility between
data in CDF and netCDF form, but NASA makes available [some
translators](http://cdf.gsfc.nasa.gov/html/dtws.html) between various
scientific data formats. For a more detailed description of differences
between CDF and netCDF, see the [CDF FAQ](http://cdf.gsfc.nasa.gov/html/FAQ.html).

----------

What is the connection between netCDF and HDF? {#What-is-the-connection-between-netCDF-and-HDF}
-----------------


The National Center for Supercomputing Applications (NCSA) originally
developed [HDF4](http://hdf.ncsa.uiuc.edu/) and made it freely
available. HDF4 is an extensible data format for self-describing files
that was developed independently of netCDF. HDF4 supports both C and
Fortran interfaces, and it has been successfully ported to a wide
variety of machine architectures and operating systems. HDF4 emphasizes
a single common format for data, on which many interfaces can be built.

NCSA implemented software that provided a netCDF-2 interface to HDF4.
With this software, it was possible to use the netCDF calling interface
to place data into an HDF4 file.

HDF5, developed and supported by The HDF Group, Inc., a non-profit
spin-off from the NCSA group, provides a richer data model, with
emphasis on efficiency of access, parallel I/O, and support for
high-performance computing. The netCDF-4 project has implemented an
enhanced netCDF interface on the HDF5 storage layer to preserve the
desirable common characteristics of netCDF and HDF5 while taking
advantage of their separate strengths: the widespread use and simplicity
of netCDF and the generality and performance of HDF5.

----------

Has anyone implemented client-server access for netCDF data? {#Has-anyone-implemented-client-server-access-for-netCDF-data}
-----------------


Yes, as part of the [OPeNDAP](http://www.opendap.org/) framework,
developers have implemented a client-server system for access to remote
data that supports use of the netCDF interface for clients. A reference
version of the software is available from the [OPeNDAP download site](http://www.opendap.org/download/index.html/). After linking your netCDF application with the OPeNDAP netCDF library, you can use URL's to access data from other sites running an OPeNDAP server. This supports accessing small subsets of large datasets remotely through the netCDF interfaces, without copying the datasets.

The 4.1 release of netCDF will include OPeNDAP client support; an
experimental version is available now in the snapshot distributions.

Other clients and servers support access through a netCDF interface to
netCDF and other kinds of data, including clients written using the
[netCDF-Java library](http://www.unidata.ucar.edu/software/netcdf-java/) and servers that use the
[THREDDS Data Server](/software/thredds/current/tds/TDS.html).

The [GrADS Data Server](http://grads.iges.org/grads/gds/) provides
subsetting and analysis services across the Internet for any
GrADS-readable dataset, including suitable netCDF datasets. The latest
version of the [PMEL Live Access Server](http://ferret.pmel.noaa.gov/LAS) uses THREDDS Data Server technology to provide flexible access to geo-referenced scientific data, including netCDF data.

----------

How do I convert between GRIB and netCDF? {#How-do-I-convert-between-GRIB-and-netCDF}
-----------------


Several programs and packages have been developed that convert between
[GRIB](http://www.wmo.ch/web/www/DPS/grib-2.html) and netCDF data:
[ncl_convert2nc](http://www.ncl.ucar.edu/Applications/grib2netCDF.shtml),
[degrib](http://www.nws.noaa.gov/mdl/NDFD_GRIB2Decoder/),
[CDAT](software.html#CDAT), [CDO](software.html#CDO),
[GDAL](http://www.gdal.org/), [GrADS](software.html#GrADS), and
[wgrib2](http://www.cpc.noaa.gov/products/wesley/wgrib2/).

The Unidata [netCDF Java Library](http://www.unidata.ucar.edu/software/netcdf-java/index.html) can
read GRIB1 and GRIB2 data (and many other data formats) through a netCDF
interface. As a command-line example, you could convert *fileIn.grib* to
*fileOut.nc* as follows:

~~~~ {.boldcode}
  java -Xmx1g -classpath netcdfAll-4.3.jar ucar.nc2.dataset.NetcdfDataset \
    -in fileIn.grib -out fileOut.nc [-isLargeFile] [-netcdf4]
~~~~

For more details on using netCDF Java, see the CDM man pages for
[nccopy](http://www.unidata.ucar.edu/software/netcdf-java/reference/manPages.html#nccopy).

----------

Problems and Bugs
-----------------

Can I recover data from a netCDF file that was not closed properly? {#Can-I-recover-data-from-a-netCDF-file-that-was-not-closed-properly}
-----------------


_I have some netcdf files which have data in them and were apparently
not properly closed. When I examine them using **ncdump** they report zero
data points, although the size is a few megabytes. Is there a way of
recovering them?_

If the files are in classic format or 64-bit offset format (if they were
created by netCDF version 3.6.3 or earlier, for example), then you can
use an editor that allows you to change binary files, such as emacs, to
correct the four-byte number of records field in the file. This is a
bigendian 4 byte integer that begins at the 4th byte in the file.

This is what the first eight bytes would look like for classic format if
you had zero records, where printable characters are specified as
US-ASCII characters within single-quotes and non-printable bytes are
denoted using a hexadecimal number with the notation '\\xDD', where each
D is a hexadecimal digit:

~~~~ {.boldcode}
  'C' 'D' 'F' \x01 \x00 \x00 \x00 \x00
~~~~

or

~~~~ {.boldcode}
  'C' 'D' 'F' \x02 \x00 \x00 \x00 \x00
~~~~

for 64-bit-offset format.

And this is what the first eight bytes should look like for classic
format if you had 500 records (500 is 01F4 in hexadecimal)

~~~~ {.boldcode}
  'C' 'D' 'F' \x01 \x00 \x01 \x0f \x04
~~~~

or

~~~~ {.boldcode}
  'C' 'D' 'F' \x02 \x00 \x01 \x0f \x04
~~~~

for 64-bit-offset format.

So if you can compute how many records should be in the file, you can
edit the second four bytes to fix this. You can find out how many
records should be in the file from the size of the file and from the
variable types and their shapes. See the [description of the netCDF format](http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#File-Format)
for classic and 64-bit offset files for how to figure out how large the
file should be for fixed sized variables of particular shapes and for a
specified number of record variables of particular shapes.

Note that if you neglected to call the appropriate netCDF close function
on a file, data in the last record written but not flushed to the disk
may also be lost, but correcting the record count should allow recovery
of the other records.

----------

Is there a list of reported problems and workarounds? {#Is-there-a-list-of-reported-problems-and-workarounds}
-----------------


Yes, the document [Known problems with the netCDF Distribution](known_problems.html) describes reported problems and workarounds in the latest version and some earlier releases.

----------

How do I make a bug report? {#How-do-I-make-a-bug-report}
-----------------


If you find a bug, send a description to
support-netcdf@unidata.ucar.edu. This is also the address to use for
questions or discussions about netCDF that are not appropriate for the
entire netcdfgroup mailing list.

----------

How do I search through past problem reports? {#How-do-I-search-through-past-problem-reports}
-----------------


A search link is available at the bottom of the [netCDF homepage](http://www.unidata.ucar.edu/software/netcdf/), providing a full-text search of the
support questions and answers about netCDF provided by Unidata support
staff.

----------

Programming with NetCDF {#Programming-with-NetCDF}
================

Which programming languages have netCDF interfaces? {#Which-programming-languages-have-netCDF-interfaces}
-----------------

The netCDF distribution comes with interfaces for C, Fortran77,
Fortran90, and C++. Other languages for which interfaces are available
separately include:

-   [Ada](http://freshmeat.net/projects/adanetcdf/)
-   [IDL](software.html#IDL)
-   [Java](software.html#Java%20interface)
-   [MATLAB](software.html#MATLAB)
-   [Perl](software.html#Perl)
-   [Python](software.html#Python)
-   [R](software.html#R)
-   [Ruby](software.html#Ruby)
-   [Tcl/Tk](software.html#Tcl/Tk)

----------

Are the netCDF libraries thread-safe? {#Are-the-netCDF-libraries-thread-safe}
-----------------

The C-based libraries are *not* thread-safe. C-based libraries are those
that depend on the C library, which currently include all language
interfaces except for the Java interface. The Java interface is
thread-safe when a few simple rules are followed, such as each thread
getting their handle to a file.

----------

How does the C++ interface differ from the C interface? {#How-does-the-Cpp-interface-differ-from-the-C-interface}
-----------------

It provides all the functionality of the C interface (except for the
generalized mapped access of ncvarputg() and ncvargetg()) and is
somewhat simpler to use than the C interface. With the C++ interface, no
IDs are needed for netCDF components, there is no need to specify types
when creating attributes, and less indirection is required for dealing
with dimensions. However, the C++ interface is less mature and
less-widely used than the C interface, and the documentation for the C++
interface is less extensive, assuming a familiarity with the netCDF data
model and the C interface. Recently development of the C++ interface has
languished as resources have been redirected to enhancing the Java
interface.

----------

How does the Fortran interface differ from the C interface? {#How-does-the-Fortran-interface-differ-from-the-C-interface}
-----------------

It provides all the functionality of the C interface. The Fortran
interface uses Fortran conventions for array indices, subscript order,
and strings. There is no difference in the on-disk format for data
written from the different language interfaces. Data written by a C
language program may be read from a Fortran program and vice-versa. The
Fortran-90 interface is much smaller than the FORTRAN 77 interface as a
result of using optional arguments and overloaded functions wherever
possible.

----------

How do the Java, Perl, Python, Ruby, ... interfaces differ from the C interface? {#How-do-the-Java-Perl-Python-Ruby-interfaces-differ-from-the-C-interface}
-----------------

They provide all the functionality of the C interface, using appropriate
language conventions. There is no difference in the on-disk format for
data written from the different language interfaces. Data written by a C
language program may be read from programs that use other language
interfaces, and vice-versa.

----------

How do I handle errors in C? {#How-do-I-handle-errors-in-C}
-----------------

For clarity, the NetCDF C Interface Guide contains examples which use a
function called handle\_err() to handle potential errors like this:

~~~~ {.boldcode}
     status = nc_create("foo.nc", NC_NOCLOBBER, &ncid);
     if (status != NC_NOERR) handle_error(status);
~~~~

Most developers use some sort of macro to invoke netCDF functions and
test the status returned in the calling context without a function call,
but using such a macro in the User's Guides arguably makes the examples
needlessly complex. For example, some really excellent developers define
an "ERR" macro and write code like this:

~~~~ {.boldcode}
         if (nc_create(testfile, NC_CLOBBER, &ncid)) ERR;
~~~~

where Err is defined in a header file:

~~~~ {.boldcode}
/* This macro prints an error message with line number and name of
 * test program. */
#define ERR do { \
fflush(stdout); /* Make sure our stdout is synced with stderr. */ \
err++; \
fprintf(stderr, "Sorry! Unexpected result, %s, line: %d\n", \
        __FILE__, __LINE__);                                \
} while (0)
~~~~

Ultimately, error handling depends on the application which is calling
netCDF functions. However we strongly suggest that some form of error
checking be used for all netCDF function calls.

----------


CMake-Related Frequently Asked Questions {#cmake_faq}
========================================

Below are a list of commonly-asked questions regarding NetCDF and CMake.

How can I see the options available to CMake? {#listoptions}
---------------------------------------------

        $ cmake [path to source tree] -L	- This will show the basic options.
        $ cmake [path to source tree] -LA	- This will show the basic and advanced options.


How do I specify how to build a shared or static library? {#sharedstatic}
--------------------------------------------------------

    This is controlled with the internal `cmake` option, `BUILD_SHARED_LIBS`.

        $ cmake [Source Directory] -DBUILD_SHARED_LIBS=[ON/OFF]


Can I build both shared and static libraries at the same time with cmake? {#sharedstaticboth}
-------------------------------------------------------------------------

Not at this time; it is required to instead build first one version, and then the other, if you need both.

How can I specify linking against a particular library? {#partlib}
-------------------------------------------------------

It depends on the library.  To specify a custom `ZLib`, for example, you would do the following:

        $ cmake [Source Directory] -DZLIB_LIBRARY=/path/to/my/zlib.lib


`HDF5` is more complex, since it requires both the `hdf5` and `hdf5_hl` libraries. You would specify custom `HDF5` libraries as follows:

        $ cmake [Source Directory] -DHDF5_LIB=/path/to/hdf5.lib \
            -DHDF5_HL_LIB=/path/to/hdf5_hl.lib \
            -DHDF5_INCLUDE_DIR=/path/to/hdf5/include


Alternatively, you may specify:

        $ cmake [Source Directory] \
            -DHDF5_LIBRARIES="/path/to/hdf5.lib;/path/to/hdf5_hl.lib" \
            -DHDF5_INCLUDE_DIRS=/path/to/hdf5/include/


What if I want to link against multiple libraries in a non-standard location {#nonstdloc}
----------------------------------------------------------------------------

    You can specify the path to search when looking for dependencies and header files using the `CMAKE_PREFIX_PATH` variable:

* Windows:

        $ cmake [Source Directory] -DCMAKE_PREFIX_PATH=c:\shared\libs\


* Linux/Unix/OSX:

        $ cmake [Source Directory] -DCMAKE_PREFIX_PATH=/usr/custom_library_locations/

How can I specify a Parallel Build using HDF5 {#parallelhdf}
----------------------------------------------

If cmake is having problems finding the parallel `HDF5` install, you can specify the location manually:


        $ cmake [Source Directory] -DENABLE_PARALLEL=ON \
            -DHDF5_LIB=/usr/lib64/openmpi/lib/libhdf5.so \
            -DHDF5_HL_LIB=/usr/lib64/openmpi/lib/libhdf5.hl.so \
            -DHDF5_INCLUDE_DIR=/usr/include/openmpi-x86_64 \

You will, of course, need to use the location of the libraries specific to your development environment.

----------------

Plans {#Plans}
================

What other future work on netCDF is planned? {#What-other-future-work-on-netCDF-is-planned}
-----------------

Issues, bugs, and plans for netCDF are maintained in the Unidata issue
tracker sites for
[netCDF-C](https://www.unidata.ucar.edu/jira/browse/NCF), [Common Data Model / NetCDF-Java](https://www.unidata.ucar.edu/jira/browse/CDM),
[netCDF-Fortran](https://www.unidata.ucar.edu/jira/browse/NCFORTRAN),
and [netCDF-CXX4](https://www.unidata.ucar.edu/jira/browse/NCXXF), and
[old netCDF-C++
(deprecated)](https://www.unidata.ucar.edu/jira/browse/NCCPP).
