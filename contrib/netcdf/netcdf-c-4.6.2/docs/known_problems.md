# Known Problems with netCDF {#known_problems}


Known Problems with netCDF 4.3.0
--------------------------------

-   [clang compiler (default on OSX 10.9 Mavericks) detects error
    building ncgen3](#clang-ncgen3)

Known Problems with netCDF 4.2
------------------------------

-   [Fortran options of nc-config utility (--fflags, --flibs, --has-f90)
    don't work correctly](#nc-config-fortran)
-   [Using "--with-hdf5=..." configure option doesn't seem to
    work](#with-hdf5)

Known Problems with netCDF 4.1.3
--------------------------------

-   [nccopy -d and -c options for compression and chunking don't work on
    netCDF-4 input files](#nccopy-compression-bug)
-   [Debug statement left in F90 source](#f90-debug-segfault)
-   [Building with Intel Fortran on Mac OS X](#intel-fortran-macosx)
-   [Accessing OPeNDAP servers using a constraint expression](#dap-413)
-   [Configuring with "--enable-benchmarks"
    option](#enabling-benchmarks)

Known Problems with netCDF 4.1.2
--------------------------------

-   [Building with Intel Fortran on Mac OS X](#intel-fortran-macosx)
-   [Problem with disabling fill mode when using Lustre (or other large
    blksize file system)](#lustre)

Known Problems with netCDF 4.1.1
--------------------------------

-   [Ncgen is known to produce bad output and should not
    be used.](#bad-ncgen)
-   [Building with Intel Fortran on Mac OS X](#intel-fortran-macosx)
-   ["make check" fails when linked with HDF5-1.8.6](#incompat-411-186)
-   [Make tries to regenerate documentation with texi2dvi
    command](#texi2dvi)
-   [Accessing a multidimensional variable with more than 4 billion
    values on a 32-bit platform](#big-mvar-32bit)

------------------------------------------------------------------------

### The clang compiler (default on OSX 10.9 Mavericks) detects error building ncgen3

Building the netCDF C library with the clang C compiler, the default
/usr/bin/cc on OSX 10.9 Mavericks, detects an error in compiling
ncgen3/load.c. A fix is to insert the line

    #include <config.h>

above the "`#include <stdlib.h>`" statement near the beginning of
ncgen3/genlib.h.

This fix will be in the next release.

### Fortran options of nc-config utility (--fflags, --flibs, --has-f90) don't work correctly

Beginning with version 4.2 of the C-based netCDF software, the netCDF
Fortran library is built from an independent netcdf-fortran release with
its own nf-config utility. In netCDF-4.2 the nc-config utility doesn't
detect whether nf-config is installed and make use of its output to
preserve backward compatibility with nc-config from previous releases.
This problem is fixed in netCDF-4.2.1-rc1 and later releases.

### Using "--with-hdf5=..." configure option doesn't seem to work

With releases of netCDF-4 after version 4.1.2 (this includes 4.1.3, 4.2,
4.2.1, ...) you don't use "--with-hdf5" to specify the location of the
HDF5 libraries, you use CPPFLAGS and LDFLAGS, as in

    CPPFLAGS=-I/usr/local/hdf5/include LDFLAGS=-L/usr/local/hdf5/lib ./configure

The reason for this change is explained
[here](https://www.unidata.ucar.edu/jira/browse/NCF-20).

### nccopy -d and -c options for compression and chunking don't work on netCDF-4 input files

Due to a bug in nccopy, the "-d n" and "-c" options only work for
classic and 64-bit input files, producing netCDF-4 classic model output
files. These options are also useful for netCDF-4 files to compress or
recompress files and to chunk or rechunk variables. The bug description
and its resolution have been
[entered](https://www.unidata.ucar.edu/jira/browse/NCF-79) into the
issue tracker.

The bug has been fixed in all releases since 4.1.3, including the
netcdf-4.2-rc1 release candidate.

### Debug statement left in F90 source

The debugging statement

            print *, values(1, 1), values(1, 2), values(1, 3), values(1, 4)

was inadvertently left in the file f90/netcdf\_expanded.f90 at line 734,
and should be removed. If the variable has a second dimension less than
4, this can cause a segfault. The problem has been fixed in the
subsequent netcdf-fortran-4.2 release.

### Ncgen is known to produce bad output.

Dave Allured at NOAA has reported that the ncgen for 4.1.1 produces bad
.nc files under circumstances. We recommend that this version of ncgen
should not be used.

### Building with Intel Fortran on Mac OS X

Setting the environment variable **lt\_cv\_ld\_force\_load=no** before
invoking the configure script is a workaround to successfully build
netCDF version 4.1.3 with the Intel ifort Fortran compiler.
Specifically, the following works on Mac OS X 10.7.x (Lion) for building
C and Fortran libraries and passing all tests on Lion:

    $ lt_cv_ld_force_load=no FC="ifort" CC="cc" CXX="" \
      LDFLAGS=-L/WHERE_HDF5_IS_INSTALLED/lib \
      CPPFLAGS=-I/WHERE_HDF5_IS_INSTALLED/include ./configure
      make check

(The CXX environment variable is set to "" in this example to disable
building and testing the legacy netCDF-3 C++ API, because of an as yet
unsolved error that's not relevant to this Fortran problem.)

### Accessing OPeNDAP servers using a constraint expression

The use of subsetting by specifying a URL with subsetting information to
dap-enabled netCDF is broken for stable release 4.1.3. This can be
demonstrated by using the 4.1.3 version of ncdump to access data from an
OPeNDAP server, using a constraint expression in the URL, which results
in the error message

    NetCDF: Malformed or inaccessible DAP DDS

This bug is fixed in 4.2 releases after 2011-09-11, as well as by fixing
the 4.1.3 release using the 3 replacement source files in [this tar
file](http://www.unidata.ucar.edu/downloads/netcdf/ftp/4.1.3-fix.tar).

### Configuring with "--enable-benchmarks" option

Using the "--enable-benchmarks" option to the configure script fails
with a compile error on some platforms, and has been reported to fail in
other ways on some other platforms.

### Problem with disabling fill mode when using Lustre (or other large blksize file system)

Joerg Henrichs has reported a bug when writing netCDF classic format
data with fill mode disabled, on a file system such as Lustre that uses
a large disk block size (for example 2 MB). On such a system, tests run
by "make check" may all pass, but some write operations will fail
without reporting an error.

A fix is available. in netCDF-4.1.3-beta1 or later and also in [this
patch](/software/netcdf/patches/nofill-bug.patch) for the file
libsrc/posixio.c in netCDF-4.1.2 and earlier versions.

An example that failed depended on all the following circumstances:

-   The disk block size must be within a range of sizes, in this case
    one of the 58080 values between 2091953 and 2150032 inclusive.
    Equivalently, the file must be created using the "double underbar"
    function nc\_\_create() with a I/O block size in this range.
-   No-fill mode must be set, for example by calling nc\_set\_fill()
    with NC\_NOFILL, or by using a software package such as NCO that
    uses no-fill mode by default.
-   A data variable at the end of a file being created must be written
    in reverse order from how it is stored on disk, for example writing
    a multidimensional variable by slices in reverse order along one of
    its more slowly varying dimensions. This must result in writing at
    least two disk blocks beyond the end of the file.

Workarounds include avoiding use of no-fill mode (NC\_NOFILL), enabling
share mode (NC\_SHARE), changing the order of writes of a
multidimensional variable written in reverse order, or creating the file
using nc\_\_create with a blocksize outside the range in which erroneous
writes occur. Some of these workarounds slow the write performance of
netCDF.

### "make check" fails when linked with HDF5-1.8.6

When built with HDF5 version 1.8.6, version 4.1.1 fails one of the tests
invoked by "make check":

       ...
      *** Checking that one var, two dimscales, one att file can still be read by HDF5...ok.
      *** Creating a HDF5 file with one var and no dimension scales...ok.
      HDF5-DIAG: Error detected in HDF5 (1.8.6) thread 0:
        #000: H5O.c line 717 in H5Oget_info_by_idx(): group not found
          major: Symbol table
          minor: Object not found
        #001: H5Gloc.c line 591 in H5G_loc_find_by_idx(): can't find object
          major: Symbol table
          minor: Object not found
       ...
      PASS: tst_endian_fill
      ================================================
      1 of 59 tests failed
      Please report to support-netcdf@unidata.ucar.edu
      ================================================
      make[2]: *** [check-TESTS] Error 1

Currently the workarounds are to either

-   Use the earlier HDF5 1.8.5-patch1 release for building netCDF 4.1.1
-   Use netCDF-4.1.2-beta2 or later with HDF5 1.8.6

The HDF5 1.8.5-patch1 release is available from the HDF5 site at
<http://www.hdfgroup.org/ftp/HDF5/prev-releases/> or from the netCDF-4
ftp site at <ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4>.

### Make tries to regenerate documentation with texi2dvi command

After building netCDF-4.1.1, invoking "make clean", and then building it
again with "make all" or "make check", a failure to find the texi2dvi
command is reported:

    make[1]: Entering directory `/usr/local/netcdf/netcdf-4.1.1/man4'
    TEXINPUTS=".:$TEXINPUTS" \
        MAKEINFO='/bin/sh /usr/local/netcdf/netcdf-4.1.1/missing --run makeinfo   -I .' \
        texi2dvi -s  --pdf --batch netcdf.texi
    make[1]: Leaving directory `/usr/local/netcdf/netcdf-4.1.1/man4'
    /bin/sh: texi2dvi: command not found

This results from a bug where "make clean" erroneously deleted the
documentation generated for the release, so make tries to regenerate the
documentation from source, which will fail if you don't happen to have
the "texi2dvi" program installed (which you shouldn't need).

This is fixed in the current snapshot and in the upcoming release 4.1.2,
but a workaround is to get a new copy of the 4.1.1 source and rebuild
from that with the same settings you used to get to the above message,
without invoking "make clean" until after the software and documentation
is successfully installed. An alternative workaround is to just invoke
"make install" after the error above and use online documentation.

### Accessing a multidimensional variable with more than 4 billion values on a 32-bit platform

Kari Hoijarvi has reported a bug in implementation of large variable
support that has been in netCDF software since at least 1997, and which
is fixed in netCDF snapshots after 2010-05-13 as well as in the upcoming
netCDF-4.1.2 release. The bug occurs when all of the following
conditions are met:

-   A 32-bit version of the netCDF library is used.
-   The file format is either classic or 64-bit offset format.
-   Values are written to or read from a fixed-size variable that has
    more than 2^32^ (4,294,967,296) values, or a record variable with
    more than 2^32^ values per record. The variable must be the last
    fixed-size variable in a file with no record variables or the last
    record variable, because otherwise it would violate the format
    constraints for netCDF classic or 64-bit offset formats described
    [here](http://www.unidata.ucar.edu/netcdf/docs/netcdf/NetCDF-Classic-Format-Limitations.html).
    Note that the bug involves number of values, not bytes, so if the
    variable is of type integer or float, for example, it would require
    more than 17 Gbytes.
-   The variable must have 2 or more dimensions.
-   The values to be written or read must begin after the first 2^32^
    values of the variable.

In this case an undetected integer overflow occurred in calculating the
file offset, and the values were written to or read from the wrong
location in the file, overwriting data stored at that location in the
case of a write.

The fix is a one-line change to a line in the libsrc/putget.m4 file,
from which the libsrc/putget.c file is generated, replacing the
statement

                lcoord += *up * *ip;

with

                lcoord += (off_t)(*up) * (off_t)(*ip);

Known Problems with netCDF 4.0.1
--------------------------------

-   [Including mpi.h before netcdf.h breaks MPI](#include_mpi_order)
-   [With Sun C compiler, 64-bit ncdump fails](#suncc-m64-ncdump)
-   [Portland Group compilers can't build shared fortran 90 or shared
    C++ library, or 64-bit C++ API](#PG_shared)
-   [Intel 10.1 64-bit C++ compiler problem](#intel_10.1_64_bit_c++)
-   [Intel 9.1 C++ compiler problem doesn't build C++
    API](#intel_9.1__c++)
-   [ncgen/ncdump test failure with Intel version 11
    compilers](#intel_11-ncgen)
-   ["ncdump -v group/var" reports "group not found"](#ncdump-v)

### Including mpi.h before netcdf.h breaks MPI

Luis Kornblueh reports a subtle bug in netcdf 4.0.1. In the netcdf.h
header file, the following mpi entities are defined:

     /* These defs added by netCDF configure because parallel HDF5 is not 
     present. */
     #define MPI_Comm int
     #define MPI_Info int
     #define MPI_COMM_WORLD 0
     #define MPI_INFO_NULL 0

If mpi.h is included before netcdf.h, these defines (may) break the MPI
implementation.

### With Sun C compiler, 64-bit ncdump fails

As identified by Udo Grabowski, using the "-m64" option to build netCDF
with the Sun C compiler results in a failed test when running "make
check" in the ncdump directory:

    *** checking that test1.cdl and test2.nc are the same...
    62,63c62,63
    <   8.88178419700125e-16, 1.11022302462516e-15, 1.33226762955019e-15,
    <   1.55431223447522e-15, 1.77635683940025e-15, 222044604925031 ;
    ---
    >   1.97215226305253e-31, 2.46519032881567e-31, 2.9582283945788e-31,
    >   3.45126646034193e-31, 3.94430452610506e-31, 0.0493038065763132 ;
    FAIL: run_tests.sh

This bug is fixed in recent Sun C compiler releases, for example "Sun C
5.11 SunOS\_i386 Aten 2010/05/10".

Short of upgrading the compiler, other workarounds include specifying

    CFLAGS="-O0 -m64" 

before rerunning the configure script, to turn off optimization, or just
install an ncdump built without "-m64". Because ncdump reads only a
little data at a time, there is no benefit to a 64-bit ncdump. The
32-bit ncdump handles classic, 64-bit offset, and netCDF-4 files
correctly even if they are larger than 4 GiB.

### Portland Group compilers can't build shared fortran 90 library or shared C++ library

The portland group compilers can't build netCDF shared fortran 90
library. They fail with this error:

    pgf90  -I../fortran -I../f90      -I../libsrc -I../fortran   -I../f90
    -g -c -o tst_f90.o tst_f90.f90
    /bin/sh ../libtool   --mode=link pgf90  -I../fortran -I../f90 --
                 ---I../libsrc -I../fortran   -I../f90 -g -L/lib --
                 ---o tst_f90 tst_f90.o ../fortran/libnetcdff.la --
                 ---lm     ../libsrc/libnetcdf.la  
    libtool: link: pgf90 -I../fortran -I../f90 -I../libsrc -I../fortran
    -I../f90 -g -o .libs/tst_f90 tst_f90.o  -L/lib
    ../fortran/.libs/libnetcdff.so -lm ../libsrc/.libs/libnetcdf.so
    -Wl,-rpath
    -Wl,/machine/netcdf/n362_test_9456/netcdf-3.6.3-snapshot2008081305/install/lib
    tst_f90.o:(.debug_info+0x135d): undefined reference to
    `..Dm_typesizes'
    tst_f90.o:(.debug_info+0x136e): undefined reference to `..Dm_netcdf'

If anyone could shed some light on this, it would be most appreciated.
Send comments to support-netcdf@unidata.ucar.edu.

The C++ compiler chokes on the netCDF C++ tests on a shared build:

           pgCC -DHAVE_CONFIG_H -I. -I.. -I../fortran -I../libdap
           -I../libsrc   -g -c -o tst_failure.o tst_failure.cpp
    /bin/sh ../libtool --tag=CXX   --mode=link pgCC  -g     -o tst_failure
           tst_failure.o ../cxx/libnetcdf_c++.la  ../libsrc/libnetcdf.la
    libtool: link: pgCC -g -o .libs/tst_failure tst_failure.o
           ../cxx/.libs/libnetcdf_c++.so ../libsrc/.libs/libnetcdf.so
           -Wl,--rpath -Wl,/usr/local/lib
    make[2]: Leaving directory `/machine/shecky/n4_new2/cxx'
    make  check-TESTS
    make[2]: Entering directory `/machine/shecky/n4_new2/cxx'
    C++ runtime abort: internal error: static object marked for
           destruction more than once
    /bin/sh: line 4:  8445 Aborted                 ${dir}$tst
    FAIL: nctst
    C++ runtime abort: internal error: static object marked for
           destruction more than once
    /bin/sh: line 4:  8468 Aborted                 ${dir}$tst
    XFAIL: tst_failure

There is a problem with the pgCC compiler noted as "Fixed in version
6.2.1" described as

      C++ runtime abort: internal error: static object marked for
      destruction more than once

here as [Technical Problem Report
3809](http://www.pgroup.com/support/tprs_62.htm#t3809).

This bug was also previously [reported by a
user](http://www.unidata.ucar.edu/support/help/MailArchives/netcdf/msg03783.html).

### Intel 10.1 64-bit C++ compiler problem

On my test machine, the intel 10.1 C++ compiler cannot build the netCDF
C++ API in 64-bit mode. I get an error like this:

    make[1]: Entering directory
    `/machine/netcdf/n362_test_16030/netcdf-3.6.3-snapshot2008081312/cxx'
    depbase=`echo netcdf.lo | sed 's|[^/]*$|.deps/&|;s|\.lo$||'`;\
    /bin/sh ../libtool --tag=CXX   --mode=compile icpc -DHAVE_CONFIG_H
    -I. -I.. -I../fortran -I../libdap      -I../libsrc
    -I/opt/intel/cce/10.1.015/include -MT netcdf.lo -MD -MP -MF
    $depbase.Tpo -c -o netcdf.lo netcdf.cpp &&\
    mv -f $depbase.Tpo $depbase.Plo
    libtool: compile:  icpc -DHAVE_CONFIG_H -I. -I.. -I../fortran
    -I../libdap -I../libsrc -I/opt/intel/cce/10.1.015/include -MT
    netcdf.lo -MD -MP -MF .deps/netcdf.Tpo -c netcdf.cpp -o netcdf.o
    /usr/include/c++/4.3.0/x86_64-redhat-linux/bits/c++locale.h(94):
    error: argument of type "__va_list_tag *" is incompatible with
    parameter of type "char *"
          const int __ret = __builtin_vsnprintf(__out, __size, __fmt,
          __args);
                                                                      ^

    compilation aborted for netcdf.cpp (code 2)
    make[1]: *** [netcdf.lo] Error 1
    make[1]: Leaving directory
    `/machine/netcdf/n362_test_16030/netcdf-3.6.3-snapshot2008081312/cxx'
    make: *** [check-recursive] Error 1

This is because the Intel C++ compiler has not caught up to the GNU C++
compiler, and for some reason that is not clear to me, it is using the
header files from gcc.

To solve this problem, install an older version of gcc (4.1.2 works in
testing at Unidata World Test Center, located at the bottom of a
six-mile deep mine shaft.) Put the bin directory at the beginning of
your PATH, and the lib (or lib64) directory at the beginning at the
LD\_LIBRARY\_PATH. Then rebuild.

### Intel 9.1 C++ compiler problem doesn't build C++ API

On my test machine, the intel 9.1 C++ compile fails like this:

    make  nctst tst_failure
    make[2]: Entering directory
    `/machine/netcdf/n362_test_16035/netcdf-3.6.3-snapshot2008081312/cxx'
    depbase=`echo nctst.o | sed 's|[^/]*$|.deps/&|;s|\.o$||'`;\
    icpc -DHAVE_CONFIG_H -I. -I.. -I../fortran -I../libdap
    -I../libsrc   -I/opt/intel/cc/9.1.047/include/c++/ -MT nctst.o -MD -MP
    -MF $depbase.Tpo -c -o nctst.o nctst.cpp &&\
    mv -f $depbase.Tpo $depbase.Po
    /bin/sh ../libtool --tag=CXX   --mode=link icpc
    -I/opt/intel/cc/9.1.047/include/c++/     -o nctst nctst.o
    ../cxx/libnetcdf_c++.la  ../libsrc/libnetcdf.la  
    libtool: link: icpc -I/opt/intel/cc/9.1.047/include/c++/ -o nctst
    nctst.o  ../cxx/.libs/libnetcdf_c++.a ../libsrc/.libs/libnetcdf.a
    nctst.o: In function `main':
    nctst.cpp:(.text+0x22e): undefined reference to
    `std::ios_base::clear(std::_Iosb::_Iostate, bool)'
    nctst.cpp:(.text+0x290): undefined reference to
    `std::ios_base::clear(std::_Iosb::_Iostate, bool)'
    nctst.cpp:(.text+0x3b1): undefined reference to
    `std::ios_base::clear(std::_Iosb::_Iostate, bool)'
    nctst.cpp:(.text+0x520): undefined reference to
    `std::ios_base::clear(std::_Iosb::_Iostate, bool)'
    nctst.cpp:(.text+0x582): undefined reference to
    `std::ios_base::clear(std::_Iosb::_Iostate, bool)'
    nctst.o:nctst.cpp:(.text+0x767): more undefined references to
    `std::ios_base::clear(std::_Iosb::_Iostate, bool)' follow
    nctst.o: In function `__sti__$E':
    nctst.cpp:(.text+0x2cb0): undefined reference to
    `std::_Winit::_Winit()'
    nctst.cpp:(.text+0x2cbf): undefined reference to `std::_Winit::~_Winit()'

Anyone who can shed light on this should send email to
support-netcdf@unidata.ucar.edu.

### ncgen/ncdump test failure with Intel version 11 compilers

Ed Anderson reports that the tests of the netcdf-4.0 (and presumable
4.0.1 and 3.6.3) package fail with the recently released version 11 of
the Intel compilers, producing the error message:

    *** creating UTF-8 test file tst_utf8.nc...Sorry! Unexpected result, tst_utf8.c, line: 63
    Sorry! Unexpected result, tst_utf8.c, line: 68
    Sorry! Unexpected result, tst_utf8.c, line: 72
    Sorry! Unexpected result, tst_utf8.c, line: 79
    Sorry! Unexpected result, tst_utf8.c, line: 90
    Sorry! Unexpected result, tst_utf8.c, line: 92
    Sorry! Unexpected result, tst_utf8.c, line: 114
    Sorry! Unexpected result, tst_utf8.c, line: 117
    Sorry! Unexpected result, tst_utf8.c, line: 119
    Sorry! Unexpected result, tst_utf8.c, line: 122
    /bin/sh: line 1:  5216 Segmentation fault      ${dir}$tst
    FAIL: tst_utf8

    *** Testing ncgen and ncdump for UTF8 support...
    *** creating classic offset file with utf8 characters...
    ncgen: NetCDF: Name contains illegal characters
    ncgen: NetCDF: Invalid dimension ID or name
    ncgen: NetCDF: Name contains illegal characters
    FAIL: run_utf8_tests.sh
    =========================================
    2 of 8 tests failed

Ed also reports this is a compiler problem (which has been reported) and
that there is a workaround:

    ... in libsrc/string.c the test

    if(ch <= 0x7f) {

    can be changed (in two places) to

    if(ch < 0x7f || ch == 0x7f) {

    This was the only change I needed to pass the netcdf-4 tests with Intel
    version 11.

### "ncdump -v group/var" reports "group not found"

John Storrs reported a bug using ncdump -v applied to netCDF-4 files, in
which an erroneous 'group not found' message was displayed for valid
group/var names. This is fixed in the next release, and the fix is also
in the [current snapshot
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/snapshot/).

Known Problems with netCDF 4.0
------------------------------

-   [Ncdump assumes default fill value for unsigned byte
    data](#ncdump_ubyte_fill)
-   [Ncdump of compound type with array field](#compound_array_field)
-   [Memory leak with VLEN attributes](#mem_leak)
-   [Error dyld: Symbol not found:
    _H5P_CLS_FILE_ACCESS_g](#o_problem_mac)
-   [Bug with multiple unlimited dimensions in one
    var](#multiple_unlimited)
-   [Fortran90 interface Using Intel ifort under
    Cygwin](#ifort-f90-cygwin)
-   [ncdump bug for filenames beginning with a numeric
    character](#ncdump-numeric-filename)
-   [ncgen/ncdump test failure with Intel version 11
    compilers](#intel_11-ncgen)

### Ncdump assumes default fill value for unsigned byte data

The ncdump utility incorrectly assumes a default fill value of "255" for
data of unsigned byte type, although no default fill value is assumed
for data of type signed byte. There should be no default fill values
when reading any byte type, signed or unsigned, because the byte ranges
are too small to assume one of the values should appear as a missing
value unless a \_FillValue attribute is set explicitly. This bug is
fixed in the current snapshot distribution.

### Ncdump of compound type with array field

Running the ncdump utility on a file with a compound type with an array
field may result in a segmentation violation. A fix is in the current
netCDF-4.0 snapshot distribution.

### Memory leak with VLEN attributes

We believe there are some memory leaks associated with VLEN attributes
in HDF5 1.8.1. This is being addressed by the HDF5 team, and will be
fixed by the next HDF5 release.

### Error dyld: Symbol not found: _H5P_CLS_FILE_ACCESS_g

On some Macintosh systems here at NetCDF World Test Center, on the
hundreth floor of UCAR Tower \#2, the following build error occurs:

    *** Checking HDF5 enum types.
    *** Checking simple HDF5 enum type...ok.
    *** Checking HDF5 enum type missing values...ok.
    *** Tests successful!
    PASS: tst_h_enums
    dyld: Symbol not found: _H5P_CLS_FILE_ACCESS_g
      Referenced from:
      /tmp/n4_sid/netcdf-4.0-snapshot2008042320/libsrc4/.libs/libnetcdf.5.dylib
      Expected in: flat namespace

    FAIL: tst_lists
    dyld: Symbol not found: _H5P_CLS_FILE_ACCESS_g
      Referenced from:
      /tmp/n4_sid/netcdf-4.0-snapshot2008042320/libsrc4/.libs/libnetcdf.5.dylib
      Expected in: flat namespace

    FAIL: tst_dims
    dyld: Symbol not found: _H5P_CLS_FILE_ACCESS_g
      Referenced from:
      /tmp/n4_sid/netcdf-4.0-snapshot2008042320/libsrc4/.libs/libnetcdf.5.dylib
      Expected in: flat namespace

    etc.

This can be caused by the configure script failing to add "-lhdf5" to
the link flags in the generated Makefiles. Set LDFLAGS to include
"-L/WHERE/HDF5/IS/INSTALLED/lib -lhdf5" and try again.

------------------------------------------------------------------------

Bug with multiple unlimited dimensions in one var

There is a bug in the 4.0 release related to the lengths of dimensions
when more than one unlimited dimension is used in the same variable.

The bug is fixed in the latest [netCDF-4 snapshot
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/snapshot/netcdf-4-daily.tar.gz).

### Fortran90 interface Using Intel ifort under Cygwin

Chris Dallimore reports success in getting the Fortran 90 interface of
Version 4.0 to compile under CYGWIN using the Intel ifort compile;

    1 - Download and unpack netcdf-4.0.tar.gz

    2 - In configure replace conftest.o and conftestf.o with conftest. 
    $ac_objext and conftest.$ac_objext, I'm Not sure why autoconf doesn't  
    do this.

    3 -
    Save http://msinttypes.googlecode.com/svn/trunk/inttypes.h as libsrc/ 
    inttypes_msvc.h
    Save ttp://msinttypes.googlecode.com/svn/trunk/stdint.h as libsrc/ 
    stdint_msvc.h

    modify line 43 of libsrc/inttypes_msvc.h
    from
    #include 
    to
    #include 

    4 - in libsrc utf8proc.h at line 79 replaces

    #include 

    with

    #ifndef _MSC_VER
    #include 
    #else
    #include 
    typedef long ssize_t;
    typedef unsigned int uint32_t;
    #endif

    It looks like configure is checking for ssize_t so there is probably a  
    better way to do this.

    5 -
    in libsrc/posixio.c at line 18 replace

    #ifdef _MSC_VER /* Microsoft Compilers */
    #include 
    #else
    #include 
    #endif

    with

    #ifdef _MSC_VER /* Microsoft Compilers */
    #include 
    typedef long ssize_t;
    typedef unsigned int uint32_t;
    #else
    #include 
    #endif

    6 -
    in putget.m4 at line 24 added

    #ifdef _MSC_VER
    #include 
    #endif // _MSC_VER ]

    Run
    ./configure --prefix=/cygdrive/z/cwr/Software/Eclipse/CWRModelSource -- 
    disable-examples  --disable-cxx --disable-utilities

    Relevant environment variables
    FC=ifort
    F90=ifort
    FFLAGS=/debug:full /traceback  /nologo
    FCFLAGS=/debug:full /traceback  /nologo
    FCFLAGS_f90=/debug:full /traceback  /nologo
    FLIBS=
    CXX=
    CPPFLAGS=/D AbsoftProFortran /D _MSC_VER /nologo

    IFORT_COMPILER10_POSIX=/cygdrive/c/Program Files/Intel/Compiler/ 
    Fortran/10.1.013
    VSINSTALLDIR=C:\Program Files\Microsoft Visual Studio 8
    INTEL_LICENSE_FILE=C:\Program Files\Common Files\Intel\Licenses
    PROCESSOR_IDENTIFIER=x86 Family 6 Model 15 Stepping 6, GenuineIntel
    TERM=cygwin
    WINDIR=C:\WINDOWS
    MAKEFLAGS=w -- F90=ifort
    VS80COMNTOOLS=C:\Program Files\Microsoft Visual Studio 8\Common7\Tools\
    VSINSTALLDIR_POSIX=/cygdrive/c/Program Files/Microsoft Visual Studio 8
    F90FLAGS=/debug:full /traceback  /nologo
    OS=CYGWIN
    MAKEOVERRIDES=${-*-command-variables-*-}
    USER=dallimor
    VCINSTALLDIR=C:\Program Files\Microsoft Visual Studio 8\VC
    INTEL_SHARED=C:\Program Files\Common Files\Intel\Shared Files
    LIB=C:\Program Files\Intel\Compiler\Fortran\10.1.013\Ia32\Lib;C: 
    \Program Files\Microsoft Visual Studio 8\VC\atlmfc\lib;C:\Program Files 
    \Microsoft Visual Studio 8\VC\lib;C:\Program Files\Microsoft Visual  
    Studio 8\VC\PlatformSDK\lib;C:\Program Files\Microsoft Visual Studio 
    \DF98\LIB;C:\Program Files\Microsoft Visual Studio\VC98\LIB
    IFORT_COMPILER10=C:\Program Files\Intel\Compiler\Fortran\10.1.013
    MFLAGS=-w
    VCINSTALLDIR_POSIX=/cygdrive/c/Program Files/Microsoft Visual Studio 8/ 
    VC
    PROCESSOR_LEVEL=6
    PATH=/cygdrive/c/Program Files/Intel/Compiler/Fortran/10.1.013/Ia32/ 
    Bin:/cygdrive/c/Program Files/Common Files/Intel/Shared Files/Ia32/ 
    Bin:/cygdrive/c/Program Files/Microsoft Visual Studio 8/Common7/IDE:/ 
    cygdrive/c/Program Files/Microsoft Visual Studio 8/VC/bin:/cygdrive/c/ 
    Program Files/Microsoft Visual Studio 8/Common7/Tools:/cygdrive/c/ 
    Program Files/Microsoft Visual Studio 8/Common7/Tools/bin:/cygdrive/c/ 
    Program Files/Microsoft Visual Studio 8/VC/PlatformSDK/bin:/cygdrive/c/ 
    apache-ant-1.7.0/bin:/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin:/ 
    cygdrive/c/Program Files/Microsoft Visual Studio/Common/Tools:/ 
    cygdrive/c/Program Files/Microsoft Visual Studio/Common/Msdev98/BIN:/ 
    cygdrive/c/Program Files/Microsoft Visual Studio/DF98/BIN:/cygdrive/c/ 
    Program Files/Microsoft Visual Studio/VC98/BIN:/cygdrive/c/Sun/SDK/jdk/ 
    bin:C:/oraclexe/app/oracle/product/10.2.0/server/bin:/cygdrive/c/ 
    WINDOWS/system32:/cygdrive/c/WINDOWS:/cygdrive/c/WINDOWS/System32/ 
    Wbem:/cygdrive/c/Program Files/jEdit:/cygdrive/c/Program Files/ 
    Microsoft SQL Server/90/Tools/binn/:/cygdrive/c/Program Files/Intel/ 
    Compiler/Fortran/10.1.013/IA32/Lib:/cygdrive/c/Program Files/Intel/ 
    Compiler/Fortran/10.1.013/EM64T/Lib:/cygdrive/c/Program Files/MATLAB/ 
    R2008a/bin:/cygdrive/c/Program Files/MATLAB/R2008a/bin/win32
    CPU=i386
    FP_NO_HOST_CHECK=NO
    PROCESSOR_ARCHITECTURE=x86
    PATHEXT=.COM;.EXE;.BAT;.CMD;.VBS;.VBE;.JS;.JSE;.WSF;.WSH
    INTEL_SHARED_POSIX=/cygdrive/c/Program Files/Common Files/Intel/Shared  
    Files
    MAKE_MODE=unix
    INFOPATH=/usr/local/info:/usr/share/info:/usr/info:
    PROGRAMFILES=C:\Program Files
    CC=cl
    INCLUDE=C:\Program Files\Intel\Compiler\Fortran 
    \10.1.013\Ia32\Include;C:\Program Files\Microsoft Visual Studio 8\VC 
    \atlmfc\include;C:\Program Files\Microsoft Visual Studio 8\VC 
    \include;C:\Program Files\Microsoft Visual Studio 8\VC\PlatformSDK 
    \include;C:\Program Files\Microsoft Visual Studio\DF98\INCLUDE;C: 
    \Program Files\Microsoft Visual Studio\VC98\INCLUDE


    Now configure works

    make will compile but fails on linking

    libtool: link: ( cd ".libs" && rm -f "libnetcdf2.la" && ln -s "../ 
    libnetcdf2.la" "libnetcdf2.la" )
    /bin/sh ../libtool --tag=CC   --mode=link /cygdrive/z/cwr/Software/ 
    Eclipse/CWRModelSource/src/external/netcdf_src/netcdf-4.0/compile  
    cl      -version-info 4:0:0   -o libnetcdf.la -rpath /cygdrive/z/cwr/ 
    Software/Eclipse/CWRModelSource/lib attr.lo ncx.lo putget.lo dim.lo  
    error.lo libvers.lo nc.lo string.lo v1hpg.lo var.lo utf8proc.lo   
    posixio.lo libnetcdf2.la ../fortran/libnetcdff.la
    libtool: link: warning: undefined symbols not allowed in i686-pc- 
    cygwin shared libraries
    libtool: link: (cd .libs/libnetcdf.lax/libnetcdf2.lib && ar x "/ 
    cygdrive/z/cwr/Software/Eclipse/CWRModelSource/src/external/netcdf_src/ 
    netcdf-4.0/libsrc/./.libs/libnetcdf2.lib")
    libtool: link: (cd .libs/libnetcdf.lax/libnetcdff.lib && ar x "/ 
    cygdrive/z/cwr/Software/Eclipse/CWRModelSource/src/external/netcdf_src/ 
    netcdf-4.0/libsrc/../fortran/.libs/libnetcdff.lib")
    .libs/libnetcdff.lax/libnetcdff90.lib/typeSizes.obj: No such file or  
    directory

    It looks like the Microsoft LInker doesn't like the GNU lib format.

    I was however able to compile and link using some static (ie non  
    automake) makefiles that are part of our overall model build  
    environment.

### ncdump bug for filenames beginning with a numeric character

The ncdump utility in releases 4.0 and 3.6.3 rejects filenames starting
with the digits 0,1 and 2 with an error message such as:

      ncdump: name begins with space or control-character: 2

This bug is fixed in the daily snapshot release and in 4.0.1-beta
releases, but a one-line patch to ncdump/dumplib.c (for either
netCDF-4.0 or netCDF-3.6.3) is to replace the line

     if((*cp >= 0x01 && *cp <= 0x32) || (*cp == 0x7f))

with the following line instead

     if((*cp >= 0x00 && *cp <= 0x20) || (*cp == 0x7f))

------------------------------------------------------------------------

Known Problems with netCDF 3.6.3
--------------------------------

-   [Building shared libraries on Macintosh with
    g95 fails.](#g95_mac_shared)
-   [Building Fortran/C++ shared libraries on AIX fails.](#AIX_shared)
-   [Building shared libraries on HPUX with native tools results in only
    static libraries.](#HPUX_shared)
-   [Can't build shared library with F90 API on IRIX.](#IRIX_f90_shared)
-   [ncdump bug for filenames beginning with a numeric
    character](#ncdump-numeric-filename)

------------------------------------------------------------------------

### Can't build shared library with F90 API on IRIX

When building shared libraries on out IRIX test system, I got the
following error:

    ld32: FATAL   12 : Expecting n32 objects: /lib/libc.so.1 is o32.

Obviously there is some ABI confusion here, but we don't know how to
resolve it. Any user who can solve this should email
support-netcdf@unidata.ucar.edu so that we can share the method with
other users.

------------------------------------------------------------------------

Known Problems with netCDF 3.6.2
--------------------------------

-   [Setting ARFLAGS does not work.](#ARFLAGS)
-   [Bugs in support for variables larger than 4 GiB](#large_vars_362)
-   [Bug in C++ interface prevents creating 64-bit offset format
    files](#cxx_64-bit)
-   [Shared libraries do not work with the NAG
    fortran compiler.](#nag_362)
-   [The tests in nf_test fail with seg fault with the Absoft Version
    10.0 fortran compiler.](#absoft10)
-   [The documented --enable-64bit option doesn't work](#enable-64bit)
-   [Building netCDF-3.6.2 with gfortran version
    4.3.x fails.](#gfortran_43)
-   [Building shared libraries on Macintosh with
    g95 fails.](#g95_mac_shared)
-   [Building shared libraries on HPUX with native tools results in only
    static libraries.](#HPUX_shared)
-   [Building Fortran/C++ shared libraries on AIX fails.](#AIX_shared)
-   [Building with older versions of g++ fails.](#old_gpp)
-   [The .NET build files are not included in the
    3.6.2 release.](#NET_3_6_2)
-   [Snapshot .NET build files do not work for Visual Studio 8.0
    beta releases.](#NET_80_362)
-   [The -disable-v2 option causes the fortran build to fail with some
    fortran compilers.](#disable-v2_3_6_2)
-   [The --disable-c option does not work.](#disable-c_3_6_2)

------------------------------------------------------------------------

### Setting ARFLAGS does not work

Sometimes when building netCDF, flags to the ar utility need to be set.
Setting ARFLAGS does not work.

(Note: If you are doing this to build 64-bit libraries on an AIX
platform, the most fool-proof way to built 64-bit applications under AIX
is to set the OBJECT\_MODE environment variable to 64. If you still feel
you must setr flags for ar, read on.)

Try the build again, setting AR\_FLAGS instead of ARFLAGS.

### Bugs in support for variables larger than 4 GiB

As first reported by Mario Emmenlauer, there is a bug in netCDF-3.6.2
(and earlier versions) in the code for creating byte and short type
variables greater than 4 GiB in size. This problem resulted in an
assertion violation or arithmetic exception that would have caused a
program to halt, rather than writing bad data or corrupting existing
data.

A fix is available as a [patch](../patches/large-vars-362-patch) to the
file libsrc/var.c in the netcdf-3.6.2 distribution. The bug is also
fixed in releases 3.6.3 and later.

1.  On 32-bit platforms (with size\_t an unsigned 32-bit type):
    -   For a short variable, if the product of dimensions (not counting
        the record dimension, if any) is greater than 2^31^ (that's
        2147483648), the following assertion violation occurs

            Assertion failed: remaining > 0, file putget.c, line 347

    -   For any type of variable, if the product of dimensions (not
        counting the record dimension, if any) is exactly 2^32^
        (that's 4294967296) or any multiple of 2^32^, a divide by zero
        occurs

            Arithmetic Exception(coredump)

2.  On 64-bit platforms (with size\_t an unsigned 64-bit type):
    -   For a byte variable, if the product of dimensions (not counting
        the record dimension, if any) is greater than 2^32^, an
        assertion violation occurs

            Assertion failed: *ulp <= X_SIZE_MAX, file ncx.c, line 1810

    -   For a short variable, if the product of dimensions (not counting
        the record dimension, if any) is greater than 2^31^, the same
        assertion violation occurs

            Assertion failed: *ulp <= X_SIZE_MAX, file ncx.c, line 1810

### Bug in C++ interface prevents creating 64-bit offset format files

As reported by Jos Verdoold, a bug in the netCDF 3.6.2 (and earlier
versions) C++ interface prevents creating new files in Offset64Bits mode
using the C++ API. The fix is to change two lines (378 and 393) in the
file src/cxx/netcdf.cpp, changing "=" to "|=" in each line, then
rebuild:

    378c378
    <   mode = NC_WRITE;
    ---
    >   mode |= NC_WRITE;
    393c393
    <   mode = NC_NOCLOBBER;
    ---
    >   mode |= NC_NOCLOBBER;

This fix has been incorporated into netCDF 3.6.3 and later versions.

### The tests in nf\_test fail with seg fault with the Absoft Version 10.0 fortran compiler.

The absoft fortran compiler, version 10.0, changes the way that a C
function returning string is called from fortran.

This causes the absoft fortran settings in cfortran.h to no longer work,
and this is reflected in a segmentation fault in the test program
nf\_test/nf\_test.

As a workaround, users with absoft version 10 can get the latest
netCDF-3 snapshot and build it with the --enable-absoft10-hack option
set.

Get the snapshot, and see the working output, on the [netCDF-3
snapshot](http://www.unidata.ucar.edu/software/netcdf/builds/snapshot/)
page.

### Shared libraries do not work with the NAG fortran compiler.

We have reports that the shared library build does not work with the NAG
fortran compiler. The NAG compiler is not one of the compilers we
current support (and test on) at Unidata. The only known work around is
to build without the --enable-shared option.

Any user who can debug this problem with the NAG compiler should send
the resuts to support-netcdf@unidata.ucar.edu, so that it can be
incorporated into the netCDF distribution.

Interested users may also wish to subscribe to the [netcdf-porting
mailing
list](http://www.unidata.ucar.edu/mailing_lists/archives/netcdf-porting/).

### The documented --enable-64bit option doesn't work.

The --enable-64bit option appeared in the 3.6.1 release, and was--
removed for the 3.6.2 release.

Unfortunately, the documentation was not updated, so that the 3.6.2
documentation still mentions the enable-64bit option. Sorry about that.

The documentation has been corrected for the [netCDF-3
snapshot](http://www.unidata.ucar.edu/software/netcdf/builds/snapshot/)
and the [netCDF-4
snapshot](http://www.unidata.ucar.edu/software/netcdf/builds/snapshot/index_4.html)
documentation.

### Building netCDF-3.6.2 with gfortran version 4.2.x or 4.3.x fails.

Something changed in gfortran version 4.3 relating to how fortran
functions can call C functions.

In netCDF, the interface between C and Fortran is handled by the
cfortran.h package, which requires a pre-processor define describing the
type of fortran you are using.

For gfortran up to version 4.1.x, the netCDF distribution builds cleanly
with the "gFortran" preprocessor symbol set. For gfortran 4.2.x and
greater, the "pgiFortran" preprocessor symbol works.

The 3.6.2 build uses "gFortran", unless you specifically set the
CPPFLAGS environmental variable to "-DpgiFortran".

This works in my bash shell:

    FC=gfortran CPPFLAGS=-DpgiFortran ./configure && make check

This problem has been fixed in the [netCDF-3
snapshot](../builds/snapshot). Now configure checks the version of
gfortran before setting the appropriate flag.

### Building shared libraries on Macintosh with g95 fails.

Building shared libraries on the Macintosh fails

     
    *** Testing netCDF-3 Fortran 90 API.
    dyld: lazy symbol binding failed: Symbol not found: __g95_size
      Referenced from:
      /tmp/n3_mort/netcdf-3.6.2-snapshot2007052501/fortran/.libs/libnetcdff.4.dylib
      Expected in: flat namespace

    dyld: Symbol not found: __g95_size
      Referenced from:
      /tmp/n3_mort/netcdf-3.6.2-snapshot2007052501/fortran/.libs/libnetcdff.4.dylib
      Expected in: flat namespace

    FAIL: tst_f90
    =========================================
    1 of 5 tests failed
    Please report to support@unidata.ucar.edu
    =========================================

### Building shared libraries on HPUX with native tools results in only static libraries.

On the only HPUX machine I have access to for testing, the
--enable-shared still results in only the static library being linked.

This may be because of and old C++ compiler on this particular platform.
Any HPUX use who can provide information about this should send email to
support-netcdf@unidata.ucar.edu. bash-2.04\$ uname -a HP-UX tweety
B.11.00 A 9000/785 2004553471

### Building shared libraries on AIX fails.

On the Unidata AIX platforms, the shared netCDF build fails with either
the Fortran or C++ compilers, like this:

    make  check-TESTS
    *** Creating fills.nc.
    *** SUCCESS!
    PASS: create_fills.sh
    /bin/sh: 16012 Segmentation fault(coredump)
    FAIL: nf_test
    /bin/sh: 16018 Segmentation fault(coredump)
    FAIL: tst_f77_v2
    /bin/sh: 16024 Segmentation fault(coredump)
    FAIL: ftest
    =========================================
    3 of 4 tests failed
    Please report to support@unidata.ucar.edu
    =========================================

If built without fortran or C++, the build succeeds, but shared
libraries are not built. (Static libraries are).

I don't know what is causing this. If any AIX users can shed any light,
that would be most helpful.

Shared builds also fail the same way when using GNU compilers.

### Building with older versions of g++ fails.

The build fails like this:

    libtool: compile:  g++ -DHAVE_CONFIG_H -I. -I. -I.. -I../fortran -DDEBUG -I../libsrc -g -O2 -MT netcdf.lo -MD -MP -MF .deps/netcdf.Tpo -c netcdf.cpp -o netcdf.o
    In file included from /opt/gnu/gcc/include/c++/3.2/powerpc-ibm-aix5.0.0.0/bits/c++io.h:35,
                     from /opt/gnu/gcc/include/c++/3.2/bits/fpos.h:44,
                     from /opt/gnu/gcc/include/c++/3.2/iosfwd:46,
                     from /opt/gnu/gcc/include/c++/3.2/ios:44,
                     from /opt/gnu/gcc/include/c++/3.2/ostream:45,
                     from /opt/gnu/gcc/include/c++/3.2/iostream:45,
                     from netcdf.cpp:13:
    /opt/gnu/gcc/include/c++/3.2/cstdio:108: `fgetpos' not declared
    /opt/gnu/gcc/include/c++/3.2/cstdio:110: `fopen' not declared
    /opt/gnu/gcc/include/c++/3.2/cstdio:115: `freopen' not declared
    /opt/gnu/gcc/include/c++/3.2/cstdio:118: `fsetpos' not declared
    netcdf.cpp: In member function `NcBool NcVar::set_cur(long int, long int, long 
       int, long int, long int)':

This happens in old versions of g++ when large files are used. To fix
this, either upgrade your g++ compiler, or else use --disable-largefile
with configure, to turn off large file handling.

### The .NET build files are not included in the 3.6.2 release.

The netCDF 3.6.2 release does not contain the .NET build files. Whoops!
Sorry about that.

.NET users should use the latest snapshot, or continue to use the 3.6.1
release.

This is now fixed in the netCDF-3 snapshot. Get the snapshot from the
[netCDF-3 snapshot build page](../builds/snapshot).

### Snapshot .NET build files do not work for Visual Studio 8.0 beta.

A user has reported that Visual Studio .NET version 8.0 beta does not
build with the netCDF .NET build files in win32/NET.

Interested users may also wish to subscribe to the [netcdf-porting
mailing
list](http://www.unidata.ucar.edu/mailing_lists/archives/netcdf-porting/).

### The -disable-v2 option causes the fortran build to fail with some fortran compilers.

The netCDF version 2 API is maintained for backward compatibility. We
are committed to maintaining the V2 API for all future releases of
netCDF. However, the --disable-v2 option is provided -- for users who
wish to compile the netCDF libraries without the V2 API.

The --disable-v2 option will cause the fortran build to fail on some
fortran 95 compilers because the netcdf.inc file still includes the
names of the V2 API functions. (Other fortran 90 compilers ignore
these).

If your compiler fails with --disable-v2, you can either refrain from
using this option (that is, build the v2 API as well as the V3 API), or
you can get the netCDF-3 [snapshot](../builds/snapshot).

This is fixed for future releases of netCDF.

### The --disable-c option does not work.

The --disable-c option should turn off the building of the netCDF C
library for use with --enable-separate-fortran (to save a small amount
of tme building and testing. However this does not happen. The C library
is built in any case.

Users may ignore this option in version 3.6.2. It is fixed in the
netCDF-3 [snapshot](../builds/snapshot) and for future releases of
netCDF.

Known Problems with netCDF 3.6.1
--------------------------------

[Building on IBM Bluegene login node (SUSE Linux)](#login_node_3_6_1)

[Linux x86 Fedora4 with Intel ifort 9.0 compiler](#ifort_3_6_1)

### Building on IBM Bluegene login node (SUSE Linux)

Michael McCracken reports the following:

     Hi, I have been working with netcdf and parallel netcdf on
    bluegene at SDSC.  I have no problems building a version to link with
    WRF and run on the ppc32 compute nodes.

    When I need to inspect a data file, I can't run the ncdump that is
    built with the cross-compiling configure without running it on a
    compute node. It works, but I end up waiting in the queue (and being
    charged for using 64 CPUs).

    I'd like to build an ncdump that works on the login node, and I was
    wondering if anyone had done that on a bluegene yet. The local staff
    at SDSC haven't.

    Here's how I did it:

    After some debugging, I can configure and compile with the following
    commands.

    I added -DIBMR2Fortran because apparently that's not getting set on
    bluegene (probably because it's linux and not AIX), and without it
    compiling in the fortran/ subdir barfs.

    setenv CC "xlc"
    setenv CFLAGS "-O3 -qstrict -DIBMR2Fortran"
    setenv CPP "xlc -E"
    setenv CXX "xlC"
    setenv CXXFLAGS "-O3 -qstrict"
    setenv CXXCPP "xlC -E"
    setenv F77 "xlf"
    setenv FC "xlf"
    setenv F90 "xlf90"
    setenv FFLAGS "-O3 -qstrict"

    ./configure --disable-flag-setting

    I get a clean build, and ncdump works for me...

### Linux x86 Fedora4 with Intel ifort 9.0 compiler

For netCDF version 3.6.1, Jonathan Rougier contributes the following
work around for an intel fortran compiler bug.

There is a bug (which may have been fixed: see the ifort pages), which
generates an error message along the lines of:

    IPO link: can not find "("
    ifort: error: problem during multi-file optimization compilation (code
    1)

The documented cludge, which worked for me, is:

    # normal flags

    setenv FC ifort
    setenv FFLAGS  "-mp -recursive"
    setenv CPPFLAGS "-DNDEBUG -DpgiFortran"

    # cludge

    echo null > \(; echo null > AS_NEEDED
    echo null > nf_test/\(; echo null > nf_test/AS_NEEDED
    echo null > ncgen/\(; echo null > ncgen/AS_NEEDED

which creates files called '(' and AS\_NEEDED in the appropriate places.

------------------------------------------------------------------------

\

Known Problems with netCDF 3.6.0
--------------------------------

-   [nctest fails on IRIX platform](#irix-nctest)
-   [C++ API doesn't build on Irix](#irix-CXX-build)
-   [Potentially serious bug with 64-bit offset files](#cdf2-bug)
-   [Cygwin build doesn't work](#bad-cygwin)
-   [Windows DLL doesn't include F77 API](#dll-fortran)
-   [F90 tests fail with Portland F90 compiler](#portland-f90)
-   [Config doesn't find working F77 or F90 compiler on
    AIX](#aix-config)
-   [F90 functions not added to library on AIX](#aix-make)
-   [Problems with fortran compile because of -Df2cFortran being added
    by configure](#fortran-config)
-   [Message: "ncgenyy.c is out-of-date with respect to
    ncgen.l"](#ncgen-timestamp)
-   [Configure help specifies FCFLAGS instead of FFLAGS](#fcflags)
-   [Specifying an edge-length of zero returns error instead of no
    data](#zeroedge)
-   [C++ library doesn't build under Cygwin](#cygwincpp)
-   [Large file problems in Visual C++ compile](#visualcpp_largefile)
-   [When using TEMP_LARGE, need a trailing slash](#temp_large)

------------------------------------------------------------------------

### nctest fails on IRIX platform

It has been reported (by Atro Tossavainen) that nctest fails on some
Irix builds. (We cannot duplicate this problem at netCDF World HQ).

The nctest fails when comparing the test-generated out with a saved copy
of the output.

This problem was fixed in the 3.6.1 release.

### C++ API doesn't build on Irix

On Irix systems without a recent version of the C++ compiler, the C++
API won't build. The solution is to set CXXFLAGS to -LANG:std.

### Potentially serious bug with 64-bit offset files

Kevin Thomas of the University of Oklahoma has reported a potentially
serious bug in using the new large file support in netCDF 3.6.0. Users
of the new large file facilities are cautioned to either apply [this
one-line patch to netCDF
3.6.0](/software/netcdf/patches/patch-3.6.0-cdf2) or to upgrade from
version 3.6.0 to the current release version 3.6.0-p1, available from
[netcdf.tar.Z](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf.tar.Z) or
[netcdf.tar.gz](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf.tar.gz).
Until you can upgrade, avoid rewriting in place any large (&gt; 2 GiB)
netCDF files that use the new 64-bit offset format under the conditions
described below.

The bug occurs following this sequence of steps:

A large netCDF file is created using the new 64-bit offset format
variant (also known as the version 2 format) and the file is closed.

Later the file is opened for writing, followed by either of the
following operations:

-   enter define mode (calling nc\_redef(), nf\_redef(), or
    nf90\_redef(), in C, Fortran-77, or Fortran-90 interface,
    for example) to add a new dimension, variable, or attribute; or
-   write a new value for an existing attribute (either a global or a
    variable-specific attribute).

Under these conditions, after you leave define mode or close the file,
the file header is written out in the "classic" (version 1) netCDF
format, chopping the leading bits off any variable offsets that were
large enough to require more than 32 bits. If there were no such huge
variable offsets, the file is undamaged and remains readable as a
classic netCDF file. If there were any huge variable offsets (&gt; 2
GiB), data for the first such variable and all subsequent variables will
not be accessed correctly. It is possible to restore the header for such
a file to the correct 64-bit offset form so that the data can
subsequently be accessed correctly, if no data values have been
overwritten since the file header was changed to classic format. Feel
free to [contact us](mailto:support-netcdf@unidata.ucar.edu) for help
restoring the file headers if this applies to you. If you have any large
64-bit offset format netCDF files that might have mistakenly been
rewritten with classic format headers, please be careful not to write
any more data into them, as it could overwrite data that could not
subsequently be recovered.

If you want to know how to tell if a 64-bit offset file has been
converted by this bug into a classic format file, see the answer to the
FAQ [How can I tell if a netCDF file uses the classic format or new
64-bit offset
format?](/software/netcdf/faq.html#Large%20File%20Support5).

### Cygwin Build Doesn't Work

To build on Cygwin, you must get the [latest 3.6.1 beta
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-beta.tar.gz).

### Windows DLL doesn't include F77 API

The netCDF windows DLL doesn't include the Fortran API. We are working
on this problem for the next release. Meanwhile, if you need the fortran
API in your DLL, you'll have to use the [netCDF 3.5.1
DLL](ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/netcdf-3.5.1-win32dll.zip).

### F90 tests fail with Portland F90 compiler

On some versions of the Portland Group F90 compiler, the F90 tests fail,
looking something like this:


    *** Failure ***
    *** example_good.cdl    2000-04-05 21:33:14.000000000 +0200
    --- example.cdl 2005-01-11 10:21:31.624884000 +0100
    ***************
    *** 34,43 ****
        953, 954, 955,
        956, 957, 958,
        959, 960, 961,
    !   962, 963, 964,
    !   965, 966, 967,
    !   968, 969, 970,
    !   971, 972, 973 ;

       lat = -90, -87.5, -85, -82.5 ;

    --- 34,43 ----
        953, 954, 955,
        956, 957, 958,
        959, 960, 961,
    !   950, 951, 952,
    !   953, 954, 955,
    !   956, 957, 958,
    !   959, 960, 961 ;

This problem is caused by a bug in the Portland F90 compiler. Upgrade to
the latest version of the compiler or get the free patch from Portland
Group to fix this.

### Config doesn't find working F77 or F90 compiler on AIX

On AIX systems, the configure step can't find either the F90 or the F77
compiler. On AIX system, you must set the environment variables FC and
F90 to xlf and xlf95.

### xlf90 fails to compile test program during configure on AIX

On AIX systems, the F90 option -qsuffix=f=f90 is required in F90FLAGS.
Configure should automatically detect and add this to F90FLAGS if it's
not already there, but it doesn't.

FIX: Make sure that -qsuffix=f=f90 is set in the F90FLAGS before running
configure.

This will be fixed in the next beta release.

### F90 functions not added to library on AIX

On AIX systems, the F90 functions may not be added to the library. This
is due to a quirk of AIX make.

Before doing "make install", change to the Fortran90 directory (src/f90)
and do "make". Then proceed to "make install".

### Problems with fortran compile because of -Df2cFortran being added by configure"

With some fortran compilers, such as Absoft, the configure script
stupidly adds a -Df2cFortran to the C preprocessor flags, which causes
the fortran tests in nf\_test to fail to link.

This problem is fixed in the 3.6.1 beta release. Get the [3.6.1 beta
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-beta.tar.Z).

### Message: "ncgenyy.c is out-of-date with respect to ncgen.l"

On some platforms (HP-UX 11.00, maybe others), make fails with an error
message like:

      Warning: ncgenyy.c is out-of-date with respect to ncgen.l
      Warning: It should be recreated via flex on a SunOS 5 system

and then fails if the "flex" command is not found.

The problem is that the modification time on the source file
src/ncgen/ncgenyy.c is being interpreted as earlier than the
modification time on the source file src/ncgen/ncgen.l, even though
ncgenyy.c was actually created after ncgen.l was modified. To workaround
this problem on a Unix system, run the following command from the netCDF
src/ directory to update the modification time of the derived file:

      touch ncgen/ncgenyy.c

Then rerun the make command.

### Configure help specifies FCFLAGS instead of FFLAGS

If you run "configure --help", it suggests setting "FCFLAGS" for the
fortran compiler flags, but "FFLAGS" is actually used for the Fortran
compiler flags. "FCFLAGS" is ignored when compiling.

This problem will be is fixed in the next beta release. Until then, use
FFLAGS, not FCFLAGS.

### Specifying a count length of zero returns an error instead of no data

For access to array sections, strided access, or mapped access, you need
to specify both a start index vector and a count vector, where the count
vector specifies the number of slices along each edge to access. If the
start index vector specifies the maximum dimension size and the
corresponding count vector is zero, the library should just return no
data, but instead it returns an error status indicating "Index exceeds
dimension bound". This problem has been present in all versions of
netCDF, and the test programs even verify that in this case an error is
returned rather than gracefully accessing no data.

This will be fixed in the next minor version.

### C++ library doesn't build under Cygwin

Running configure on Cygwin fails to find GNU C++ compiler, even if it
is present on the platform. As a result, the C++ interface is never
built.

This problem is fixed in the 3.6.1 beta release. Cygwin users interested
in the C++ interface should get the [3.6.1 beta
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-3.6.1-beta1.tar.Z).

### Large file problems in Visual C++ compile

The use of large files, and an 8-byte off\_t type, is not handled
correctly in the 3.6.0 release of the code and project files needed to
compile the netCDF library with Visual C++.NET.

This problem is fixed in the 3.6.1 beta release. Users interested in
building their own DLL should get the [3.6.1 beta
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-3.6.1-beta1.tar.Z).
The DLL offered on the binary release page is 3.6.1 beta.

### When using TEMP\_LARGE, need a trailing slash

When using the environment variable TEMP\_LARGE during the netCDF 3.6.0
make extra\_test phase, the directory name must be followed by a slash
to work. For example, use 'setenv TEMP\_LARGE /tmp/' instead of 'setenv
TEMP\_LARGE /tmp', as one would usually expect, and as the documentation
describes.

This problem is fixed in the [3.6.1 beta
release](ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-3.6.1-beta1.tar.Z).
Users of 3.6.0 should specify the trailing slash to use the TEMP\_LARGE
environment variable in make extra\_test.

------------------------------------------------------------------------

Reported problems and workarounds are also available for some previous
releases: [version 3.5](/software/netcdf/known_problems_35.html),
[version 3.4](/software/netcdf/known_problems_34.html), [version
3.3.1](/software/netcdf/known_problems_331.html), [version
3.3](/software/netcdf/known_problems_330.html), [version
2.4.3](/software/netcdf/known_problems_243.html), and [version
2.4.2](/software/netcdf/known_problems_242.html).
