NetCDF-4 Filter Support
============================
<!-- double header is needed to workaround doxygen bug -->

NetCDF-4 Filter Support {#compress}
=================================

[TOC]

Introduction {#compress_intro}
==================

The HDF5 library (1.8.11 and later) 
supports a general filter mechanism to apply various
kinds of filters to datasets before reading or writing. 
The netCDF enhanced (aka netCDF-4) library inherits this
capability since it depends on the HDF5 library.

Filters assume that a variable has chunking
defined and each chunk is filtered before
writing and "unfiltered" after reading and
before passing the data to the user.

The most common kind of filter is a compression-decompression
filter, and that is the focus of this document.

HDF5 supports dynamic loading of compression filters using the following
process for reading of compressed data.

1. Assume that we have a dataset with one or more variables that
were compressed using some algorithm. How the dataset was compressed
will be discussed subsequently.

2. Shared libraries or DLLs exist that implement the compress/decompress
algorithm. These libraries have a specific API so that the HDF5 library
can locate, load, and utilize the compressor.
These libraries are expected to installed in a specific
directory.

Enabling A Compression Filter {#Enable}
=============================

In order to compress a variable, the netcdf-c library
must be given three pieces of information:
(1) some unique identifier for the filter to be used,
(2) a vector of parameters for
controlling the action of the compression filter, and
(3) a shared library implementation of the filter.

The meaning of the parameters is, of course,
completely filter dependent and the filter
description [3] needs to be consulted. For
bzip2, for example, a single parameter is provided
representing the compression level.
It is legal to provide a zero-length set of parameters.
Defaults are not provided, so this assumes that
the filter can operate with zero parameters.

Filter ids are assigned by the HDF group. See [4]
for a current list of assigned filter ids.
Note that ids above 32767 can be used for testing without
registration.

The first two pieces of  information can be provided in one of three ways:
using __ncgen__, via an API call, or via command line parameters to __nccopy__.
In any case, remember that filtering also requires setting chunking, so the
variable must also be marked with chunking information.

Using The API {#API}
-------------
The necessary API methods are included in __netcdf.h__ by default.
One API method is for setting the filter to be used
when writing a variable. The relevant signature is
as follows.
````
int nc_def_var_filter(int ncid, int varid, unsigned int id, size_t nparams, const unsigned int* parms);
````
This must be invoked after the variable has been created and before
__nc_enddef__ is invoked.

A second API methods makes it possible to query a variable to
obtain information about any associated filter using this signature.
````
int nc_inq_var_filter(int ncid, int varid, unsigned int* idp, size_t* nparams, unsigned int* params);

````
The filter id will be returned in the __idp__ argument (if non-NULL),
the number of parameters in __nparamsp__ and the actual parameters in
__params__.  As is usual with the netcdf API, one is expected to call
this function twice. The first time to get __nparams__ and the
second to get the parameters in client-allocated memory.

Using ncgen {#NCGEN}
-------------

In a CDL file, compression of a variable can be specified
by annotating it with the following attribute:

* ''_Filter'' -- a string containing a comma separated list of
constants specifying (1) the filter id to apply, and (2)
a vector of constants representing the
parameters for controlling the operation of the specified filter.
See the section on the <a href="#Syntax">parameter encoding syntax</a>
for the details on the allowable kinds of constants.

This is a "special" attribute, which means that
it will normally be invisible when using
__ncdump__ unless the -s flag is specified.

Example CDL File (Data elided)
------------------------------
````
netcdf bzip2 {
dimensions:
  dim0 = 4 ; dim1 = 4 ; dim2 = 4 ; dim3 = 4 ;
variables:
  float var(dim0, dim1, dim2, dim3) ;
    var:_Filter = "307,9" ;
    var:_Storage = "chunked" ;
    var:_ChunkSizes = 4, 4, 4, 4 ;
data:
...
}
````

Using nccopy {#NCCOPY}
-------------
When copying a netcdf file using __nccopy__ it is possible
to specify filter information for any output variable by
using the "-F" option on the command line; for example:
````
nccopy -F "var,307,9" unfiltered.nc filtered.nc
````
Assume that __unfiltered.nc__ has a chunked but not bzip2 compressed
variable named "var". This command will create that variable in
the __filtered.nc__ output file but using filter with id 307
(i.e. bzip2) and with parameter(s) 9 indicating the compression level.
See the section on the <a href="#Syntax">parameter encoding syntax</a>
for the details on the allowable kinds of constants.

The "-F" option can be used repeatedly as long as the variable name
part is different. A different filter id and parameters can be
specified for each occurrence.

As a rule, any input filter on an input variable will be applied
to the equivalent output variable -- assuming the output file type
is netcdf-4. It is, however, sometimes convenient to suppress
output compression either totally or on a per-variable basis.
Total suppression of output filters can be accomplished by specifying
a special case of "-F", namely this.
````
nccopy -F "none" input.nc output.nc
````
Suppression of output filtering for a specific variable can be accomplished
using this format.
````
nccopy -F "var,none" input.nc output.nc
````
where "var" is the fully qualified name of the variable.

The rules for all possible cases of the "-F" flag are defined
by this table.

<table>
<tr><th>-F none<th>-Fvar,...<th>Input Filter<th>Applied Output Filter
<tr><td>true<td>unspecified<td>NA<td>unfiltered
<tr><td>true<td>-Fvar,none<td>NA<td>unfiltered
<tr><td>true<td>-Fvar,...<td>NA<td>use output filter
<tr><td>false<td>unspecified<td>defined<td>use input filter
<tr><td>false<td>-Fvar,none<td>NA<td>unfiltered
<tr><td>false<td>-Fvar,...<td>NA<td>use output filter
</table> 

Parameter Encoding {#ParamEncode}
==========

The parameters passed to a filter are encoded internally as a vector
of 32-bit unsigned integers. It may be that the parameters
required by a filter can naturally be encoded as unsigned integers.
The bzip2 compression filter, for example, expects a single
integer value from zero thru nine. This encodes naturally as a
single unsigned integer.

Note that signed integers and single-precision (32-bit) float values
also can easily be represented as 32 bit unsigned integers by
proper casting to an unsigned integer so that the bit pattern
is preserved. Simple integer values of type short or char
(or the unsigned versions) can also be mapped to an unsigned
integer by truncating to 16 or 8 bits respectively and then
zero extending.

Machine byte order (aka endian-ness) is an issue for passing
some kinds of parameters. You might define the parameters when
compressing on a little endian machine, but later do the
decompression on a big endian machine. Byte order is not an
issue for 32-bit values because HDF5 takes care of converting
them between the local machine byte order and network byte
order.

Parameters whose size is larger than 32-bits present a byte order problem.
This typically includes double precision floats and (signed or unsigned)
64-bit integers. For these cases, the machine byte order must be
handled by the compression code. This is because HDF5 will treat,
for example, an unsigned long long as two 32-bit unsigned integers
and will convert each to network order separately. This means that
on a machine whose byte order is different than the machine in which
the parameters were initially created, the two integers are out of order
and must be swapped to get the correct unsigned long long value.
Consider this example. Suppose we have this little endian unsigned long long.

    1000000230000004

In network byte order, it will be stored as two 32-bit integers.

    20000001 40000003

On a big endian machine, this will be given to the filter in that form.

    2000000140000003

But note that the proper big endian unsigned long long form is this.

4000000320000001

So, the two words need to be swapped.

But consider the case when both original and final machines are big endian.

1. 4000000320000001
2. 40000003 20000001
3. 40000003 20000001

where #1 is the original number, #2 is the network order and
#3 is the what is given to the filter. In this case we do not
want to swap words.

The solution is to forcibly encode the original number using some
specified endianness so that the filter always assumes it is getting
its parameters in that order and will always do swapping as needed.
This is irritating, but one needs to be aware of it. Since most
machines are little-endian. We choose to use that as the endianness
for handling 64 bit entities.

Filter Specification Syntax {#Syntax}
==========

Both of the utilities
<a href="#NCGEN">__ncgen__</a>
and
<a href="#NCCOPY">__nccopy__</a>
allow the specification of filter parameters.
These specifications consist of a sequence of comma
separated constants. The constants are converted
within the utility to a proper set of unsigned int
constants (see the <a href="#ParamEncode">parameter encoding section</a>).

To simplify things, various kinds of constants can be specified
rather than just simple unsigned integers. The utilities will encode
them properly using the rules specified in 
the <a href="#ParamEncode">parameter encoding section</a>.

The currently supported constants are as follows.
<table>
<tr halign="center"><th>Example<th>Type<th>Format Tag<th>Notes
<tr><td>-17b<td>signed 8-bit byte<td>b|B<td>Truncated to 8 bits and zero extended to 32 bits
<tr><td>23ub<td>unsigned 8-bit byte<td>u|U b|B<td>Truncated to 8 bits and zero extended to 32 bits
<tr><td>-25S<td>signed 16-bit short<td>s|S<td>Truncated to 16 bits and zero extended to 32 bits
<tr><td>27US<td>unsigned 16-bit short<td>u|U s|S<td>Truncated to 16 bits and zero extended to 32 bits
<tr><td>-77<td>implicit signed 32-bit integer<td>Leading minus sign and no tag<td>
<tr><td>77<td>implicit unsigned 32-bit integer<td>No tag<td>
<tr><td>93U<td>explicit unsigned 32-bit integer<td>u|U<td>
<tr><td>789f<td>32-bit float<td>f|F<td>
<tr><td>12345678.12345678d<td>64-bit double<td>d|D<td>Network byte order
<tr><td>-9223372036854775807L<td>64-bit signed long long<td>l|L<td>Network byte order
<tr><td>18446744073709551615UL<td>64-bit unsigned long long<td>u|U l|L<td>Network byte order
</table>
Some things to note.

1. In all cases, except for an untagged positive integer,
   the format tag is required and determines how the constant
   is converted to one or two unsigned int values.
   The positive integer case is for backward compatibility.
2. For signed byte and short, the value is sign extended to 32 bits
   and then treated as an unsigned int value.
3. For double, and signed|unsigned long long, they are converted
   to network byte order and then treated as two unsigned int values.
   This is consistent with the <a href="#ParamEncode">parameter encoding</a>.

Dynamic Loading Process {#Process}
==========

The documentation[1,2] for the HDF5 dynamic loading was (at the time
this was written) out-of-date with respect to the actual HDF5 code
(see HDF5PL.c). So, the following discussion is largely derived
from looking at the actual code. This means that it is subject to change.

Plugin directory {#Plugindir}
----------------

The HDF5 loader expects plugins to be in a specified plugin directory.
The default directory is:
   * "/usr/local/hdf5/lib/plugin” for linux/unix operating systems (including Cygwin)
   * “%ALLUSERSPROFILE%\\hdf5\\lib\\plugin” for Windows systems, although the code
     does not appear to explicitly use this path.

The default may be overridden using the environment variable
__HDF5_PLUGIN_PATH__.

Plugin Library Naming {#Pluginlib}
---------------------

Given a plugin directory, HDF5 examines every file in that
directory that conforms to a specified name pattern
as determined by the platform on which the library is being executed.
<table>
<tr halign="center"><th>Platform<th>Basename<th>Extension
<tr halign="left"><td>Linux<td>lib*<td>.so*
<tr halign="left"><td>OSX<td>lib*<td>.so*
<tr halign="left"><td>Cygwin<td>cyg*<td>.dll*
<tr halign="left"><td>Windows<td>*<td>.dll
</table>

Plugin Verification {#Pluginverify}
-------------------
For each dynamic library located using the previous patterns,
HDF5 attempts to load the library and attempts to obtain information
from it. Specifically, It looks for two functions with the following
signatures.

1. __H5PL_type_t H5PLget_plugin_type(void)__ --
This function is expected to return the constant value
__H5PL_TYPE_FILTER__ to indicate that this is a filter library.
2. __const void* H5PLget_plugin_info(void)__ --
This function returns a pointer to a table of type __H5Z_class2_t__.
This table contains the necessary information needed to utilize the
filter both for reading and for writing. In particular, it specifies
the filter id implemented by the library and if must match that id
specified for the variable in __nc_def_var_filter__ in order to be used.

If plugin verification fails, then that plugin is ignored and
the search continues for another, matching plugin.

Debugging {#Debug}
-------
Debugging plugins can be very difficult. You will probably
need to use the old printf approach for debugging the filter itself.

One case worth mentioning is when you have a dataset that is
using an unknown filter. For this situation, you need to
identify what filter(s) are used in the dataset. This can
be accomplished using this command.
````
ncdump -s -h <dataset filename>
````
Since ncdump is not being asked to access the data (the -h flag), it
can obtain the filter information without failures. Then it can print
out the filter id and the parameters (the -s flag).

Test Case {#TestCase}
-------
Within the netcdf-c source tree, the directory
__netcdf-c/nc_test4__ contains a test case (__test_filter.c__) for
testing dynamic filter writing and reading using
bzip2. Another test (__test_filter_misc.c__) validates
parameter passing.  These tests are disabled if __--enable-shared__
is not set or if __--enable-netcdf-4__ is not set.

Example {#Example}
-------
A slightly simplified version of the filter test case is also
available as an example within the netcdf-c source tree
directory __netcdf-c/examples/C. The test is called __filter_example.c__
and it is executed as part of the __run_examples4.sh__ shell script.
The test case demonstrates dynamic filter writing and reading.

The files __example/C/hdf5plugins/Makefile.am__
and  __example/C/hdf5plugins/CMakeLists.txt__
demonstrate how to build the hdf5 plugin for bzip2.

Notes
==========

Memory Allocation Issues
-----------

Starting with HDF5 version 1.10.x, the plugin code MUST be
careful when using the standard *malloc()*, *realloc()*, and
*free()* function.

In the event that the code is allocating, reallocating, for
free'ing memory that either came from or will be exported to the
calling HDF5 library, then one MUST use the corresponding HDF5
functions *H5allocate_memory()*, *H5resize_memory()*,
*H5free_memory()* [5] to avoid memory failures.

Additionally, if your filter code leaks memory, then the HDF5 library
generates a failure something like this.
````
H5MM.c:232: H5MM_final_sanity_check: Assertion `0 == H5MM_curr_alloc_bytes_s' failed.
````

One can look at the the code in plugins/H5Zbzip2.c and H5Zmisc.c to see this.

SZIP Issues
-----------
The current szip plugin code in the HDF5 library
has some behaviors that can catch the unwary.
Specifically, this filter may do two things.

1. Add extra parameters to the filter parameters: going from
   the two parameters provided by the user to four parameters
   for internal use. It turns out that the two parameters provided
   when calling nc_def_var_filter correspond to the first two
   parameters of the four parameters returned by nc_inq_var_filter.
2. Change the values of some parameters: the value of the
   __options_mask__ argument is known to add additional flag bits,
   and the __pixels_per_block__ parameter may be modified.

The reason for these changes is has to do with the fact that
the szip API provided by the underlying H5Pset_szip function
is actually a subset of the capabilities of the real szip implementation.
Presumably this is for historical reasons.

In any case, if the caller uses the __nc_inq_var_szip__, then
the values returned may differ from those originally specified.
If one used the __nc_inq_var_filter__ API calls, it may be the case that
both the number of parameters and the values will differ from the original
call to __nc_def_var_filter__.

Supported Systems
-----------------
The current matrix of OS X build systems known to work is as follows.
<table>
<tr><th>Build System<th>Supported OS
<tr><td>Automake<td>Linux, Cygwin
<tr><td>Cmake<td>Linux, Cygwin, Visual Studio
</table>

Generic Plugin Build
--------------------
If you do not want to use Automake or Cmake, the following
has been known to work.
````
gcc -g -O0 -shared -o libbzip2.so <plugin source files>  -L${HDF5LIBDIR} -lhdf5_hl -lhdf5 -L${ZLIBDIR} -lz
````

Appendix A. Byte Swap Code {#AppendixA}
==========
Since in some cases, it is necessary for a filter to
byte swap from little-endian to big-endian, This appendix
provides sample code for doing this. It also provides
a code snippet for testing if the machine the
endianness of a machine.

Byte swap an 8-byte chunk of memory
-------
````
static void
byteswap8(unsigned char* mem)
{
    register unsigned char c;
    c = mem[0];
    mem[0] = mem[7];
    mem[7] = c;
    c = mem[1];
    mem[1] = mem[6];
    mem[6] = c;
    c = mem[2];
    mem[2] = mem[5];
    mem[5] = c;
    c = mem[3];
    mem[3] = mem[4];
    mem[4] = c;
}

````

Test for Machine Endianness
-------
````
static const unsigned char b[4] = {0x0,0x0,0x0,0x1}; /* value 1 in big-endian*/
int endianness = (1 == *(unsigned int*)b); /* 1=>big 0=>little endian
````
References {#References}
========================

1. https://support.hdfgroup.org/HDF5/doc/Advanced/DynamicallyLoadedFilters/HDF5DynamicallyLoadedFilters.pdf
2. https://support.hdfgroup.org/HDF5/doc/TechNotes/TechNote-HDF5-CompressionTroubleshooting.pdf
3. https://portal.hdfgroup.org/display/support/Contributions#Contributions-filters
4. https://support.hdfgroup.org/services/contributions.html#filters
5. https://support.hdfgroup.org/HDF5/doc/RM/RM_H5.html

Point of Contact
================

__Author__: Dennis Heimbigner<br>
__Email__: dmh at ucar dot edu
__Initial Version__: 1/10/2018<br>
__Last Revised__: 2/5/2018

