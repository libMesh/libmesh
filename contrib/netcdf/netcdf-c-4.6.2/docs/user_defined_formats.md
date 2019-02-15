User-Defined Formats for NetCDF {#user_defined_formats}
===============================

[TOC]

User-Defined Formats {#udf_user_defined_formats}
=====================================

## Introduction {#udf_Introduction}

User-defined formats allow users to write their own adaptors for the
netCDF C library, so that it can read and (optionally) write a
proprietary format through the netCDF API.

This capability is currently experimental. It involves the exposing of internal
netcdf interfaces and data structures that were previously invisible to users.
This means that it is unstable and the exposed interfaces are subject to change.
Use with caution.

User-defined format code is packaged into a separate library, the
user-defined format dispatch library. This library, when linked with
the netCDF library, will allow user programs to read their proprietary
format through the netCDF API. The proprietary format is treated as if
it were one of the netCDF C library native binary formats.

Coding the user-defined format dispatch library requires knowledge of
the netCDF library internals. User-defined format dispatch libraries
must be written in C.

### Magic Numbers

Some file formats use the first few bytes of the file as an identifier
for format. For example, HDF5 files have "HDF5" as the fist 4 bytes,
and netCDF classic files have "CDF1" as the first four bytes. This is
called the "magic number" of the file.

User-defined formats can optionally support magic numbers. If the
user-defined format uses a magic number, and that magic number is
associated with the user-defined format, then netCDF will be able to
correctly identify those files from nc_open(). It will not be
necessary for the user to know or specify the underlying format.

## Using User-Defined Formats from C Programms {#udf_With_C}

A user-defined format can be added dynamically in the case of C programs.

```
      /* Add our test user defined format. */
      if (nc_def_user_format(NC_UDF0, &tst_dispatcher, NULL)) ERR;
```

The file can now be opened by netCDF:

```
      if (nc_open(FILE_NAME, NC_UDF0, &ncid)) ERR;
```

If a magic number is used in the file, that may be passed to
nc_def_user_format(). In that case, specifying the NC_UDF0 mode flag
to nc_open() is optional. The nc_open() will check the file and find
the magic number, and automatically associate the file with
NC_UDF0. The user will not need to know the format in order to open
the file with nc_open().

## Building NetCDF C Library with a User-Defined Format Library {#udf_Build_NetCDF_With_UDF}

Once a user-defined format library is created, it may built into a
netCDF install. This allows the netCDF Fortran APIs, and the netCDF
utilities (ncdump, ncgen, nccopy) to natively use the user-defined
format.

First the user-defined dispatch library must be built and installed.

Then the netcdf-c package must be (re-)built. When building netcdf-c,
add the location of the user-defined format dispatch library include
file to the CPPFLAGS, and the location of the user-defined format
dispatch library in LDFLAGS.

Configure netcdf-c with the option ````--with-udf0=<udf_lib_name>````.

If a magic number is associated with the user-defined format, it can
be specified with the --with-udf0-magic-number= argument.

## Creating a User-Defined Format {#udf_Create_UDF}

Creators of user-defined format libraries will have to become familar
with the internals of the netCDF-4 code.

### Read-Only User-Defined Formats

Many users will find that a read-only user-defined formats meets most
of their needs. With a read-only user-defined format, netCDF will be
able to read files of the user-defined format. Tools like ncdump and
nccopy can work on the files.

A read-only user-defined format can be implemented with only 6
functions. The code in libhdf4 is an excellent example of a read-only
dispatch layer.

## Examples {#udf_Examples}

The most simple-case example of a user-defined format is provided in
test nc_test4/tst_udf.c.

A slightly more complex example, including the required
autoconf/automake files to build a user-defined format library, can be
found at the [sample user-defined format
library](https://github.com/NOAA-GSD/sample-netcdf-dispatch). In this
example, the HDF4 SD reader is re-implemented as an external
user-defined format. (This is unnecessary if you just want to read
HDF4 SD files, since the netCDF C library already includes an HDF4 SD
file reader. This user-defined format library uses the same code. It
is repackaged as a user-defined library to provide a working sample.)

