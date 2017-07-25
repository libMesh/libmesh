# Unidata NetCDF

[![Build Status](https://travis-ci.org/Unidata/netcdf-c.svg?branch=master)](https://travis-ci.org/Unidata/netcdf-c)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/157/badge.svg)](https://scan.coverity.com/projects/157)

### About
The Unidata network Common Data Form (**netCDF**) is an interface for
scientific data access and a freely-distributed software library that
provides an implementation of the interface.  The netCDF library also
defines a machine-independent format for representing scientific data.
Together, the interface, library, and format support the creation,
access, and sharing of scientific data.  The current netCDF software
provides C interfaces for applications and data.  Separate software
distributions available from Unidata provide Java, Fortran, Python,
and C++ interfaces.  They have been tested on various common
platforms.

#### Properties
NetCDF files are self-describing, network-transparent, directly
accessible, and extendible.  `Self-describing` means that a netCDF file
includes information about the data it contains.  `Network-transparent`
means that a netCDF file is represented in a form that can be accessed
by computers with different ways of storing integers, characters, and
floating-point numbers.  `Direct-access` means that a small subset of a
large dataset may be accessed efficiently, without first reading through
all the preceding data.  `Extendible` means that data can be appended to
a netCDF dataset without copying it or redefining its structure.

#### Use
NetCDF is useful for supporting access to diverse kinds of scientific
data in heterogeneous networking environments and for writing
application software that does not depend on application-specific
formats.  For information about a variety of analysis and display
packages that have been developed to analyze and display data in
netCDF form, see

* [Software for Manipulating or Displaying NetCDF Data](http://www.unidata.ucar.edu/netcdf/software.html)

##### More information
For more information about netCDF, see

* [Unidata Network Common Data Form (NetCDF)](http://www.unidata.ucar.edu/netcdf/)

### Latest releases
You can obtain a copy of the latest released version of netCDF
software for various languages:

* [C library and utilities](http://github.com/Unidata/netcdf-c)
* [Fortran](http://github.com/Unidata/netcdf-fortran)
* [Java](http://www.unidata.ucar.edu/downloads/netcdf/netcdf-java-4/)
* [Python](http://github.com/Unidata/netcdf4-python)
* [C++](http://github.com/Unidata/netcdf-cxx4)

### Copyright
Copyright and licensing information can be found [here](http://www.unidata.ucar.edu/software/netcdf/copyright.html), as well as in the COPYRIGHT file accompanying the software

### Installation
To install the netCDF-C software, please see the file INSTALL in the
netCDF-C distribution, or the (usually more up-to-date) document:

* [Building NetCDF](http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)

### Documentation
A language-independent User's Guide for netCDF, and some other
language-specific user-level documents are available from:

* [Language-independent User's Guide](http://www.unidata.ucar.edu/software/netcdf/docs/user_guide.html)
* [NetCDF-C Tutorial](http://www.unidata.ucar.edu/software/netcdf/docs/tutorial_8dox.html)
* [Fortran-90 User's Guide](http://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf-f90/index.html#Top)
* [Fortran-77 User's Guide](http://www.unidata.ucar.edu/software/netcdf/documentation/historic/netcdf-f77/index.html#Top)
* [netCDF-Java/Common Data Model library](http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/)
* [netCDF4-python](http://unidata.github.io/netcdf4-python/)

A mailing list, netcdfgroup@unidata.ucar.edu, exists for discussion of
the netCDF interface and announcements about netCDF bugs, fixes, and
enhancements.  For information about how to subscribe, see the URL

* [Unidata netCDF Mailing-Lists](http://www.unidata.ucar.edu/netcdf/mailing-lists.html)

### Feedback
We appreciate feedback from users of this package.  Please send comments, suggestions, and bug reports to <support-netcdf@unidata.ucar.edu>.  
