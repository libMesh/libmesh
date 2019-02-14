# This is part of Unidata's netCDF package.  Copyright 2005-2011, see
# the COPYRIGHT file for more information.

# This automake file generates the Makefile to build netCDF-4.
# Ed Hartnett

include $(top_srcdir)/lib_flags.am

libnetcdf4_la_CPPFLAGS = ${AM_CPPFLAGS}

# This is our output. The netCDF-4 convenience library.
noinst_LTLIBRARIES = libnetcdf4.la
libnetcdf4_la_SOURCES = nc4dispatch.c nc4attr.c nc4dim.c nc4grp.c	\
nc4internal.c nc4type.c nc4var.c ncfunc.c error4.c nc4printer.c		\
ncindex.c

EXTRA_DIST = CMakeLists.txt





