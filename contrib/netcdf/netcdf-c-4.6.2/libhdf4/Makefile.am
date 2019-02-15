# This is part of Unidata's netCDF package. Copyright 2018, see the
# COPYRIGHT file for more information.

# Build the HDF4 dispatch layer.

# Ed Hartnett

# Get AM_CPPFLAGS.
include $(top_srcdir)/lib_flags.am
libnetcdf4_la_CPPFLAGS = ${AM_CPPFLAGS}

# This is our output. The HDF4 convenience library.
noinst_LTLIBRARIES = libnchdf4.la

# The source files.
libnchdf4_la_SOURCES = hdf4dispatch.c hdf4file.c hdf4var.c hdf4func.c

# Package this for cmake build.
EXTRA_DIST = CMakeLists.txt





