## This is an automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.

# This automake file is in charge of building the libsrc directory,
# which contains the classic library code.

include $(top_srcdir)/lib_flags.am

libnetcdfp_la_CPPFLAGS = ${AM_CPPFLAGS}

# This is the code for a dispatch table for pnetcdf
# (which has CDF5 as its magic number)

# Turn on a pre-processor flag when building a DLL for windows.
if BUILD_DLL
libnetcdfp_la_CPPFLAGS += -DDLL_EXPORT
endif # BUILD_DLL

# These files comprise the pnetcdf dispatch library code.
libnetcdfp_la_SOURCES = ncpdispatch.c

noinst_LTLIBRARIES = libnetcdfp.la

EXTRA_DIST = CMakeLists.txt 
