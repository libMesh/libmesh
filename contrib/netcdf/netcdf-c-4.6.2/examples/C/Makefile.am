# This is a automake file, part of Unidata's netCDF package.
# Copyright 2006, see the COPYRIGHT file for more information.

# This file builds the C examples.

# Ed Hartnett, Ward Fisher, Dennis Heimbigner

# Un comment to use a more verbose test driver
#SH_LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose

LDADD = -lm
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS =

# These are the netCDF-3 examples.
check_PROGRAMS = simple_xy_wr simple_xy_rd sfc_pres_temp_wr	\
sfc_pres_temp_rd pres_temp_4D_wr pres_temp_4D_rd

TESTS = run_examples.sh

# To build netcdf-4, or not to build netcdf-4, that is the question...
if USE_NETCDF4
# These are the extra netCDF-4 examples.
check_PROGRAMS += simple_nc4_wr simple_nc4_rd simple_xy_nc4_wr	\
simple_xy_nc4_rd

if ENABLE_FILTER_TESTING
# filter_example.c should be same as nc_test4/test_filter.c
check_PROGRAMS += filter_example
TESTS += run_filter.sh
endif

endif #USE_NETCDF4

if USE_PNETCDF
# These are the extra netCDF-4 examples.
check_PROGRAMS += parallel_vara
TESTS += parallel_vara
endif #USE_PNETCDF

AM_CPPFLAGS += -I$(top_builddir)/liblib
AM_LDFLAGS += ${top_builddir}/liblib/libnetcdf.la

if USE_NETCDF4
TESTS += run_examples4.sh
endif #USE_NETCDF4

# These files are created by the tests.
CLEANFILES = *.nc

EXTRA_DIST = CMakeLists.txt run_examples.sh run_examples4.sh run_filter.sh

DISTCLEANFILES = findplugin.sh
