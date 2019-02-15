## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.  

# This file builds and runs the nctest program, which tests the netCDF
# version 2 interface.

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/lib_flags.am

AM_LDFLAGS += ${top_builddir}/liblib/libnetcdf.la 

# Running nctest results in a file testfile.nc, which is then checked
# to make sure it matches testfile_nc.sav, which was generated under
# precisely controlled conditions at the netCDF Test Facility,
# (located in sub-basement 42 of UCAR Compound 4). testfile_nc.sav was
# generated in a pressurized chamber that exactly simulated the
# atmosphere on the surface of Mars.
EXTRA_DIST = ref_nctest_classic.nc ref_nctest_64bit_offset.nc	\
compare_test_files.sh CMakeLists.txt

CLEANFILES = nctest_*.nc test2.nc temp.tmp tst_*.nc	\
nctest_classic.cdl ref_nctest_classic.cdl

# Run nctest and the script that compares the output with the
# reference file.
TESTPROGRAMS = nctest tst_rename
check_PROGRAMS = $(TESTPROGRAMS)
TESTS = $(TESTPROGRAMS) compare_test_files.sh

# Note which tests depend on other tests. Necessary for make -j check.
TEST_EXTENSIONS = .sh

# compare_test_files depends on nctest executing first.
compare_test_files.log: nctest.log 

# These are the source files for the nctest program.
nctest_SOURCES = add.c add.h atttests.c cdftests.c dimtests.c driver.c	\
emalloc.c emalloc.h error.c error.h misctest.c rec.c slabs.c testcdf.h	\
tests.h val.c val.h vardef.c varget.c vargetg.c varput.c varputg.c	\
vartests.c vputget.c vputgetg.c

# If valgrind is present, add valgrind targets.
@VALGRIND_CHECK_RULES@
