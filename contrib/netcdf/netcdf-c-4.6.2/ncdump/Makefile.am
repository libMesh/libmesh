## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.
# This file builds and runs the ncdump program.

# Ed Hartnett, Dennis Heimbigner, Ward Fisher

#SH_LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#sh_LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#TESTS_ENVIRONMENT = export SETX=1;

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/lib_flags.am
LDADD = ${top_builddir}/liblib/libnetcdf.la
noinst_PROGRAMS=
# Note which tests depend on other tests. Necessary for make -j check.
TEST_EXTENSIONS = .sh
XFAIL_TESTS=""

# This is the program we're building, and it's sources.
bin_PROGRAMS = ncdump
ncdump_SOURCES = ncdump.c vardata.c dumplib.c indent.c nctime0.c        \
ncdump.h vardata.h dumplib.h indent.h isnan.h nctime0.h cdl.h utils.h   \
utils.c nciter.h nciter.c nccomps.h

# Another utility program that copies any netCDF file using only the
# netCDF API
bin_PROGRAMS += nccopy
nccopy_SOURCES = nccopy.c nciter.c nciter.h chunkspec.h chunkspec.c     \
utils.h utils.c dimmap.h dimmap.c list.c list.h

# A simple netcdf-4 metadata -> xml printer. Do not install.
if USE_NETCDF4
noinst_PROGRAMS += nc4print
nc4print_SOURCES = nc4print.c
endif

# Conditionally build the ocprint program, but do not install
if ENABLE_DAP
noinst_PROGRAMS += ocprint
ocprint_SOURCES = ocprint.c
endif

# This is the man page.
man_MANS = ncdump.1 nccopy.1

if BUILD_TESTSETS
# C programs needed by shell scritps for classic tests.
check_PROGRAMS = rewrite-scalar ref_ctest ref_ctest64 ncdump tst_utf8   \
bom tst_dimsizes nctrunc

# Tests for classic and 64-bit offset files.
TESTS = tst_inttags.sh run_tests.sh tst_64bit.sh ref_ctest	\
ref_ctest64 tst_output.sh tst_lengths.sh tst_calendars.sh	\
run_utf8_tests.sh tst_nccopy3.sh tst_nccopy3_subset.sh		\
tst_charfill.sh tst_iter.sh tst_formatx3.sh tst_bom.sh		\
tst_dimsizes.sh run_ncgen_tests.sh tst_ncgen4_classic.sh

# The tst_nccopy3.sh test uses output from a bunch of other
# tests. This records the dependency so parallel builds work.
tst_nccopy3.log: tst_calendars.log run_utf8_tests.log tst_output.log	\
tst_64bit.log run_tests.log tst_lengths.log

TESTS += tst_null_byte_padding.sh
if USE_STRICT_NULL_BYTE_HEADER_PADDING
XFAIL_TESTS += tst_null_byte_padding.sh
endif

if LARGE_FILE_TESTS
TESTS += tst_iter.sh
endif

if USE_NETCDF4
# NetCDF-4 has some extra C programs to build. These will be run by
# the shell script tests.
check_PROGRAMS += tst_fileinfo tst_create_files tst_h_rdc0	\
tst_group_data tst_enum_data tst_opaque_data tst_string_data	\
tst_vlen_data tst_comp tst_comp2 tst_nans tst_special_atts	\
tst_unicode tst_fillbug tst_compress tst_chunking tst_h_scalar

check_PROGRAMS += tst_vlen_demo

# Tests for netCDF-4 behavior.
TESTS += tst_fileinfo.sh tst_hdf5_offset.sh tst_inttags4.sh		\
tst_netcdf4.sh tst_fillbug.sh tst_netcdf4_4.sh tst_nccopy4.sh		\
tst_nccopy5.sh tst_grp_spec.sh tst_mud.sh tst_h_scalar.sh tst_formatx4.sh		\
run_utf8_nc4_tests.sh run_back_comp_tests.sh run_ncgen_nc4_tests.sh	\
tst_ncgen4.sh

# Record interscript dependencies so parallel builds work.
tst_nccopy4.log: run_ncgen_tests.log tst_output.log tst_ncgen4.log	\
tst_fillbug.log tst_netcdf4_4.log tst_h_scalar.log
tst_nccopy5.log: tst_nccopy4.log
endif #!USE_NETCDF4

TESTS += tst_inmemory_nc3.sh
if USE_NETCDF4
TESTS += tst_inmemory_nc4.sh
endif

else # when BUILD_TESTSETS is false, we still need this rule
ctest.c:
	cp $(top_srcdir)/ncdump/ref_ctest.c $(top_srcdir)/ncdump/ctest.c

ctest64.c:
	cp $(top_srcdir)/ncdump/ref_ctest64.c $(top_srcdir)/ncdump/ctest64.c

endif BUILD_TESTSETS_FALSE

# These files all have to be included with the distribution.
EXTRA_DIST = run_tests.sh tst_64bit.sh tst_output.sh test0.cdl		\
ref_ctest1_nc4.cdl ref_ctest1_nc4c.cdl ref_tst_solar_1.cdl		\
ref_tst_solar_2.cdl tst_netcdf4.sh tst_netcdf4_4.sh ref_tst_small.cdl	\
tst_lengths.sh tst_ncml.cdl ref1.ncml ref_tst_group_data.cdl		\
ref_tst_enum_data.cdl ref_tst_opaque_data.cdl ref_tst_string_data.cdl	\
ref_tst_vlen_data.cdl ref_tst_comp.cdl ref_tst_unicode.cdl		\
ref_tst_nans.cdl small.cdl small2.cdl $(man_MANS) run_utf8_tests.sh	\
ref_tst_utf8.cdl ref_tst_fillbug.cdl tst_fillbug.sh tst_calendars.cdl	\
tst_calendars.sh ref_times.cdl ref_tst_special_atts.cdl			\
ref_tst_noncoord.cdl ref_tst_compounds2.nc ref_tst_compounds2.cdl	\
ref_tst_compounds3.nc ref_tst_compounds3.cdl ref_tst_compounds4.nc	\
ref_tst_compounds4.cdl ref_tst_group_data_v23.cdl tst_mslp.cdl		\
tst_bug321.cdl ref_tst_format_att.cdl ref_tst_format_att_64.cdl		\
tst_nccopy3.sh tst_nccopy4.sh tst_nccopy5.sh ref_nc_test_netcdf4_4_0.nc		\
run_back_comp_tests.sh ref_nc_test_netcdf4.cdl				\
ref_tst_special_atts3.cdl tst_brecs.cdl ref_tst_grp_spec0.cdl		\
ref_tst_grp_spec.cdl tst_grp_spec.sh ref_tst_charfill.cdl		\
tst_charfill.cdl tst_charfill.sh tst_iter.sh tst_mud.sh			\
ref_tst_mud4.cdl ref_tst_mud4-bc.cdl ref_tst_mud4_chars.cdl		\
inttags.cdl inttags4.cdl ref_inttags.cdl ref_inttags4.cdl		\
ref_tst_ncf213.cdl tst_h_scalar.sh run_utf8_nc4_tests.sh		\
tst_formatx3.sh tst_formatx4.sh ref_tst_utf8_4.cdl			\
ref_tst_nc4_utf8_4.cdl tst_inttags.sh tst_inttags4.sh CMakeLists.txt	\
XGetopt.c tst_bom.sh tst_inmemory_nc3.sh tst_dimsizes.sh		\
tst_inmemory_nc4.sh tst_fileinfo.sh run_ncgen_tests.sh			\
ref_test_360_day_1900.nc ref_test_365_day_1900.nc			\
ref_test_366_day_1900.nc ref_test_360_day_1900.cdl			\
ref_test_365_day_1900.cdl ref_test_366_day_1900.cdl			\
tst_hdf5_offset.sh run_ncgen_nc4_tests.sh tst_nccopy3_subset.sh		\
ref_nccopy3_subset.nc ref_test_corrupt_magic.nc tst_ncgen_shared.sh	\
tst_ncgen4.sh tst_ncgen4_classic.sh tst_ncgen4_diff.sh			\
tst_ncgen4_cycle.sh tst_null_byte_padding.sh				\
ref_null_byte_padding_test.nc ref_tst_irish_rover.nc ref_provenance_v1.nc

# The L512.bin file is file containing exactly 512 bytes each of value 0.
# It is used for creating hdf5 files with varying offsets for testing.
EXTRA_DIST += L512.bin

# CDL files and Expected results
SUBDIRS = cdl expected

CLEANFILES = tst_*.nc tmp*.nc test*.nc iter.* tmp*.cdl			\
tst_output_*.cdl tst_output_*.c tst_utf8_*.cdl tst_tst8.cdl		\
tst_netcdf4_*.cdl test1_ncdump.cdl test2_ncdump.cdl test1.cdl		\
ctest1.cdl test1_cdf5.cdl test2_cdf5.cdl test1_offset.cdl		\
test2_offset.cdl ctest0.nc ctest0_64.nc c1.cdl c1_4.cdl ctest1_64.cdl	\
c0.nc c0_4.nc small.nc small2.nc c0tmp.nc c1.ncml utf8.cdl		\
utf8_64.cdl utf8.nc utf8_64.nc nc4_utf8.cdl nc4_utf8.nc			\
tst_unicode.cdl tst_group_data.cdl tst_compounds2.cdl tst_comp.cdl	\
tst_enum_data.cdl tst_small.cdl tst_times.cdl tst_solar_2.cdl		\
tst_string_data.cdl tst_fillbug.cdl tst_opaque_data.cdl			\
tst_compounds4.cdl tst_utf8.cdl tst_compounds3.cdl			\
tst_special_atts.cdl tst_nans.cdl tst_format_att_64.cdl			\
tst_vlen_data.cdl tst_solar_1.cdl tst_format_att.cdl			\
tst_nc_test_netcdf4_4_0.cdl tst_mud4.cdl tst_mud4-bc.cdl		\
tst_ncf213.cdl tst_h_scalar.cdl tst_mud4_chars.cdl inttags.nc		\
inttags4.nc tst_inttags.cdl tst_inttags4.cdl nc4_fileinfo.nc		\
hdf5_fileinfo.hdf nccopy3_subset_out.nc c5.nc				\
compound_datasize_test.nc compound_datasize_test2.nc ncf199.nc		\
tst_c0.cdl tst_c0_4.cdl tst_c0_4c.cdl tst_c0_64.cdl			\
tst_compound_datasize_test.cdl tst_compound_datasize_test2.cdl		\
tst_ncf199.cdl tst_tst_gattenum.cdl tst_tst_usuffix.cdl ctest.c		\
ctest64.c nccopy3_subset_out.nc camrun.c
