# This is part of the netCDF package. Copyright 2005-2018 University
# Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
# conditions of use.
#
# This directory holds tests for netCDF-4. It is skipped if netCDF-4
# is not enabled.
#
# Ed Hartnett, Ward Fisher

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/lib_flags.am

# Un comment to use a more verbose test driver
#SH_LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
TEST_EXTENSIONS = .sh

# Note which tests depend on other tests. necessary for make -j check
extradir =

# Link to our assembled library.
AM_LDFLAGS += ${top_builddir}/liblib/libnetcdf.la
LDADD = ${top_builddir}/liblib/libnetcdf.la

# These are netCDF-4 C test programs which are built and run.
NC4_TESTS = tst_dims tst_dims2 tst_dims3 tst_files tst_files4 		\
tst_vars tst_varms tst_unlim_vars tst_converts tst_converts2 tst_grps	\
tst_grps2 tst_compounds tst_compounds2 tst_compounds3 tst_opaques	\
tst_strings tst_strings2 tst_interops tst_interops4 tst_interops5	\
tst_interops6 tst_enums tst_coords tst_coords2 tst_coords3 tst_vars3	\
tst_vars4 tst_chunks tst_chunks2 tst_utf8 tst_fills tst_fills2		\
tst_fillbug tst_xplatform tst_xplatform2 tst_endian_fill tst_atts	\
t_type cdm_sea_soundings tst_camrun tst_vl tst_atts1 tst_atts2		\
tst_vars2 tst_files5 tst_files6 tst_sync tst_h_scalar tst_rename	\
tst_rename2 tst_h5_endians tst_atts_string_rewrite			\
tst_hdf5_file_compat tst_fill_attr_vanish tst_rehash tst_filterparser	\
tst_bug324 tst_types tst_atts3 tst_put_vars tst_elatefill tst_udf

# Temporary I hope
if !ISCYGWIN
NC4_TESTS += tst_h_strbug tst_h_refs
endif

check_PROGRAMS =
TESTS =

# Build test programs plus programs used in test scripts.
if USE_NETCDF4
check_PROGRAMS += $(NC4_TESTS) renamegroup tst_empty_vlen_unlim
TESTS += $(NC4_TESTS)
endif

# Add these if large file tests are turned on.
if LARGE_FILE_TESTS
check_PROGRAMS += tst_large tst_large2
TESTS += tst_large tst_large2
endif

if USE_NETCDF4
if BUILD_UTILITIES
TESTS += run_grp_rename.sh tst_misc.sh
endif

TESTS += run_empty_vlen_test.sh
endif

# If the v2 API was built, add its test program.
if BUILD_V2
check_PROGRAMS += tst_v2
TESTS += tst_v2
endif # BUILD_V2

# If benchmarks were turned on, build and run a bunch more tests.
if BUILD_BENCHMARKS
check_PROGRAMS += tst_create_files bm_file tst_chunks3 tst_ar4	\
tst_ar4_3d tst_ar4_4d bm_many_objs tst_h_many_atts bm_many_atts	\
tst_files2 tst_files3 tst_mem tst_knmi bm_netcdf4_recs

bm_netcdf4_recs_SOURCES = bm_netcdf4_recs.c tst_utils.c
bm_many_atts_SOURCES = bm_many_atts.c tst_utils.c
bm_many_objs_SOURCES = bm_many_objs.c tst_utils.c
tst_ar4_3d_SOURCES = tst_ar4_3d.c tst_utils.c
tst_ar4_4d_SOURCES = tst_ar4_4d.c tst_utils.c
tst_files2_SOURCES = tst_files2.c tst_utils.c
tst_ar4_SOURCES = tst_ar4.c tst_utils.c
tst_h_many_atts_SOURCES = tst_h_many_atts.c tst_utils.c
bm_file_SOURCES = bm_file.c tst_utils.c
tst_knmi_SOURCES = tst_knmi.c tst_utils.c

TESTS += tst_ar4_3d tst_create_files run_bm_test1.sh run_bm_elena.sh	\
run_bm_test2.sh run_tst_chunks.sh tst_files2 tst_files3 tst_mem		\
run_knmi_bm.sh

# tst_create_files creates files for other tests.
run_bm_test1.log: tst_create_files.log
run_bm_test2.log: tst_create_files.log
run_bm_elena.log: tst_create_files.log

# This will run a parallel I/O benchmark for parallel builds.
if TEST_PARALLEL4
TESTS += run_par_bm_test.sh
# This benchmark depends on tst_create_files being run.
run_par_bm_test.log: tst_create_files.log
endif # TEST_PARALLEL4
endif # BUILD_BENCHMARKS

# Szip Tests (requires ncdump)
if USE_SZIP
if BUILD_UTILITIES
check_PROGRAMS += test_szip h5testszip
TESTS += tst_szip.sh
endif
endif

# Filter Tests (requires ncdump and ncgen)
if ENABLE_FILTER_TESTING
if BUILD_UTILITIES
extra_PROGRAMS = test_filter test_filter_misc
TESTS += tst_filter.sh
endif
endif

# This are extra tests that will only be run if netcdf-4 is configured
# with --enable-parallel-tests.
if TEST_PARALLEL4
check_PROGRAMS += tst_mpi_parallel tst_parallel tst_parallel3		\
tst_parallel4 tst_parallel5 tst_nc4perf tst_mode tst_simplerw_coll_r	\
tst_mode
TESTS += run_par_test.sh
endif

if ENABLE_METADATA_PERF
check_PROGRAMS += bigmeta openbigmeta tst_attsperf
TESTS += tst_attsperf perftest.sh
endif

EXTRA_DIST = run_par_test.sh.in run_par_bm_test.sh.in run_bm_test1.sh	\
run_bm_test2.sh run_bm_elena.sh run_tst_chunks.sh			\
ref_tst_compounds.nc ref_tst_xplatform2_1.nc ref_tst_xplatform2_2.nc	\
ref_tst_dims.nc ref_tst_interops4.nc run_knmi_bm.sh CMakeLists.txt	\
run_grp_rename.sh tst_h5_endians.c tst_atts_string_rewrite.c		\
tst_put_vars_two_unlim_dim.c tst_empty_vlen_unlim.c			\
run_empty_vlen_test.sh ref_hdf5_compat1.nc ref_hdf5_compat2.nc		\
ref_hdf5_compat3.nc tst_misc.sh tdset.h5 tst_szip.sh ref_szip.h5	\
ref_szip.cdl tst_filter.sh bzip2.cdl filtered.cdl unfiltered.cdl	\
ref_bzip2.c findplugin.in perftest.sh

CLEANFILES = tst_mpi_parallel.bin cdm_sea_soundings.nc bm_chunking.nc	\
tst_floats_1D.cdl floats_1D_3.nc floats_1D.cdl tst_*.nc			\
tst_floats2_*.cdl tst_ints2_*.cdl tst_shorts2_*.cdl tst_elena_*.cdl	\
tst_simple*.cdl tst_chunks.cdl pr_A1.* tauu_A1.* usi_01.* thetau_01.*	\
tst_*.h5 tst_grp_rename.cdl tst_grp_rename.dmp ref_grp_rename.cdl	\
foo1.nc tst_*.h4 test.nc testszip.nc test.h5 szip_dump.cdl		\
perftest.txt bigmeta.nc bigvars.nc run_par_test.sh *.gz MSGCPP_*.nc	\
floats*.nc floats*.cdl shorts*.nc shorts*.cdl ints*.nc ints*.cdl

DISTCLEANFILES = findplugin.sh

# If valgrind is present, add valgrind targets.
@VALGRIND_CHECK_RULES@
