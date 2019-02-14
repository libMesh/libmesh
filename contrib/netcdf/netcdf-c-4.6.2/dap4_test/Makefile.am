## This is a automake file, part of Unidata's netCDF package.
# Copyright 2011, see the COPYRIGHT file for more information.

# This file builds and runs DAP4 tests.

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/lib_flags.am

#TEST_EXTENSIONS = .sh

#SH_LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#sh_LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#LOG_DRIVER = $(SHELL) $(top_srcdir)/test-driver-verbose
#TESTS_ENVIRONMENT = export SETX=1;
	
LDADD = ${top_builddir}/liblib/libnetcdf.la
AM_CPPFLAGS += -I$(top_srcdir)/libdap4

# Set up the tests; do the .sh first, then .c
check_PROGRAMS =
TESTS =

if ENABLE_DAP4

# WARNING: these are unit tests, so they will not
# appear in CMakeLists.txt
check_PROGRAMS += test_parse test_meta test_data
test_parse_SOURCES = test_parse.c test_common.h
test_meta_SOURCES = test_meta.c test_common.h
test_data_SOURCES = test_data.c test_common.h
TESTS += test_parse.sh

if BUILD_UTILITIES
# These rely on ncdump
TESTS += test_raw.sh
TESTS += test_meta.sh
TESTS += test_data.sh
TESTS += test_fillmismatch.sh
endif

# Note which tests depend on other tests. Necessary for make -j check.
test_raw.log: test_parse.log
test_meta.log: test_raw.log
test_data.log: test_meta.log

if ENABLE_DAP_REMOTE_TESTS
# Note: This program name was changed to findtestserver4
# to avoid cmake complaint about duplicate targets.
check_PROGRAMS += findtestserver4
findtestserver4_SOURCES = findtestserver4.c
if BUILD_UTILITIES
  # relies on ncdump
  TESTS += test_remote.sh
endif
endif

endif #ENABLE_DAP4

EXTRA_DIST = test_parse.sh test_meta.sh test_data.sh \
             test_raw.sh test_remote.sh test_hyrax.sh test_fillmismatch.sh \
             tst_curlopt.sh d4test_common.sh \
	     daptestfiles dmrtestfiles cdltestfiles nctestfiles misctestfiles \
	     baseline baselineraw baselineremote CMakeLists.txt

CLEANFILES = *.exe
# This should only be left behind if using parallel io
CLEANFILES += tmp_*

DISTCLEANFILES = findtestserver4.c

# One last thing
BUILT_SOURCES = .daprc

.daprc:
	echo "#DAPRC" >.daprc

clean-local: clean-local-check

.PHONY: clean-local-check

clean-local-check:
	-rm -rf results
	-rm -f .dodsrc .daprc

# The shell file maketests.sh is used to build the testdata
# for dap4 testing. It creates and fills the directories
# daptestfiles, dmrtestfiles, and cdltestfiles
# As a rule, this will only be invoked if there was a change
# in the GenerateRaw.java program in the thredds tree.
maketests::
	bash ./maketests.sh
