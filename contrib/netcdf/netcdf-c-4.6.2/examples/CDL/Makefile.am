## This is a automake file, part of Unidata's netCDF package.
# Copyright 2006, see the COPYRIGHT file for more information.
# This file builds CDL examples.
# Ed Hartnett

# Note which tests depend on other tests. Necessary for make -j check.
TEST_EXTENSIONS = .sh

# All we do is run two shell scripts, which create the sample files
# and then compare them with the C versions of the examples.
TESTS = do_comps.sh

# Ship the scripts needed to create the sample files and compare them.
EXTRA_DIST = do_comps.sh simple_xy.cdl	\
sfc_pres_temp.cdl pres_temp_4D.cdl CMakeLists.txt

# Clean up files created during the process.
CLEANFILES = *.nc
