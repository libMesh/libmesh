## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.
# This file builds and runs the ncgen program.

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/lib_flags.am
LDADD = ${top_builddir}/liblib/libnetcdf.la

# Build ncgen from the listed sources.
bin_PROGRAMS = ncgen

ncgen_SOURCES=generate.c main.c cdata.c bindata.c genchar.c cvt.c data.c debug.c	\
escapes.c genc.c genbin.c generr.c genlib.c getfill.c odom.c	\
semantics.c dump.c util.c bytebuffer.c list.c data.h		\
debug.h generr.h genlib.h includes.h ncgen.h odom.h  dump.h	\
util.h bytebuffer.h list.h genf77.c f77data.c genj.c jdata.c nc_iter.h	\
nc_iter.c \
ncgeny.c ncgeny.h

# This is the man page.
man_MANS = ncgen.1

# These files all need to be distributed.
EXTRA_DIST = ncgen.y ncgen.l ncgenl.c $(man_MANS) internals.html	\
c0.cdl c0_4.cdl ref_camrun.cdl \
ncf199.cdl CMakeLists.txt XGetopt.c c5.cdl \
compound_datasize_test.cdl compound_datasize_test2.cdl \
tst_gattenum.cdl tst_usuffix.cdl

# This shell script causes ncgen to build a classic and a 64-bit
# offset file from a cdl file shipped with the distribution.
#if !BUILD_DLL
##
# Moved to ncdump as these tests rely on ncdump having been built.
##
#TESTS = run_tests.sh
#if USE_NETCDF4
#TESTS += run_nc4_tests.sh
#endif # USE_NETCDF4
#endif # !BUILD_DLL

CLEANFILES = c0.nc c0_64.nc c0_4.nc c0_4c.nc ref_camrun.c \
		ncf199.nc c5.nc compound_datasize_test.nc compound_datasize_test2.nc \
		tst_compound_datasize_test.cdl tst_compound_datasize_test2.cdl tst_ncf199.cdl \
		tst_tst_gattenum.cdl tst_gattenum.nc \
		tst_tst_usuffix.cdl tst_usuffix.nc \
		tst_c0.cdl tst_c0_4.cdl tst_c0_4c.cdl tst_c0_64.cdl

# These rules are used if someone wants to rebuild ncgenl.c or ncgeny.c
# Otherwise never invoked, but records how to do it.
# BTW: note that renaming is essential because otherwise
# autoconf will forcibly delete files of the name *.tab.*
# Note also that this should be built under linux or cygwin
# using bison version 3 or later.
# Note also that this code is in a shell script if you do
# not want to other to ./configure.

makeparser::
	flex -Pncg -8 ncgen.l
	rm -f ncgenl.c
	sed -e s/lex.ncg.c/ncgenl.c/g <lex.ncg.c >ncgenl.c
	bison -pncg -t -d ncgen.y
	rm -f ncgeny.c ncgeny.h
	sed -e s/ncgen.tab.c/ncgeny.c/g -e s/ncgen.tab.h/ncgeny.h/g <ncgen.tab.c >ncgeny.c
	sed -e s/ncgen.tab.c/ncgeny.c/g -e s/ncgen.tab.h/ncgeny.h/g <ncgen.tab.h >ncgeny.h
