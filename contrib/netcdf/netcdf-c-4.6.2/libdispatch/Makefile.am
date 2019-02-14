## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.

# This Makefile controls the building of the dispatch layer of the
# netCDF library. The dispatch layer decides whether to call the
# netcdf-classic code, netcdf-4 code, nc3 dap code, or nc4 dap
# code. It also contains code that sit above the dispatch layer, like
# the v2 API.

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/lib_flags.am

# This is our output, the dispatch convenience library.
noinst_LTLIBRARIES = libdispatch.la
libdispatch_la_CPPFLAGS = ${AM_CPPFLAGS}

# The source files.
libdispatch_la_SOURCES = dparallel.c dcopy.c dfile.c ddim.c datt.c	\
dattinq.c dattput.c dattget.c derror.c dvar.c dvarget.c dvarput.c	\
dvarinq.c dinternal.c ddispatch.c dutf8.c nclog.c dstring.c ncuri.c	\
nclist.c ncbytes.c nchashmap.c nctime.c nc.c nclistmgr.c drc.c		\
dauth.c doffsets.c dwinpath.c dutil.c dreadonly.c dnotnc4.c dnotnc3.c	\
crc32.c crc32.h daux.c

# Add the utf8 codebase
libdispatch_la_SOURCES += utf8proc.c utf8proc.h

# Add functions only found in netCDF-4.
if USE_NETCDF4
libdispatch_la_SOURCES += dgroup.c dvlen.c dcompound.c dtype.c denum.c	\
dopaque.c dfilter.c
endif # USE_NETCDF4

# Add V2 API convenience library if needed.
if BUILD_V2
noinst_LTLIBRARIES += libnetcdf2.la
libnetcdf2_la_SOURCES = dv2i.c
libnetcdf2_la_CPPFLAGS = ${AM_CPPFLAGS} -DDLL_EXPORT
endif # BUILD_V2

EXTRA_DIST=CMakeLists.txt ncsettings.hdr utf8proc_data.c

# Build ncsettings.c as follows:
# 1. copy ncsettings.hdr to ncsettings.c
# 2. append libnetcdf.settings to ncsettings.c after
#    processing it as follows:
#    1. convert tabs and cr to blanks
#    2. convert embedded double quote (") to escaped form (\").
#    3. append newline (\n) to each line
#    4. surround each line with double quotes.
# 3. finally, add a semicolon to the end of ncsettings.c
#    to complete the string constant.

ncsettings.c: $(top_srcdir)/libnetcdf.settings ncsettings.hdr
	rm -f ncsettings.c
	cat ncsettings.hdr > ncsettings.c
	tr '\t\r' '  ' <${top_srcdir}/libnetcdf.settings | \
	sed -e 's/"/\\"/g' | \
	sed -e 's/\(.*\)/\"\1\\n\"/' | \
	cat >> ncsettings.c
	echo ';' >> ncsettings.c

# Unit tests
test_ncuri_SOURCES = test_ncuri.c ncuri.c ncbytes.c nclist.c
test_pathcvt_SOURCES = test_pathcvt.c dwinpath.c
check_PROGRAMS = test_ncuri test_pathcvt
TESTS = test_ncuri test_pathcvt
