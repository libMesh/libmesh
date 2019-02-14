## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.

# This automake file generates the Makefile to build netCDF-4. The
# generated makefile is not run unless the user selected to build
# netCDF-4.

# $Id: Makefile.am,v 1.1 2010/05/23 21:05:33 dmh Exp $

# Cause C preprocessor to search current and parent directory.
AM_CPPFLAGS =  -I.. -I$(top_srcdir) -I$(top_srcdir)/include

# OC Sources; include the dapy.[ch] to avoid the need for bison by user
SRC=oc.c \
daplex.c dapparse.c dapy.c \
occompile.c occurlfunctions.c \
ocdata.c ocdebug.c ocdump.c  \
ocinternal.c ocnode.c \
ochttp.c \
ocread.c ocutil.c \
xxdr.c

HDRS=oc.h ocx.h \
dapparselex.h dapy.h \
occompile.h occonstraints.h occurlfunctions.h \
ocdata.h ocdatatypes.h ocdebug.h ocdump.h \
ocinternal.h ocnode.h \
ochttp.h ocread.h ocutil.h \
xxdr.h

EXTRA_DIST = dap.y CMakeLists.txt auth.html.in oc.css

noinst_LTLIBRARIES = liboc.la
liboc_la_SOURCES = $(SRC) $(HDRS) 
liboc_la_CPPFLAGS = $(AM_CPPFLAGS)

# These rule are used if someone wants to rebuild the grammar files.
# Otherwise never invoked, but records how to do it.
# BTW: note that renaming is essential because otherwise
# autoconf will forcibly delete files of the name *.tab.*

.PHONEY: bison

bison:: dap.y
	rm -f dap.tab.c dap.tab.h
	bison --debug -d -p dap dap.y
	mv dap.tab.c dapy.c; mv dap.tab.h dapy.h
