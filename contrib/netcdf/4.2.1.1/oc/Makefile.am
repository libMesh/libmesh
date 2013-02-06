## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.

# This automake file generates the Makefile to build netCDF-4. The
# generated makefile is not run unless the user selected to build
# netCDF-4.

# $Id: Makefile.am,v 1.1 2010/05/23 21:05:33 dmh Exp $

# Cause C preprocessor to search current and parent directory.
AM_CPPFLAGS =  -I$(top_srcdir) -I$(top_srcdir)/include

# OC Sources; include the daptab.[ch] to avoid the need for bison by user
SRC= occurlfunctions.c dapparse.c daplex.c ocuri.c ochttp.c ocbytes.c	\
oc.c occontent.c ocdebug.c ocdrno.c ocdump.c	\
ocinternal.c oclist.c oclog.c ocnode.c ocutil.c occlientparams.c	\
ocrc.c ocread.c xxdr.c daptab.c

HDRS= occurlfunctions.h occonstraints.h dapparselex.h ocuri.h ochttp.h \
ocbytes.h occontent.h ocdata.h ocdatatypes.h ocdebug.h ocdrno.h ocdump.h \
ocinternal.h oclist.h oclog.h ocnode.h occlientparams.h ocutil.h ocrc.h	\
ocread.h xxdr.h daptab.h \
oc.h

if BUILD_DAP
noinst_LTLIBRARIES = liboc.la
liboc_la_SOURCES = $(SRC) $(HDRS) 
liboc_la_CPPFLAGS = $(AM_CPPFLAGS)
endif

# These rule are used if someone wants to rebuild the grammar files.
# Otherwise never invoked, but records how to do it.
# BTW: note that renaming is essential because otherwise
# autoconf will forcibly delete files of the name *.tab.*

.PHONEY: bison

bison:: dap.y
	rm -f dap.tab.c dap.tab.h
	bison --debug -d -p dap dap.y
	mv dap.tab.c daptab.c; mv dap.tab.h daptab.h

