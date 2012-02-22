#SUBDIRS = introduction

include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)
LIBS        = $(top_builddir)/libmesh.la $(libmesh_LIBS)

INCLUDES = $(libmesh_INCLUDE)

check_PROGRAMS =
check_SCRIPTS = 
EXTRA_DIST =

#
# Introductory examples
#
check_PROGRAMS += introduction_ex1
check_SCRIPTS  += introduction/introduction_ex1/run.sh
introduction_ex1_SOURCES = introduction/introduction_ex1/introduction_ex1.C introduction/introduction_ex1/run.sh

check_PROGRAMS += introduction_ex2
check_SCRIPTS  += introduction/introduction_ex2/run.sh
introduction_ex2_SOURCES = introduction/introduction_ex2/introduction_ex2.C introduction/introduction_ex2/run.sh





TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
