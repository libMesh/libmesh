#SUBDIRS = introduction

include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)
AM_INCLUDES = $(libmesh_INCLUDE)
LIBS        = $(libmesh_LIBS) 

check_PROGRAMS =
check_SCRIPTS = 
CLEANFILES =

#
# Introductory examples
#
check_PROGRAMS += introduction_ex1
check_SCRIPTS  += introduction/introduction_ex1/run.sh
introduction_ex1_SOURCES = introduction/introduction_ex1/introduction_ex1.C introduction/introduction_ex1/run.sh
introduction_ex1_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES += output.xda

check_PROGRAMS += introduction_ex2
check_SCRIPTS  += introduction/introduction_ex2/run.sh
introduction_ex2_SOURCES = introduction/introduction_ex2/introduction_ex2.C introduction/introduction_ex2/run.sh
introduction_ex2_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES += eqn_sys.dat

check_PROGRAMS += introduction_ex3
check_SCRIPTS  += introduction/introduction_ex3/run.sh
introduction_ex3_SOURCES = introduction/introduction_ex3/introduction_ex3.C introduction/introduction_ex3/exact_solution.C introduction/introduction_ex3/run.sh
introduction_ex3_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += introduction_ex4
check_SCRIPTS  += introduction/introduction_ex4/run.sh
introduction_ex4_SOURCES = introduction/introduction_ex4/introduction_ex4.C introduction/introduction_ex3/exact_solution.C introduction/introduction_ex4/run.sh
introduction_ex4_LDADD   = $(top_builddir)/libmesh.la




#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
