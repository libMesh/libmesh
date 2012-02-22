include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)
LIBS        = $(top_builddir)/libmesh.la $(libmesh_LIBS)

INCLUDES = $(libmesh_INCLUDE)

check_PROGRAMS =

check_PROGRAMS += introduction-ex1
introduction_ex1_SOURCES = introduction_ex1/introduction_ex1.C

#TESTS_ENVIRONMENT = 
TESTS = $(check_PROGRAMS)
