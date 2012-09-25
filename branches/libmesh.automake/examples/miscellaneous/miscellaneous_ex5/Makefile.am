include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_INCLUDE) -I$(top_builddir)/include
LIBS         = $(libmesh_LIBS)

if LIBMESH_ENABLE_FPARSER
  AM_CPPFLAGS += -I$(top_srcdir)/contrib/fparser
endif

check_PROGRAMS            = miscellaneous_ex5
miscellaneous_ex5dir      = $(examples_install_path)/miscellaneous/ex5
miscellaneous_ex5_DATA    = miscellaneous_ex5.C ../../adjoints/adjoints_ex1/lshaped.xda lshaped3D.xda miscellaneous_ex5.in
miscellaneous_ex5_SOURCES = $(miscellaneous_ex5_DATA) run.sh
miscellaneous_ex5_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = lshaped_dg.e 

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
