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

check_PROGRAMS       = adjoints_ex2
check_SCRIPTS        = run.sh
adjoints_ex2dir      = $(examples_install_path)/adjoints/ex2
adjoints_ex2_DATA    = adjoints_ex2.C ../adjoints_ex1/lshaped.xda L-shaped.C L-shaped.h femparameters.C femparameters.h general.in l-shaped.in L-qoi.C L-qoi.h
adjoints_ex2_SOURCES = $(adjoints_ex2_DATA) run.sh
adjoints_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = primal.out.* adjoint_0.out.*

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
