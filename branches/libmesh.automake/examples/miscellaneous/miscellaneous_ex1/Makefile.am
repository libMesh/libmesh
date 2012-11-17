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

check_PROGRAMS            = miscellaneous_ex1
miscellaneous_ex1dir      = $(examples_install_path)/miscellaneous/ex1
miscellaneous_ex1_DATA    = miscellaneous_ex1.C
miscellaneous_ex1_SOURCES = $(miscellaneous_ex1_DATA) run.sh
miscellaneous_ex1_LDADD   = $(top_builddir)/libmesh.la

if LIBMESH_ENABLE_INFINITE_ELEMENTS
  check_SCRIPTS = run.sh
endif

CLEANFILES = eqn_sys.dat ifems_added.e orig_mesh.e

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


