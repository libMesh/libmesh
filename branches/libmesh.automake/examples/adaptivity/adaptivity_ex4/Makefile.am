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

check_PROGRAMS         = adaptivity_ex4
check_SCRIPTS          = run.sh
adaptivity_ex4dir      = $(examples_install_path)/adaptivity/ex4
adaptivity_ex4_DATA    = adaptivity_ex4.C adaptivity_ex4.in
adaptivity_ex4_SOURCES = $(adaptivity_ex4_DATA) run.sh
adaptivity_ex4_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = 2D_hermite_adaptive.m 2D_hermite_adaptive.e 

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


