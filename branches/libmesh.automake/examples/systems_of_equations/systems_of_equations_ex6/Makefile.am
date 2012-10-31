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

check_PROGRAMS                   = systems_of_equations_ex6
check_SCRIPTS                    = run.sh
systems_of_equations_ex6dir      = $(examples_install_path)/systems_of_equations/ex6
systems_of_equations_ex6_DATA    = systems_of_equations_ex6.C
systems_of_equations_ex6_SOURCES = $(systems_of_equations_ex6_DATA) run.sh
systems_of_equations_ex6_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = displacement.e

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


