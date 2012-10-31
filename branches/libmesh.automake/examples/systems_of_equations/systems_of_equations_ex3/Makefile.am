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

check_PROGRAMS                   = systems_of_equations_ex3
check_SCRIPTS                    = run.sh
systems_of_equations_ex3dir      = $(examples_install_path)/systems_of_equations/ex3
systems_of_equations_ex3_DATA    = systems_of_equations_ex3.C
systems_of_equations_ex3_SOURCES = $(systems_of_equations_ex3_DATA) run.sh
systems_of_equations_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


