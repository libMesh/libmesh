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

check_PROGRAMS        = vector_fe_ex1
check_SCRIPTS         = run.sh
vector_fe_ex1dir      = $(examples_install_path)/vector_fe/ex1
vector_fe_ex1_DATA    = vector_fe_ex1.C exact_solution.C
vector_fe_ex1_SOURCES = $(vector_fe_ex1_DATA) run.sh
vector_fe_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.gmv out.e 

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


