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

check_PROGRAMS            = eigenproblems_ex2
eigenproblems_ex2dir      = $(examples_install_path)/eigenproblems/ex2
eigenproblems_ex2_DATA    = eigenproblems_ex2.C
eigenproblems_ex2_SOURCES = $(eigenproblems_ex2_DATA) run.sh
eigenproblems_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.e

# Only *run* the tests if SLEPc is actually enabled
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


