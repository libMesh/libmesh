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

check_PROGRAMS            = eigenproblems_ex3
eigenproblems_ex3dir      = $(examples_install_path)/eigenproblems/ex3
eigenproblems_ex3_DATA    = eigenproblems_ex3.C drum1_mesh.e drum2_mesh.e
eigenproblems_ex3_SOURCES = $(eigenproblems_ex3_DATA) run.sh
eigenproblems_ex3_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES                = drum1_evals.txt drum2_evals.txt

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
