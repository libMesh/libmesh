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

check_PROGRAMS            = reduced_basis_ex2
reduced_basis_ex2dir      = $(examples_install_path)/reduced_basis/ex2
reduced_basis_ex2_DATA    = reduced_basis_ex2.C assembly.h rb_classes.h reduced_basis_ex2.in
reduced_basis_ex2_SOURCES = $(reduced_basis_ex2_DATA) run.sh
reduced_basis_ex2_LDADD   = $(top_builddir)/libmesh.la

# Only *run* the tests if SLEPc and GLPK is actually enabled
if LIBMESH_ENABLE_SLEPC
if LIBMESH_ENABLE_GLPK
  check_SCRIPTS = run.sh
endif
endif

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
