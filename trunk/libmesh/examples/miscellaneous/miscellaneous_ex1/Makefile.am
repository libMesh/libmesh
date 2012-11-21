check_PROGRAMS            = miscellaneous_ex1
miscellaneous_ex1dir      = $(examples_install_path)/miscellaneous/ex1
miscellaneous_ex1_DATA    = miscellaneous_ex1.C
miscellaneous_ex1_SOURCES = $(miscellaneous_ex1_DATA) run.sh
miscellaneous_ex1_LDADD   = $(top_builddir)/libmesh.la

# conditionally enable this test
if LIBMESH_ENABLE_INFINITE_ELEMENTS
  check_SCRIPTS = run.sh
endif

CLEANFILES = eqn_sys.dat ifems_added.e orig_mesh.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
