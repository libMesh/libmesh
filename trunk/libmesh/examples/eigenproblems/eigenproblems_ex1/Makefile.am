check_PROGRAMS            = eigenproblems_ex1
eigenproblems_ex1dir      = $(examples_install_path)/eigenproblems/ex1
eigenproblems_ex1_DATA    = eigenproblems_ex1.C $(top_srcdir)/contrib/utils/Makefile
eigenproblems_ex1_SOURCES = $(eigenproblems_ex1_DATA) run.sh
eigenproblems_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
