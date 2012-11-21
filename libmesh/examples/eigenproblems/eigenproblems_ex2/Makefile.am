check_PROGRAMS            = eigenproblems_ex2
eigenproblems_ex2dir      = $(examples_install_path)/eigenproblems/ex2
eigenproblems_ex2_DATA    = eigenproblems_ex2.C
eigenproblems_ex2_SOURCES = $(eigenproblems_ex2_DATA) run.sh
eigenproblems_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
