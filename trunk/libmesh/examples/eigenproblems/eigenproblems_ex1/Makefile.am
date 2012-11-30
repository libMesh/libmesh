example_name = eigenproblems_ex1
install_dir  = $(examples_install_path)/eigenproblems/ex1
data         = eigenproblems_ex1.C
sources      = $(data) run.sh

CLEANFILES = out.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
