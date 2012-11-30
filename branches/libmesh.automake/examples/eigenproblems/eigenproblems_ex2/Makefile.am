example_name = eigenproblems_ex2
install_dir  = $(examples_install_path)/eigenproblems/ex2
data         = eigenproblems_ex2.C run.sh
sources      = $(data) run.sh

CLEANFILES = out.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
