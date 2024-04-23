example_name = eigenproblems_ex4
install_dir  = $(examples_install_path)/eigenproblems/ex4
data         = eigenproblems_ex4.C
sources      = $(data) run.sh

CLEANFILES = out.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
