example_name = miscellaneous_ex1
install_dir  = $(examples_install_path)/miscellaneous/ex1
data         = miscellaneous_ex1.C run.sh
sources      = $(data) run.sh

# conditionally enable this test
if LIBMESH_ENABLE_INFINITE_ELEMENTS
  check_SCRIPTS = run.sh
endif

CLEANFILES = eqn_sys.dat ifems_added.e orig_mesh.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
