example_name = miscellaneous_ex15
install_dir  = $(examples_install_path)/miscellaneous/ex15
data         = miscellaneous_ex15.C run.sh
sources      = $(data) run.sh

# conditionally enable this test
if LIBMESH_ENABLE_INFINITE_ELEMENTS
  check_SCRIPTS = run.sh
endif

CLEANFILES = ""

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
