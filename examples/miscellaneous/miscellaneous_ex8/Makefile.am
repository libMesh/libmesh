example_name  = miscellaneous_ex8
install_dir   = $(examples_install_path)/miscellaneous/ex8
data          = miscellaneous_ex8.C run.sh
sources       = $(data) run.sh

# conditionally enable this test
if LIBMESH_ENABLE_NANOFLANN
  check_SCRIPTS = run.sh
endif

CLEANFILES = 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
