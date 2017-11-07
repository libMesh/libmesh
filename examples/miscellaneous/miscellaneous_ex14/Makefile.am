example_name = miscellaneous_ex14
install_dir  = $(examples_install_path)/miscellaneous/ex14
data         = miscellaneous_ex14.C run.sh
sources      = $(data) run.sh

# conditionally enable this test
if LIBMESH_ENABLE_INFINITE_ELEMENTS
  check_SCRIPTS = run.sh
endif

CLEANFILES = U-0_inf.e *_infini_0.txt 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
