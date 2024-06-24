example_name  = miscellaneous_ex16
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/miscellaneous/ex16
data          = miscellaneous_ex16.C exact_solution.C run.sh
sources       = $(data) run.sh

CLEANFILES = out.pvtu

clean-local:
	-rm -rf out_*.vtu

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
