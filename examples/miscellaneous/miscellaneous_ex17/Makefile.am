example_name  = miscellaneous_ex17
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/miscellaneous/ex17
data          = miscellaneous_ex17.C exact_solution.C run.sh
sources       = $(data)

CLEANFILES = out.pvtu

clean-local:
	-rm -rf out_*.vtu

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
