example_name  = miscellaneous_ex4
install_dir   = $(examples_install_path)/miscellaneous/ex4
data          = miscellaneous_ex4.C run.sh
sources       = $(data) run.sh
check_SCRIPTS = run.sh

CLEANFILES = out.pvtu

clean-local:
	-rm -rf out_*.vtu

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
