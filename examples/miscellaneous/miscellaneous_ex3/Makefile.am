# need to link lshaped.xda
BUILT_SOURCES = lshaped.xda
lshaped.xda:
	$(AM_V_GEN)$(LN_S) -f $(top_srcdir)/examples/adjoints/adjoints_ex1/lshaped.xda .

example_name  = miscellaneous_ex3
install_dir   = $(examples_install_path)/miscellaneous/ex3
data          = miscellaneous_ex3.C ../../adjoints/adjoints_ex1/lshaped.xda run.sh
sources       = $(data) run.sh
check_SCRIPTS = run.sh

CLEANFILES = lshaped.xda out.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
