BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C
	$(AM_V_GEN)$(LN_S) -f $< $@

example_name  = introduction_ex5
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/introduction/ex5
data          = introduction_ex5.C exact_solution.C run.sh
sources       = $(data) run.sh
CLEANFILES    = out_0.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
