BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C
	$(AM_V_GEN)$(LN_S) -f $< $@

example_name  = introduction_ex4
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/introduction/ex4
data          = introduction_ex4.C exact_solution.C run.sh
sources       = $(data) run.sh
CLEANFILES    = gnuplot_script gnuplot_script_data out_2.e out_3.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
