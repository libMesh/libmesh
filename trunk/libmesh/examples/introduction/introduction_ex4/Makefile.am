BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS           = introduction_ex4
check_SCRIPTS            = run.sh
introduction_ex4dir      = $(examples_install_path)/introduction/ex4
introduction_ex4_DATA    = introduction_ex4.C exact_solution.C 
introduction_ex4_SOURCES = $(introduction_ex4_DATA) run.sh
introduction_ex4_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES               = gnuplot_script gnuplot_script_data out_2.e out_3.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
