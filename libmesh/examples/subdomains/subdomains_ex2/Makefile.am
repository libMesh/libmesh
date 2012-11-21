BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS         = subdomains_ex2
check_SCRIPTS          = run.sh
subdomains_ex2dir      = $(examples_install_path)/subdomains/ex2
subdomains_ex2_DATA    = subdomains_ex2.C exact_solution.C
subdomains_ex2_SOURCES = $(subdomains_ex2_DATA) run.sh
subdomains_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e gnuplot_*

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
