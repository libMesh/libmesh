BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS         = subdomains_ex1
check_SCRIPTS          = run.sh
subdomains_ex1dir      = $(examples_install_path)/subdomains/ex1
subdomains_ex1_DATA    = subdomains_ex1.C exact_solution.C run.sh $(top_srcdir)/contrib/utils/Makefile
subdomains_ex1_SOURCES = $(subdomains_ex1_DATA) run.sh
subdomains_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.gmv out_*.e gnuplot_*

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
