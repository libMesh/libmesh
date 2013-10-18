BUILT_SOURCES = mesh.xda

mesh.xda: $(top_srcdir)/examples/adaptivity/adaptivity_ex2/mesh.xda
	$(AM_V_GEN)$(LN_S) $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

example_name  = transient_ex1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/transient/ex1
data          = transient_ex1.C exact_solution.C ../../adaptivity/adaptivity_ex2/mesh.xda run.sh
sources       = $(data) run.sh

CLEANFILES = out_*.gmv transient_ex1.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
