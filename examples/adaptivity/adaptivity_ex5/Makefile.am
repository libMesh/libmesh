BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/transient/transient_ex1/exact_solution.C
	$(AM_V_GEN)$(LN_S) -f $< $@

example_name  = adaptivity_ex5
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adaptivity/ex5
data          = adaptivity_ex5.C exact_solution.C run.sh
sources       = $(data) run.sh

CLEANFILES = solution_read_in.gmv solution_read_in.e saved_solution.xdr saved_solution.gmv saved_solution.e saved_mesh.xdr out_*.e out.gmv.*

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
