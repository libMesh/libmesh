BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/transient/transient_ex1/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS         = adaptivity_ex5
check_SCRIPTS          = run.sh
adaptivity_ex5dir      = $(examples_install_path)/adaptivity/ex5
adaptivity_ex5_DATA    = adaptivity_ex5.C exact_solution.C run.sh $(top_srcdir)/contrib/utils/Makefile
adaptivity_ex5_SOURCES = $(adaptivity_ex5_DATA) run.sh
adaptivity_ex5_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = solution_read_in.gmv solution_read_in.e saved_solution.xdr saved_solution.gmv saved_solution.e saved_mesh.xdr out_*.e out.gmv.*

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
