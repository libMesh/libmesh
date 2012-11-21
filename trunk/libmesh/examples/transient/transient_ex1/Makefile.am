BUILT_SOURCES = mesh.xda

mesh.xda: $(top_srcdir)/examples/adaptivity/adaptivity_ex2/mesh.xda
	$(AM_V_GEN)$(LN_S) $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

check_PROGRAMS        = transient_ex1
check_SCRIPTS         = run.sh
transient_ex1dir      = $(examples_install_path)/transient/ex1
transient_ex1_DATA    = transient_ex1.C exact_solution.C ../../adaptivity/adaptivity_ex2/mesh.xda run.sh $(top_srcdir)/contrib/utils/Makefile
transient_ex1_SOURCES = $(transient_ex1_DATA) run.sh
transient_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.gmv

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
