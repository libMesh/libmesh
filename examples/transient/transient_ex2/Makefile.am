example_name  = transient_ex2
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/transient/ex2
data          = transient_ex2.C pipe-mesh.unv run.sh
sources       = $(data) run.sh

CLEANFILES = pressure_node.res out.*.gmv out_*.vtu  out_*.pvtu

if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f pipe-mesh.unv && $(LN_S) $(srcdir)/pipe-mesh.unv .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += pipe-mesh.unv .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
