check_PROGRAMS        = transient_ex2
check_SCRIPTS         = run.sh
transient_ex2dir      = $(examples_install_path)/transient/ex2
transient_ex2_DATA    = transient_ex2.C pipe-mesh.unv
transient_ex2_SOURCES = $(transient_ex2_DATA) run.sh
transient_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = pressure_node.res out.*.gmv

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
