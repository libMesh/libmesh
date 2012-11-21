BUILT_SOURCES = exact_solution.C
exact_solution.C: $(top_srcdir)/examples/transient/transient_ex1/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

check_PROGRAMS         = adaptivity_ex2
check_SCRIPTS          = run.sh
adaptivity_ex2dir      = $(examples_install_path)/adaptivity/ex2
adaptivity_ex2_DATA    = adaptivity_ex2.C exact_solution.C mesh.xda
adaptivity_ex2_SOURCES = $(adaptivity_ex2_DATA) run.sh
adaptivity_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = solution_read_in.gmv saved_solution.xda saved_solution.gmv saved_mesh.xda out.gmv.*

# also need to link mesh.xda for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f mesh.xda && $(LN_S) $(srcdir)/mesh.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += mesh.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
