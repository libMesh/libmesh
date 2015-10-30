BUILT_SOURCES = exact_solution.C
exact_solution.C: $(top_srcdir)/examples/transient/transient_ex1/exact_solution.C
	$(AM_V_GEN)$(LN_S) -f $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

example_name  = adaptivity_ex2
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adaptivity/ex2
data          = adaptivity_ex2.C exact_solution.C mesh.xda run.sh
sources       = $(data) run.sh

CLEANFILES = solution_read_in.gmv saved_solution.xda saved_solution.gmv saved_mesh.xda out.gmv.*

# also need to link mesh.xda for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f mesh.xda && $(LN_S) -f $(srcdir)/mesh.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += mesh.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
