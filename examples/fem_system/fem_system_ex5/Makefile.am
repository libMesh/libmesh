# need to link sources from other example
BUILT_SOURCES = elasticity_system.C
DISTCLEANFILES = $(BUILT_SOURCES)

CLEANFILES = $(BUILT_SOURCES) out.e* out.newmark.* out.steady.*

example_name  = fem_system_ex5
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/fem_system/ex5
data          = fem_system_ex5.C fem_system_ex5.in ../fem_system_ex3/elasticity_system.C elasticity_system.h run.sh
sources       = $(data) run.sh

elasticity_system.C : $(top_srcdir)/examples/fem_system/fem_system_ex3/elasticity_system.C
	$(AM_V_GEN)$(LN_S) -f $< $@

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f fem_system_ex5.in && $(LN_S) -f $(srcdir)/fem_system_ex5.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += fem_system_ex5.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
