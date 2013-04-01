BUILT_SOURCES = mesh.xda

mesh.xda: $(top_srcdir)/examples/adaptivity/adaptivity_ex2/mesh.xda
	$(AM_V_GEN)$(LN_S) $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

example_name  = subdomains_ex3
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/subdomains/ex3
data          = subdomains_ex3.C ../../adaptivity/adaptivity_ex2/mesh.xda hybrid_3d.xda run.sh
sources       = $(data) run.sh

CLEANFILES = in_*.dat out_*.dat

# also need to link hybrid_3d.xda for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f hybrid_3d.xda && $(LN_S) $(srcdir)/hybrid_3d.xda .
	$(AM_V_GEN)touch .linkstamp

  DISTCLEANFILES += hybrid_3d.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
