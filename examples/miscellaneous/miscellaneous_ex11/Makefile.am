example_name  = miscellaneous_ex11
install_dir   = $(examples_install_path)/miscellaneous/ex11
data          = miscellaneous_ex11.C run.sh square_mesh.off
sources       = $(data)
check_SCRIPTS = run.sh

CLEANFILES = without_ghosts.pvtu without_ghosts_*.vtu without_ghosts.e with_ghosts.pvtu with_ghosts_*.vtu with_ghosts.e out.pvtu out_*.vtu out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f square_mesh.off && $(LN_S) -f $(srcdir)/square_mesh.off .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += square_mesh.off .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
