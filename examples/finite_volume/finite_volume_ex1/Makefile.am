example_name  = finite_volume_ex1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/finite_volume/ex1
data          = finite_volume_ex1.C finite_volume_ex1.in run.sh
sources       = $(data)

CLEANFILES = out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f finite_volume_ex1.in && $(LN_S) -f $(srcdir)/finite_volume_ex1.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += finite_volume_ex1.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
