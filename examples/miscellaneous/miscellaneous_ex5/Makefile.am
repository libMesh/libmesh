# need to link lshaped.xda
BUILT_SOURCES = lshaped.xda
lshaped.xda:
	$(AM_V_GEN)$(LN_S) -f $(top_srcdir)/examples/adjoints/adjoints_ex1/lshaped.xda .

example_name  = miscellaneous_ex5
install_dir   = $(examples_install_path)/miscellaneous/ex5
data          = miscellaneous_ex5.C ../../adjoints/adjoints_ex1/lshaped.xda lshaped3D.xda miscellaneous_ex5.in run.sh
sources       = $(data) run.sh
check_SCRIPTS = run.sh

CLEANFILES = lshaped.xda lshaped_dg.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f miscellaneous_ex5.in && $(LN_S) -f $(srcdir)/miscellaneous_ex5.in .
	-rm -f lshaped3D.xda        && $(LN_S) -f $(srcdir)/lshaped3D.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += miscellaneous_ex5.in lshaped3D.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
