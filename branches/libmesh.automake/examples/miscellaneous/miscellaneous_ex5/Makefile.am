# need to link lshaped.xda
BUILT_SOURCES = lshaped.xda
lshaped.xda:
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/examples/adjoints/adjoints_ex1/lshaped.xda .

check_PROGRAMS            = miscellaneous_ex5
miscellaneous_ex5dir      = $(examples_install_path)/miscellaneous/ex5
miscellaneous_ex5_DATA    = miscellaneous_ex5.C ../../adjoints/adjoints_ex1/lshaped.xda lshaped3D.xda miscellaneous_ex5.in run.sh $(top_srcdir)/contrib/utils/Makefile
miscellaneous_ex5_SOURCES = $(miscellaneous_ex5_DATA) run.sh
miscellaneous_ex5_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = lshaped.xda lshaped_dg.e 

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f miscellaneous_ex5.in && $(LN_S) $(srcdir)/miscellaneous_ex5.in .
	-rm -f lshaped3D.xda        && $(LN_S) $(srcdir)/lshaped3D.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += miscellaneous_ex5.in lshaped3D.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
