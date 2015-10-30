example_name  = miscellaneous_ex8
install_dir   = $(examples_install_path)/miscellaneous/ex8
data          = miscellaneous_ex8.C meshless_interpolation_function.h run.sh unstruct.ucd.gz struct.ucd.gz
sources       = $(data) run.sh

# conditionally enable this test
if LIBMESH_ENABLE_NANOFLANN
  check_SCRIPTS = run.sh
endif

CLEANFILES = src.dat dest*.dat

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f struct.ucd.gz   && $(LN_S) -f $(srcdir)/struct.ucd.gz
	-rm -f unstruct.ucd.gz && $(LN_S) -f $(srcdir)/unstruct.ucd.gz
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += struct.ucd.gz unstruct.ucd.gz .linkstamp
endif


##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
