CLEANFILES   = # start empty, append below

example_name = miscellaneous_ex2
install_dir  = $(examples_install_path)/miscellaneous/ex2
data         = miscellaneous_ex2.C lshape.unv lshape_data.unv
sources      = $(data) run.sh

# conditionally enable this test
if LIBMESH_ENABLE_COMPLEX
  check_SCRIPTS = run.sh
endif

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f lshape.unv && $(LN_S) -f $(srcdir)/lshape.unv .
	-rm -f lshape_data.unv && $(LN_S) -f $(srcdir)/lshape_data.unv .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += lshape.unv lshape_data.unv out_*.e eqn_sys.dat .linkstamp
endif


##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
