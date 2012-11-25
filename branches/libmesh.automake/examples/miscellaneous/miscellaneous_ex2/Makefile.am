CLEANFILES   = # start empty, append beloe

check_PROGRAMS            = miscellaneous_ex2
miscellaneous_ex2dir      = $(examples_install_path)/miscellaneous/ex2
miscellaneous_ex2_DATA    = miscellaneous_ex2.C lshape.unv lshape_data.unv
miscellaneous_ex2_SOURCES = $(miscellaneous_ex2_DATA) run.sh
miscellaneous_ex2_LDADD   = $(top_builddir)/libmesh.la

# conditionally enable this test
if LIBMESH_ENABLE_COMPLEX
  check_SCRIPTS = run.sh
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
