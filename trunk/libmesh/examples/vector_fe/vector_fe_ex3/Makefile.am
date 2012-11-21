check_PROGRAMS        = vector_fe_ex3
#check_SCRIPTS         = run.sh
vector_fe_ex3dir      = $(examples_install_path)/vector_fe/ex3
vector_fe_ex3_DATA    = curl_curl_exact_solution.h vector_fe_ex3.in curl_curl_system.C solution_function.h curl_curl_system.h vector_fe_ex3.C
vector_fe_ex3_SOURCES = $(vector_fe_ex3_DATA) run.sh
vector_fe_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = 

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex3.in && $(LN_S) $(srcdir)/vector_fe_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
