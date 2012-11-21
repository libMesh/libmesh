# need to link lshaped.xda
BUILT_SOURCES = lshaped.xda
lshaped.xda:
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/examples/adjoints/adjoints_ex1/lshaped.xda .

check_PROGRAMS            = miscellaneous_ex3
miscellaneous_ex3dir      = $(examples_install_path)/miscellaneous/ex3
miscellaneous_ex3_DATA    = miscellaneous_ex3.C ../../adjoints/adjoints_ex1/lshaped.xda run.sh $(top_srcdir)/contrib/utils/Makefile
miscellaneous_ex3_SOURCES = $(miscellaneous_ex3_DATA) run.sh
miscellaneous_ex3_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = lshaped.xda out.e 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
