check_PROGRAMS            = miscellaneous_ex6
miscellaneous_ex6dir      = $(examples_install_path)/miscellaneous/ex6
miscellaneous_ex6_DATA    = miscellaneous_ex6.C
miscellaneous_ex6_SOURCES = $(miscellaneous_ex6_DATA) run.sh
miscellaneous_ex6_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = hole_3D.e delaunay_l_shaped_hole.e 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
