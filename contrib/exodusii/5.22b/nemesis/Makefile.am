pkg_cppflags = $(EXODUS_INCLUDE) $(NETCDF_INCLUDE)

pkg_sources = \
	ne_nemesisI.h \
	ne_nemesisI_int.h \
	nemesis_wrapper.c

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common

#####################################################
## support the ne_ctest_wrap function for 'make check'
#check_PROGRAMS = 
#
#if LIBMESH_DBG_MODE
#  check_PROGRAMS            += ne_ctest_wrap-dbg
#  ne_ctest_wrap_dbg_SOURCES  = ne_ctest_wrap.c
#  ne_ctest_wrap_dbg_LDADD    = libdbg.la \
#			       $(top_builddir)/contrib/exodusii/v5.22/exodus/cbind/libdbg.la \
#                               $(top_builddir)/contrib/netcdf/v4/liblib/libnetcdf.la
#endif
#
#if LIBMESH_DEVEL_MODE
#  check_PROGRAMS              += ne_ctest_wrap-devel
#  ne_ctest_wrap_devel_SOURCES  = ne_ctest_wrap.c
#  ne_ctest_wrap_devel_LDADD    = libdevel.la \
#			         $(top_builddir)/contrib/exodusii/v5.22/exodus/cbind/libdevel.la \
#                                 $(top_builddir)/contrib/netcdf/v4/liblib/libnetcdf.la
#endif
#
#if LIBMESH_OPT_MODE
#  check_PROGRAMS            += ne_ctest_wrap-opt
#  ne_ctest_wrap_opt_SOURCES  = ne_ctest_wrap.c
#  ne_ctest_wrap_opt_LDADD    = libopt.la \
#			       $(top_builddir)/contrib/exodusii/v5.22/exodus/cbind/libopt.la \
#                               $(top_builddir)/contrib/netcdf/v4/liblib/libnetcdf.la
#endif
#
#TESTS = $(check_PROGRAMS)
#CLEANFILES = ne_test.exoII

