SUBDIRS = libmesh

# special handholding for prefix_config.m4 generated files
# so that 'make clean ; make' works as does 'make distcheck'
# libmesh_config.h is made by ./configure, so it should get
# cleaned by 'make distclean'
DISTCLEANFILES = libmesh_config.h

#
# special things to do when running 'make dist'
dist-hook:
	rm -rf $(distdir)/libmesh_config.h

# where to install our headers
includedir=$(prefix)/include/libmesh

# definition of include_HEADERS - get from auto-maintained list
include include_HEADERS

# include our header checking script when doing 'make dist'
EXTRA_DIST = ../contrib/bin/test_installed_headers.sh

######################################################################
#
# What to do for 'make installcheck'
#bindir=@bindir@/@host@_@METHOD@
#libdir=@libdir@/@host@_@METHOD@
pkgconfigdir=$(libdir)/pkgconfig

# the actual header files we want to test for self-sufficiency
headers_to_test = $(include_HEADERS)
if LIBMESH_ENABLE_EXODUS
  headers_to_test += contrib/exodusII.h
endif
if LIBMESH_ENABLE_FPARSER
  headers_to_test += contrib/fparser.hh
endif

test_headers:
	@echo " "
	@echo "Checking for standalone headers in source tree $(top_builddir)/include ..."
	@echo " "
	@CXX="$(CXX)" \
	  HEADERS_TO_TEST="$(headers_to_test)" \
	  test_CXXFLAGS="$(libmesh_contrib_INCLUDES) $(libmesh_optional_INCLUDES) -I$(top_builddir)/include" \
	  $(top_srcdir)/contrib/bin/test_installed_headers.sh

installcheck-local:
	@echo " "
	@echo "Checking for standalone headers in installed tree $(prefix) ..."
	@echo " "
	@CXX="$(CXX)" \
	  HEADERS_TO_TEST="$(headers_to_test)" \
	  PKG_CONFIG=$(PKG_CONFIG) \
	  PKG_CONFIG_PATH=$(pkgconfigdir) \
	  LIBMESH_CONFIG_PATH=$(bindir) \
	  $(top_srcdir)/contrib/bin/test_installed_headers.sh
