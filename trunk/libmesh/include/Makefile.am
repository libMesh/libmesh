SUBDIRS = libmesh


# where to install our headers
includedir=$(prefix)/include/libmesh

# definition of include_HEADERS - get from auto-maintained list
include include_HEADERS

# TESTS_ENVIRONMENT = CXX="$(CXX)" INCLUDEDIR="$(srcdir)" HEADERS_TO_TEST="$(nobase_include_HEADERS)" LIBMESH_DIR=$(abs_top_srcdir)
# TESTS = $(check_SCRIPTS)

# include our header checking script when doing 'make dist'
EXTRA_DIST = ../contrib/bin/test_installed_headers.sh

######################################################################
#
# What to do for 'make installcheck'
bindir=@bindir@/@host@_@METHOD@
libdir=@libdir@/@host@_@METHOD@
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
	@echo "Checking for standalone headers in source tree..."
	@echo " "
	@CXX="$(CXX)" \
          HEADERS_TO_TEST="$(headers_to_test)" \
          test_CXXFLAGS="$(libmesh_contrib_INCLUDES) $(libmesh_optional_INCLUDES) -I$(top_builddir)/include" \
          $(top_srcdir)/contrib/bin/test_installed_headers.sh

installcheck-local:
	@echo " "
	@echo "Checking for standalone headers in installed tree..."
	@echo " "
	@CXX="$(CXX)" \
          HEADERS_TO_TEST="$(headers_to_test)" \
          PKG_CONFIG=$(PKG_CONFIG) \
          PKG_CONFIG_PATH=$(pkgconfigdir) \
          LIBMESH_CONFIG_PATH=$(bindir) \
          $(top_srcdir)/contrib/bin/test_installed_headers.sh
