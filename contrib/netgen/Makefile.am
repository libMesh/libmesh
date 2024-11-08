if LIBMESH_ENABLE_NETGEN

# Use a stamp file to make sure we don't have multiple sub-$(MAKE)
# going at once; cmake gets horribly confused if that happens.
.buildstamp:
	$(MAKE) -C build $(AM_MAKEFLAGS)
	touch .buildstamp

all-local : .buildstamp

# This doesn't seem to be cleaning...
clean-local:
	$(MAKE) -C build $(AM_MAKEFLAGS) clean
	rm -f .buildstamp


netgenincludedir=$(includedir)/netgen
nglibincludedir=$(includedir)/netgen/nglib
netgenlibdir=$(libdir)
ngliblibdir=$(libdir)

# Make sure the installed libnglib is linking to the installed
# libngcore, not to the libngcore in the build directory
#
# Install-exec is automake speak for anything platform-dependent, not
# just for executables, so it should probably be what we use here, but
# I'm not sure the best way to tell automake "this is
# platform-dependent but we're building it ourselves", so I'm calling
# these libraries _DATA for now.
#
# We manually add the install-ngliblibDATA dependency so patchelf etc.
# definitely are looking at an existing correct file when we do `make
# install` in parallel.
install-data-local: install-ngliblibDATA
	@if which install_name_tool 2>/dev/null; then \
	  echo "install_name_tool exists on OSX but is unneeded for rpath-based libnglib there"; \
	elif which patchelf 2>/dev/null; then \
	  patchelf --set-rpath "$(libdir)" "$(libdir)/libnglib$(LIBMESH_LIBRARY_SUFFIX)"; \
	elif which chrpath 2>/dev/null; then \
	  chrpath -r "$(libdir)" "$(libdir)/libnglib$(LIBMESH_LIBRARY_SUFFIX)"; \
	else \
	  echo "No tools to change RPATH when installing libnglib$(LIBMESH_LIBRARY_SUFFIX)"; \
	  echo "If linking fails we may need install_name_tool/patchelf/chrpath installed"; \
	fi

.PHONY : .buildstamp

.buildstamp:
	$(MAKE) -C build $(AM_MAKEFLAGS)
	touch .buildstamp

netgeninclude_HEADERS = build/netgen/netgen_config.hpp \
                  build/netgen/netgen_version.hpp

nglibinclude_HEADERS = netgen/nglib/nglib.h

BUILT_SOURCES = build/netgen/netgen_config.hpp \
                build/netgen/netgen_version.hpp

ngliblib_DATA = build/netgen/libnglib$(LIBMESH_LIBRARY_SUFFIX)
netgenlib_DATA = build/netgen/libsrc/core/libngcore$(LIBMESH_LIBRARY_SUFFIX)

build/netgen/netgen_config.hpp : .buildstamp
build/netgen/netgen_version.hpp : .buildstamp
build/netgen/libnglib$(LIBMESH_LIBRARY_SUFFIX): .buildstamp
build/netgen/libsrc/core/libngcore$(LIBMESH_LIBRARY_SUFFIX) : .buildstamp

endif

EXTRA_DIST = netgen

# Maybe useful later?
#              build/netgen/libnggui.so \
#              build/netgen/libsrc/core/pyngcore.cpython-310-x86_64-linux-gnu.so \
#              build/netgen/ng/libngpy.so \
#              build/netgen/ng/libngguipy.so
