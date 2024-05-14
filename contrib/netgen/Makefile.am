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

.PHONY : .buildstamp

.buildstamp:
	$(MAKE) -C build $(AM_MAKEFLAGS)
	touch .buildstamp

netgeninclude_HEADERS = build/netgen/netgen_config.hpp \
                  build/netgen/netgen_version.hpp

nglibinclude_HEADERS = netgen/nglib/nglib.h

BUILT_SOURCES = build/netgen/netgen_config.hpp \
                build/netgen/netgen_version.hpp

build/netgen/netgen_config.hpp : .buildstamp
build/netgen/netgen_version.hpp : .buildstamp
endif

#netgendir = $(libdir)
#netgen_DATA = build/netgen/libnglib.so \
#              build/netgen/libsrc/core/libngcore.so

# Maybe useful later?
#              build/netgen/libnggui.so \
#              build/netgen/libsrc/core/pyngcore.cpython-310-x86_64-linux-gnu.so \
#              build/netgen/ng/libngpy.so \
#              build/netgen/ng/libngguipy.so

#netgen_SOURCES = build/netgen/netgen_config.hpp \
#                 build/netgen/netgen_version.hpp \
#                 netgen/nglib/nglib.h
