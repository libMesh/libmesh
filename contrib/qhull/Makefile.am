if LIBMESH_ENABLE_QHULL

# Use a stamp file to make sure we don't have multiple sub-$(MAKE)
# going at once; cmake gets horribly confused if that happens.
.buildstamp:
	$(MAKE) -C build $(AM_MAKEFLAGS)
	touch .buildstamp

all-local: .buildstamp

clean-local:
	$(MAKE) -C build $(AM_MAKEFLAGS) clean
	rm -f .buildstamp

distclean-local:
	rm -rf build install
	rm -f .buildstamp

qhulllibdir = $(libdir)
qhulllib_DATA = build/libqhullstatic_r.a \
                build/libqhullcpp.a

build/libqhullstatic_r.a: .buildstamp
build/libqhullcpp.a: .buildstamp

endif

EXTRA_DIST = qhull
