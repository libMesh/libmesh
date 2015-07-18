includedir = $(prefix)/include/libmesh


# these "sources" are built from rb_data.capnp, and also will be cleaned
# note we use rb_data-capnp.c++ rather than rb_data.capnp.c++ to avoid
# an erroneous suffix rule conflict with rb_data.capnp
BUILT_SOURCES = rb_data.capnp.h rb_data-capnp.c++
CLEANFILES = $(BUILT_SOURCES)

# these are the sources that go into building the library, and also rolled up
# into the distribution
pkg_sources = $(BUILT_SOURCES) rb_data.capnp

if LIBMESH_ENABLE_CAPNPROTO

# We need to specify the source file's dependencies and how it is
# built. The capnp program generates both the header and source file,
# so we don't need to have a separate target for each.
#
# Notes:
# 1.) We use double-$ in the perl commands below to escape the single-$ sign from Make.
# 2.) We have already tested that capnp is in the user's PATH (this is required for capnp to work)
#     so there is no need to provide a PATH here.
rb_data.capnp.h: rb_data.capnp
	$(AM_V_GEN)capnp compile --src-prefix=$(srcdir) -oc++:$(builddir) $(srcdir)/rb_data.capnp
	perl -pi -e 'print "#include \"libmesh/libmesh_config.h\"\n#if defined(LIBMESH_HAVE_CAPNPROTO)\n" if $$. == 1' rb_data.capnp.h
	perl -pi -e 'eof && s/$$/\n#endif \/\/ #if defined(LIBMESH_HAVE_CAPNPROTO)/;' rb_data.capnp.h
	mv rb_data.capnp.c++ rb_data-capnp.c++

# The rule above actually generates both the .h and .c++ files.  This
# is just a dummy rule so that make knows how the c++ file is
# "generated", and makes the whole thing work with "make -j" runs.
rb_data-capnp.c++: rb_data.capnp.h

pkg_cppflags = $(CAPNPROTO_INCLUDE)

# These are the headers we actually want to install and make available
# to a user.  Consider these the public interface of libMesh.
# will get installed into '$(prefix)/include/libmesh'
include_HEADERS = rb_data.capnp.h

else !LIBMESH_ENABLE_CAPNPROTO
rb_data.capnp.h:
	$(AM_V_GEN)echo "// nothing needed without LIBMESH_ENABLE_CAPNPROTO" > $@

rb_data-capnp.c++:
	$(AM_V_GEN)echo "// nothing needed without LIBMESH_ENABLE_CAPNPROTO" > $@
endif

# Make.common is another automake file that handles compiling/linking flags and METHODs stuff.
include $(top_srcdir)/contrib/Make.common

# Ensure these don't get packaged up, with or without Cap'N Proto we want to make sure thet are rebuilt
dist-hook:
	rm -f $(distdir)/rb_data.capnp.h \
	      $(distdir)/rb_data-capnp.c++
