includedir = $(prefix)/include/libmesh

if LIBMESH_ENABLE_CAPNPROTO

# Description of the BUILT_SOURCES variable:
# A source file listed in BUILT_SOURCES is made on 'make all' or 'make
# check' (or even 'make install') before other targets are
# processed. However, such a source file is not compiled unless
# explicitly requested by mentioning it in some other _SOURCES
# variable.
BUILT_SOURCES = $(builddir)/rb_data.capnp.c++ $(builddir)/rb_data.capnp.h

# We need to specify the source file's dependencies and how it is
# built. The capnp program generates both the header and source file,
# so we don't need to have a separate target for each.
#
# Notes:
# 1.) We use double-$ in the perl commands below to escape the single-$ sign from Make.
# 2.) We have already tested that capnp is in the user's PATH (this is required for capnp to work)
#     so there is no need to provide a PATH here.
$(builddir)/rb_data.capnp.h: $(srcdir)/rb_data.capnp
	capnp compile --src-prefix=$(srcdir) -oc++:$(builddir) $(srcdir)/rb_data.capnp
	perl -pi -e 'print "#include \"libmesh/libmesh_config.h\"\n#if defined(LIBMESH_HAVE_CAPNPROTO)\n" if $$. == 1' $(builddir)/rb_data.capnp.h
	perl -pi -e 'eof && s/$$/\n#endif \/\/ #if defined(LIBMESH_HAVE_CAPNPROTO)/;' $(builddir)/rb_data.capnp.h

# The rule above actually generates both the .h and .c++ files.  This
# is just a dummy rule so that make knows how the c++ file is
# "generated", and makes the whole thing work with "make -j" runs.
$(builddir)/rb_data.capnp.c++: $(builddir)/rb_data.capnp.h

pkg_cppflags = $(CAPNPROTO_INCLUDE)

endif

# placing this here, unconditionally based on LIBMESH_ENABLE_CAPNPROTO,
# ensures that it will be bundled up by 'make dist' & 'make distcheck'
pkg_sources = rb_data.capnp

# Make.common is another automake file that handles compiling/linking flags and METHODs stuff.
include $(top_srcdir)/contrib/Make.common
