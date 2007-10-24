# $Id: Makefile,v 1.35 2004-02-08 20:36:57 benkirk Exp $
#
# This is the Makefile for the libMesh library and helper
# applications.  This file is specific to the project.
# See the file Make.common for path configurations.
# Make.common is what should be included by applications
# using the library.


# include the library options determined by configure
include Make.common


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# header files & directories
headerfiles 	:= $(wildcard include/*/*.h)

#
# source files
srcfiles 	:= $(wildcard src/*/*.C)

#
# examples source files
examplesrcfiles	:= $(wildcard examples/ex*/ex*.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

#
# logged files -- all files you might want to log information for
loggedfiles     := $(srcfiles) $(filter-out include/base/libmesh_config.h, $(headerfiles))
###############################################################################



###############################################################################
# Target:
#
target := $(mesh_library)

ifeq ($(syn-mode),on)
  all:: $(objects)
else
  all:: $(target)
endif
###############################################################################

#
# 
#
.PHONY: echo echo_cxxflags echo_include echo_ldflags \
	clean clobber distclean \
	doc upload doc_upload log cvsweb TODO

#
# static library
#
$(mesh_library_dir)/libmesh.a: $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@$(AR) rv $(mesh_library) $(objects)
	@$(MAKE) -C contrib

#
# shared library
#
$(mesh_library_dir)/libmesh.so: $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $(mesh_library) $(objects) $(LDFLAGS)
	@$(MAKE) -C contrib

#
# Build the examples
#
examples: $(mesh_library) $(examplesrcfiles)
	@$(MAKE) -C examples

#
# Run the examples
#
run_examples: $(mesh_library)
	@$(MAKE) -C examples run


#
# Useful for checking make rules
#
echo:
	@echo -e "Source Files:\n$(srcfiles)\n"
	@echo -e "Object Files:\n$(objects)\n"
	@echo -e "Target:\n$(target)\n"
	@echo -e "Examples Source Files:\n$(examplesrcfiles)\n"
	@echo -e "CFLAGS:\n$(CXXFLAGS)\n"
	@echo -e "CXXFLAGS:\n$(CXXFLAGS)\n"
	@echo -e "INCLUDE:\n$(INCLUDE)\n"
	@echo -e "LIBS:\n$(LIBS)\n"
	@echo -e "LDFLAGS:\n$(LDFLAGS)\n"
	@echo -e "DLFLAGS:\n$(DLFLAGS)\n"

#
# Print the flags used for C++ compilation, padded with whitespace
#
echo_cxxflags:
	@echo -n " " $(CXXFLAGS) " "

#
# Print C++ compiler include path, padded with whitespace
#
echo_include:
	@echo -n " " $(INCLUDE) " "

#
# Print the flags used to link, padded with whitespace
#
echo_ldflags:
	@echo -n " " $(LIBS) $(LDFLAGS) $(DLFLAGS) " "

#	
# Remove project object files for the current mode
#
clean:
	@rm -f *~ include/*~ include/*/*~ src/*/*~ src/*/*.$(obj-suffix) doc/html/*~

#
# Make clean, remove contributed objects, and remove binaries
#
clobber:
	@$(MAKE) clean
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@$(MAKE) -C examples $(MAKECMDGOALS)
	@rm -rf config.status $(targ_dir) bin/grid2grid bin/meshtool bin/testexodus bin/amr bin/compare

#
# Make clobber, remove documentation, removes all libraries
# Should restore to a pristine state, except for files you
# have added
distclean:
	@$(MAKE) clobber
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@$(MAKE) -C examples $(MAKECMDGOALS)
	@rm -rf doc/html/*.html doc/html/*.png doc/html/*.gif \
                doc/html/doxygen/*.html doc/latex \
                doc/man/man3 doc/cvshtml/*.html doc/cvshtml/diff \
	        src/*/*.o src/*/*.g.o src/*/*.pg.o \
	        lib/*_opt lib/*_dbg lib/*_pro



#
# doxygen documentation
#
doc:
	$(doxygen) ./doc/Doxyfile

#
# Upload the web page to sourceforge
#
upload:
	chmod -R g+w ./doc/html/*
	rsync -rltzve ssh ./doc/html/ $(shell cat CVS/Root | cut -d"@" -f1 | cut -d":" -f3)@libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs
	chmod -R g-w ./doc/html/*

#
# Build and upload the documentation to sourceforge
#
doc_upload:
	@$(MAKE) doc
	@$(MAKE) upload


log: $(loggedfiles)
	cvs log $(loggedfiles) > cvs_log

#
# Web-based CVS documentation
# Please note: SourceForge provides a nice web interface to the CVS
# logs for the library.  You probably shouldn't need this target for
# anything.
cvsweb:
	./contrib/bin/cvs2html -f -p -o doc/cvshtml/index.html -v -a -b -n 2 -C crono.html

#
# Standalone applications.  Anything in the ./src/apps directory that ends in .cc
# can be compiled with this rule.  For example, if ./src/apps/foo.cc contains a main()
# and is a standalone program, then make bin/foo will work.
#
bin/% : src/apps/%.cc $(mesh_library)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ $(LIBS) $(LDFLAGS) $(DLFLAGS)

#
# In the contrib/bin directory, we run the test_headers.sh shell
# script.  This is a make rule for those tests.
#
contrib/bin/%.o : contrib/bin/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

#
# Make a TODO list
#
TODO:
	@egrep -i '// *todo' $(srcfiles) $(includes) \
	| perl -pi -e 's#(.*)(TODO:?)(\[.+\])#\3 \1\2\3#i;' \
	| perl -pi -e 's#\s*//\s*TODO:\s*(\[.+\])\s*#\n\1     #i;'



#
# Dependencies
#
.depend:
	@$(perl) ./contrib/bin/make_dependencies.pl $(foreach i, $(wildcard include/*), -I./$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"



#
# Make.common target.  Tell the user to run configue.
#
Make.common:
	@echo -e "You must run ./configure first!"
	exit 1

###############################################################################
include .depend


# Local Variables:
# mode: makefile
# End:
