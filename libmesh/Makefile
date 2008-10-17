# $Id$
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
examplesrcfiles	:= $(wildcard examples/ex*/*.C)

#
# apps source files
appsrcfiles	:= $(wildcard src/apps/*.cc)

#
# apps binary files
appbinfiles	:= $(patsubst %.cc, %, $(addprefix bin/, $(notdir $(appsrcfiles))))

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
  all:: $(target) $(appbinfiles)
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
ifeq ($(findstring darwin,$(hostos)),darwin)
$(mesh_library_dir)/libmesh$(static_libext): $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@libtool -static -o $(mesh_library) $(objects)
	@$(MAKE) -C contrib

else
$(mesh_library_dir)/libmesh$(static_libext): $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@$(AR) rv $(mesh_library) $(objects)
	@$(MAKE) -C contrib
endif
# shared library
#
$(mesh_library_dir)/libmesh$(shared_libext): $(objects)
	@$(shell mkdir -p $(mesh_library_dir))
	@echo "Linking "$@
	@$(libmesh_CXX) $(libmesh_CXXSHAREDFLAG) -o $(mesh_library) $(objects) $(libmesh_LDFLAGS)
	@$(MAKE) -C contrib

#
# Build just object files
#
obj: $(objects)

#
# Build the examples
#
examples: $(mesh_library) $(examplesrcfiles)
	@$(MAKE) -C examples

#
# Only link the examples.  Useful on machines where you
# do not have permission to run programs outside of a queue.
#
link_examples: $(mesh_library)
	@$(MAKE) -C examples link


#
# Run the examples
#
run_examples: $(mesh_library)
	@$(MAKE) -C examples run


#
# Useful for checking make rules
#
echo:
	@echo -e "C++ compiler:\n$(libmesh_CXX)\n"
	@echo -e "Source Files:\n$(srcfiles)\n"
	@echo -e "Object Files:\n$(objects)\n"
	@echo -e "Target:\n$(target)\n"
	@echo -e "Examples Source Files:\n$(examplesrcfiles)\n"
	@echo -e "libmesh_CFLAGS:\n$(libmesh_CFLAGS)\n"
	@echo -e "libmesh_CXXFLAGS:\n$(libmesh_CXXFLAGS)\n"
	@echo -e "libmesh_CXXSHAREDFLAG:\n$(libmesh_CXXSHAREDFLAG)\n"
	@echo -e "libmesh_INCLUDE:\n$(libmesh_INCLUDE)\n"
	@echo -e "libmesh_LIBS:\n$(libmesh_LIBS)\n"
	@echo -e "libmesh_LDFLAGS:\n$(libmesh_LDFLAGS)\n"
	@echo -e "libmesh_DLFLAGS:\n$(libmesh_DLFLAGS)\n"
	@echo -e "EXAMPLES:\n$(examplesrcfiles)\n"

#
# Print the name of the C++ compiler, padded with whitespace
#
echo_cxx:
	@echo  " " $(libmesh_CXX) " "

#
# Print the flags used for C++ compilation, padded with whitespace
#
echo_cxxflags:
	@echo " " $(libmesh_CXXFLAGS) " "

#
# Print C++ compiler include path, padded with whitespace
#
echo_include:
	@echo " " $(libmesh_INCLUDE) " "

#
# Print the flags used to link, padded with whitespace
#
echo_ldflags:
	@echo " " $(libmesh_LIBS) $(libmesh_LDFLAGS) $(libmesh_DLFLAGS) " "

#	
# Remove object files for the current mode
#
clean:
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@$(MAKE) -C examples $(MAKECMDGOALS)
	@rm -f *~ include/*~ include/*/*~ src/*/*~ src/*/*.$(obj-suffix) doc/html/*~

#
# Make clean, remove all binaries and generated files.  Leaves libraries in-tact
#
clobber:
	@$(MAKE) clean
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@$(MAKE) -C examples $(MAKECMDGOALS)
	@rm -rf config.status $(targ_dir) $(appbinfiles)

#
# Make clobber, remove documentation, removes all libraries & object files for all modes
# Should restore to a pristine state, except for files you
# have added
distclean:
	@$(MAKE) clobber
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@$(MAKE) -C examples $(MAKECMDGOALS)
	@rm -rf doc/man/man3
	@rm -rf doc/html/doxygen/*.html # split these up, otherwise command line gets too long
	@rm -rf doc/html/doxygen/*.php
	@rm -rf doc/html/doxygen/*.png
	@rm -rf doc/html/doxygen/*.gif
	@rm -rf doc/html/doxygen/*.map
	@rm -rf doc/html/doxygen/*.md5
	@rm -rf doc/html/doxygen/*.dot
	@rm -rf doc/html/doxygen/formula.repository doc/html/doxygen/graph_legend.dot
	@rm -rf doc/latex/doxygen
	@rm -rf doc/latex/*/*.aux doc/latex/*/*~ doc/latex/*/*.log doc/latex/*/*.out
	@rm -rf src/*/*.o
	@rm -rf lib/*_opt lib/*_dbg lib/*_pro lib/*_devel



#
# doxygen documentation
#
doc:
	$(doxygen) ./doc/Doxyfile
	@rm -rf doc/html/doxygen/*.map
	@rm -rf doc/html/doxygen/*.md5
	@rm -rf doc/html/doxygen/*.dot
	@rm -rf doc/html/doxygen/formula.repository

#
# Upload the web page to sourceforge.  We need a way to specify usernames
# other than $USER when connecting to sourceforge servers... Please set the
# environment variable: $LIBMESH_SVN_USER if your sourceforge username
# is different than whatever $USER is.
#
ifeq (x$(LIBMESH_SVN_USER),x) #unset
  upload_name=$(USER),libmesh@
else
  upload_name=$(LIBMESH_SVN_USER),libmesh@
endif

upload:
	chmod -R g+w ./doc/html/* ./doc/latex/*/*
	rsync -rltzve ssh --exclude '.svn' ./doc/html/ $(upload_name)web.sourceforge.net:/home/groups/l/li/libmesh/htdocs
	rsync -rltzve ssh --exclude '.svn' ./doc/latex/howto $(upload_name)web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -rltzve ssh --exclude '.svn' ./doc/latex/xda_format $(upload_name)web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	chmod -R g-w ./doc/html/* ./doc/latex/*/*


#
# Test Uploading (rsync -n) the web page to sourceforge.  Try this
# if you want to see what *would be* uploaded with a real rsync.
#
upload_test:
	chmod -R g+w ./doc/html/* ./doc/latex/*/*
	rsync -nrltzve ssh --exclude '.svn' ./doc/html/ $(upload_name)libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs
	rsync -nrltzve ssh --exclude '.svn' ./doc/latex/howto $(upload_name)libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -nrltzve ssh --exclude '.svn' ./doc/latex/xda_format $(upload_name)libmesh.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	chmod -R g-w ./doc/html/* ./doc/latex/*/*


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
	@echo "Building $@"
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) $< -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) $(libmesh_DLFLAGS)

#
# In the contrib/bin directory, we run the test_headers.sh shell
# script.  This is a make rule for those tests.
#
contrib/bin/%.o : contrib/bin/%.cc
	$(libmesh_CXX) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $< -o $@

#
# Make a TODO list
#
TODO:
	@egrep -i '// *todo' $(srcfiles) $(headerfiles) \
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
