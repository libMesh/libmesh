# $Id: Makefile,v 1.2 2003-02-07 18:07:44 ddreyer Exp $


# The location of the mesh library
meshdir := ../..

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(meshdir)/Make.common


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ex6


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects) $(meshdir)/include/mesh_config.h
	@echo "Linking "$@"..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $@ $(LIBS) $(LDFLAGS)


# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) ifems_added.gmv out.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

run:
	@$(target) 
###############################################################################
