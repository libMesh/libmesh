# $Id: Makefile,v 1.1 2003-01-30 19:15:23 benkirk Exp $


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



###############################################################################
# Configuration: Augment libraries
#
ifeq ($(enable-shared), yes)
  LIBS := $(meshdir)/lib/$(hosttype)_$(METHOD)/libmesh.so $(LIBS)
else
  LIBS := $(meshdir)/lib/$(hosttype)_$(METHOD)/libmesh.a $(LIBS)
endif
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ex1


all:: $(target)

# Production rules:  how to make the target
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $@ $(LIBS)


# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target)

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.go *.pgo
###############################################################################
