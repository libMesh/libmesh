# $Id: Makefile,v 1.1.1.1 2009/03/19 20:28:07 dave Exp $

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common

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
target 	   :=  ./ex23-$(METHOD)


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
%-$(METHOD): $(objects)
	@echo "Linking "$@"..."
	$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)



# Useful rules.
clean:
	@rm -f $(objects) *.gmv.* *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) out*.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

run: $(target)
	@echo "***************************************************************"
	@echo "* Running " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
