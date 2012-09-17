# The location of the mesh library
LIBMESH_DIR ?= ../../..

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
target 	   :=  ./systems_of_equations_ex4-$(METHOD)


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
%-$(METHOD): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)



# Useful rules.
clean:
	@rm -f $(objects) *.gmv.* *~ .depend

clobber:
	@$(MAKE) clean
	@rm -f $(target) out*.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o .depend

run: $(target)
	@echo "***************************************************************"
	@echo "* Running " $(LIBMESH_RUN) $(target) -ksp_type cg #$(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(target) -ksp_type cg #$(LIBMESH_OPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running " $(LIBMESH_RUN) $(target) -ksp_type cg #$(LIBMESH_OPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend: $(srcfiles) $(LIBMESH_DIR)/include/*/*.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(LIBMESH_DIR)/include $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################
