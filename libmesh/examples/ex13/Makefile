# $Id$


# The location of the mesh library
LIBMESH_DIR ?= ../..

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
target 	   :=  ./ex13-$(METHOD)
#target 	   :=  ex13 ex13_nonlinear


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
%-$(METHOD): $(patsubst %.C, %.$(obj-suffix), %.C)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $< -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)



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
	@echo "* Running Example " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@#./ex13_nonlinear -snes_ls cubic -snes_max_it 10 -ksp_rtol 1.e-3 -ksp_max_it 50 -ksp_xmonitor -ksp_gmres_restart 25 -snes_mf_operator -snes_view 
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
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
