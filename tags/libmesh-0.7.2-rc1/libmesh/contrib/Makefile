#
# How to make the dependencies
#
all:
	@echo "----------------------------------------------"
	@echo "------- Building Contributed Packages --------"
	@echo "----------------------------------------------"
	@echo ""
	@echo "--- Building LASPACK -------------------------"
	@$(MAKE) -C laspack      $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building Metis ---------------------------"
	@$(MAKE) -C metis/Lib    $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building Parmetis ------------------------"
	@$(MAKE) -C parmetis/Lib $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building sfcurves ------------------------"
	@$(MAKE) -C sfcurves     $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building libgzstream ---------------------"
	@$(MAKE) -C gzstream     $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building TetGen ---------------------"
	@$(MAKE) -C tetgen       $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building triangle ---------------------"
	@$(MAKE) -C triangle     $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building gmv ---------------------"
	@$(MAKE) -C gmv          $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building netCDF ---------------------"
	@$(MAKE) -C netcdf/Lib   $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building ExodusII ---------------------"
	@$(MAKE) -C exodusii/Lib/src $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building Nemesis ---------------------"
	@$(MAKE) -C nemesis/Lib $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building libHilbert ---------------------"
	@$(MAKE) -C libHilbert $(MAKECMDGOALS)
	@echo ""
	@echo "----------------------------------------------"
	@echo "----- Done Building Contributed Packages -----"
	@echo "----------------------------------------------"



#
# How to clean the packages
#
clean:
	@rm -f *~
	@$(MAKE) -C laspack      $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib    $(MAKECMDGOALS)
	@$(MAKE) -C parmetis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves     $(MAKECMDGOALS)
	@$(MAKE) -C gzstream     $(MAKECMDGOALS)
	@$(MAKE) -C tetgen       $(MAKECMDGOALS)
	@$(MAKE) -C triangle     $(MAKECMDGOALS)
	@$(MAKE) -C gmv          $(MAKECMDGOALS)
	@$(MAKE) -C netcdf/Lib   $(MAKECMDGOALS)
	@$(MAKE) -C exodusii/Lib/src $(MAKECMDGOALS)
	@$(MAKE) -C nemesis/Lib  $(MAKECMDGOALS)
	@$(MAKE) -C libHilbert   $(MAKECMDGOALS)

#
# How to clobber the packages
#
clobber:
	@$(MAKE) clean
	@$(MAKE) -C laspack      $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib    $(MAKECMDGOALS)
	@$(MAKE) -C parmetis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves     $(MAKECMDGOALS)
	@$(MAKE) -C gzstream     $(MAKECMDGOALS)
	@$(MAKE) -C tetgen       $(MAKECMDGOALS)
	@$(MAKE) -C triangle     $(MAKECMDGOALS)
	@$(MAKE) -C gmv          $(MAKECMDGOALS)
	@$(MAKE) -C netcdf/Lib   $(MAKECMDGOALS)
	@$(MAKE) -C exodusii/Lib/src $(MAKECMDGOALS)
	@$(MAKE) -C nemesis/Lib   $(MAKECMDGOALS)
	@$(MAKE) -C libHilbert   $(MAKECMDGOALS)

#
# How to distclean the packages
#
distclean:
	@$(MAKE) clobber
	@$(MAKE) -C laspack      $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib    $(MAKECMDGOALS)
	@$(MAKE) -C parmetis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves     $(MAKECMDGOALS)
	@$(MAKE) -C gzstream     $(MAKECMDGOALS)
	@$(MAKE) -C tetgen       $(MAKECMDGOALS)
	@$(MAKE) -C triangle     $(MAKECMDGOALS)
	@$(MAKE) -C gmv          $(MAKECMDGOALS)
	@$(MAKE) -C netcdf/Lib   $(MAKECMDGOALS)
	@$(MAKE) -C exodusii/Lib/src $(MAKECMDGOALS)
	@$(MAKE) -C nemesis/Lib   $(MAKECMDGOALS)
	@$(MAKE) -C libHilbert   $(MAKECMDGOALS)
	@rm -rf lib/*_dbg lib/*_opt lib/*_prof lib/*_oprof lib/*_devel
