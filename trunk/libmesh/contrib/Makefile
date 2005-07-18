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
	@$(MAKE) -C tetgen     $(MAKECMDGOALS)
	@echo ""
	@echo "--- Building triangle ---------------------"
	@$(MAKE) -C triangle     $(MAKECMDGOALS)
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
	@rm -rf lib/*_dbg lib/*_opt lib/*_pro lib/*_devel
