#
# How to make the dependencies
#
all:
	@echo "----------------------------------------------"
	@echo "------- Building Contributed Packages --------"
	@echo "----------------------------------------------"
	@$(MAKE) -C laspack   $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves  $(MAKECMDGOALS)
	@echo "----------------------------------------------"
	@echo "----- Done Building Contributed Packages -----"
	@echo "----------------------------------------------"



#
# How to clean the packages
#
clean:
	@rm -f *~
	@$(MAKE) -C laspack   $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves  $(MAKECMDGOALS)

#
# How to clobber the packages
#
clobber:
	@$(MAKE) clean
	@$(MAKE) -C laspack   $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves  $(MAKECMDGOALS)

#
# How to distclean the packages
#
distclean:
	@$(MAKE) clobber
	@$(MAKE) -C laspack   $(MAKECMDGOALS)
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves  $(MAKECMDGOALS)
	@rm -rf lib/*_dbg lib/*_opt lib/*_pro
