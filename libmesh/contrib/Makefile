#
# How to make the dependencies
#
all:
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves $(MAKECMDGOALS)



#
# How to clean the packages
#
clean:
	@rm -f *~
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves $(MAKECMDGOALS)

#
# How to clobber the packages
#
clobber:
	@$(MAKE) clean
	@$(MAKE) -C metis/Lib $(MAKECMDGOALS)
	@$(MAKE) -C sfcurves $(MAKECMDGOALS)
