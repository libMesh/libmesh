all:: examples

#
# Build the example programs
examples:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)


.PHONY: clean clobber distclean

#
# Clean up the example programs
clean:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)


#
# Really clean up the example programs
clobber:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)

#
# Make the example program directories look
# like a clean distribution
distclean:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)


