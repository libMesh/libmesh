all:: examples

#
# Build the example programs
examples:
	@$(MAKE) -C ex1
	@$(MAKE) -C ex2
	@$(MAKE) -C ex3
	@$(MAKE) -C ex4
	@$(MAKE) -C ex5
	@$(MAKE) -C ex7


.PHONY: clean clobber distclean

#
# Clean up the example programs
clean:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)
	@$(MAKE) -C ex3 $(MAKECMDGOALS)
	@$(MAKE) -C ex4 $(MAKECMDGOALS)
	@$(MAKE) -C ex5 $(MAKECMDGOALS)
	@$(MAKE) -C ex7 $(MAKECMDGOALS)


#
# Really clean up the example programs
clobber:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)
	@$(MAKE) -C ex3 $(MAKECMDGOALS)
	@$(MAKE) -C ex4 $(MAKECMDGOALS)
	@$(MAKE) -C ex5 $(MAKECMDGOALS)
	@$(MAKE) -C ex7 $(MAKECMDGOALS)

#
# Make the example program directories look
# like a clean distribution
distclean:
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)
	@$(MAKE) -C ex3 $(MAKECMDGOALS)
	@$(MAKE) -C ex4 $(MAKECMDGOALS)
	@$(MAKE) -C ex5 $(MAKECMDGOALS)
	@$(MAKE) -C ex7 $(MAKECMDGOALS)



run_print_line:
	@echo ""
	@echo "********************************************************************************"
	@echo "* FINISHED RUNNING EXAMPLE, PROCEED TO THE NEXT                                *"
	@echo "********************************************************************************"
	@echo ""


#
# Run the example programs to see if shared lib's get properly loaded
run: examples
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) run_print_line
	@$(MAKE) -C ex2 $(MAKECMDGOALS)
	@$(MAKE) run_print_line
	@$(MAKE) -C ex3 $(MAKECMDGOALS)
	@$(MAKE) run_print_line
	@$(MAKE) -C ex4 $(MAKECMDGOALS)
	@$(MAKE) run_print_line
	@$(MAKE) -C ex5 $(MAKECMDGOALS)
	@$(MAKE) run_print_line
	@$(MAKE) -C ex7 $(MAKECMDGOALS)
