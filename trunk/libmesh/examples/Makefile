all:: examples

#
# Build the example programs
examples:
	@echo "----------------------------------------------"
	@echo "--------- Building Example Programs ----------"
	@echo "----------------------------------------------"
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) -C ex2 $(MAKECMDGOALS)
	@$(MAKE) -C ex3 $(MAKECMDGOALS)
	@$(MAKE) -C ex4 $(MAKECMDGOALS)
	@$(MAKE) -C ex5 $(MAKECMDGOALS)
	@$(MAKE) -C ex7 $(MAKECMDGOALS)
	@echo "----------------------------------------------"
	@echo "------- Done Building Example Programs -------"
	@echo "----------------------------------------------"


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



Example%:
	@echo ""
	@echo "----------------------------------------------"
	@echo "---- Finished running $(MAKECMDGOALS). Proceed  -----"
	@echo "----------------------------------------------"
	@echo ""


#
# Run the example programs to see if shared lib's get properly loaded
run: $(wildcard examples/ex?/ex?.C)
	@$(MAKE) -C ex1 $(MAKECMDGOALS)
	@$(MAKE) Example1
	@$(MAKE) -C ex2 $(MAKECMDGOALS)
	@$(MAKE) Example2
	@$(MAKE) -C ex3 $(MAKECMDGOALS)
	@$(MAKE) Example3
	@$(MAKE) -C ex4 $(MAKECMDGOALS)
	@$(MAKE) Example4
	@$(MAKE) -C ex5 $(MAKECMDGOALS)
	@$(MAKE) Example5
	@$(MAKE) -C ex7 $(MAKECMDGOALS)
	@echo ""
	@echo "----------------------------------------------"
	@echo "--------- Finished running Example7. ---------
	@echo "----------------------------------------------"
	@echo ""

