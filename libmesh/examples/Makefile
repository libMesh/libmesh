all:: examples

#
# Build the example programs
examples:
	@echo "----------------------------------------------"
	@echo "--------- Building Example Programs ----------"
	@echo "----------------------------------------------"
	@$(MAKE) -C ex1  $(MAKECMDGOALS)
	@$(MAKE) -C ex2  $(MAKECMDGOALS)
	@$(MAKE) -C ex3  $(MAKECMDGOALS)
	@$(MAKE) -C ex4  $(MAKECMDGOALS)
	@$(MAKE) -C ex5  $(MAKECMDGOALS)
	@$(MAKE) -C ex6  $(MAKECMDGOALS)
	@$(MAKE) -C ex7  $(MAKECMDGOALS)
	@$(MAKE) -C ex8  $(MAKECMDGOALS)
	@$(MAKE) -C ex9  $(MAKECMDGOALS)
	@$(MAKE) -C ex10 $(MAKECMDGOALS)
	@$(MAKE) -C ex11 $(MAKECMDGOALS)
	@$(MAKE) -C ex12 $(MAKECMDGOALS)
	@$(MAKE) -C ex13 $(MAKECMDGOALS)
	@$(MAKE) -C ex14 $(MAKECMDGOALS)
	@echo "----------------------------------------------"
	@echo "------- Done Building Example Programs -------"
	@echo "----------------------------------------------"


.PHONY: clean clobber distclean

#
# Clean up the example programs
clean:
	@$(MAKE) -C ex1  $(MAKECMDGOALS)
	@$(MAKE) -C ex2  $(MAKECMDGOALS)
	@$(MAKE) -C ex3  $(MAKECMDGOALS)
	@$(MAKE) -C ex4  $(MAKECMDGOALS)
	@$(MAKE) -C ex5  $(MAKECMDGOALS)
	@$(MAKE) -C ex6  $(MAKECMDGOALS)
	@$(MAKE) -C ex7  $(MAKECMDGOALS)
	@$(MAKE) -C ex8  $(MAKECMDGOALS)
	@$(MAKE) -C ex9  $(MAKECMDGOALS)
	@$(MAKE) -C ex10 $(MAKECMDGOALS)
	@$(MAKE) -C ex11 $(MAKECMDGOALS)
	@$(MAKE) -C ex12 $(MAKECMDGOALS)
	@$(MAKE) -C ex13 $(MAKECMDGOALS)
	@$(MAKE) -C ex14 $(MAKECMDGOALS)

#
# Really clean up the example programs
clobber:
	@$(MAKE) -C ex1  $(MAKECMDGOALS)
	@$(MAKE) -C ex2  $(MAKECMDGOALS)
	@$(MAKE) -C ex3  $(MAKECMDGOALS)
	@$(MAKE) -C ex4  $(MAKECMDGOALS)
	@$(MAKE) -C ex5  $(MAKECMDGOALS)
	@$(MAKE) -C ex6  $(MAKECMDGOALS)
	@$(MAKE) -C ex7  $(MAKECMDGOALS)
	@$(MAKE) -C ex8  $(MAKECMDGOALS)
	@$(MAKE) -C ex9  $(MAKECMDGOALS)
	@$(MAKE) -C ex10 $(MAKECMDGOALS)
	@$(MAKE) -C ex11 $(MAKECMDGOALS)
	@$(MAKE) -C ex12 $(MAKECMDGOALS)
	@$(MAKE) -C ex13 $(MAKECMDGOALS)
	@$(MAKE) -C ex14 $(MAKECMDGOALS)
#
# Make the example program directories look
# like a clean distribution
distclean:
	@$(MAKE) -C ex1  $(MAKECMDGOALS)
	@$(MAKE) -C ex2  $(MAKECMDGOALS)
	@$(MAKE) -C ex3  $(MAKECMDGOALS)
	@$(MAKE) -C ex4  $(MAKECMDGOALS)
	@$(MAKE) -C ex5  $(MAKECMDGOALS)
	@$(MAKE) -C ex6  $(MAKECMDGOALS)
	@$(MAKE) -C ex7  $(MAKECMDGOALS)
	@$(MAKE) -C ex8  $(MAKECMDGOALS)
	@$(MAKE) -C ex9  $(MAKECMDGOALS)
	@$(MAKE) -C ex10 $(MAKECMDGOALS)
	@$(MAKE) -C ex11 $(MAKECMDGOALS)
	@$(MAKE) -C ex12 $(MAKECMDGOALS)
	@$(MAKE) -C ex13 $(MAKECMDGOALS)
	@$(MAKE) -C ex14 $(MAKECMDGOALS)

#
# Link the example programs
link: $(wildcard examples/ex*/ex*.C)
	@$(MAKE) -C ex1  ex1  $(MAKECMDGOALS)
	@$(MAKE) -C ex2  ex2  $(MAKECMDGOALS)
	@$(MAKE) -C ex3  ex3  $(MAKECMDGOALS)
	@$(MAKE) -C ex4  ex4  $(MAKECMDGOALS)
	@$(MAKE) -C ex5  ex5  $(MAKECMDGOALS)
	@$(MAKE) -C ex6  ex6  $(MAKECMDGOALS)
	@$(MAKE) -C ex7  ex7  $(MAKECMDGOALS)
	@$(MAKE) -C ex8  ex8  $(MAKECMDGOALS)
	@$(MAKE) -C ex9  ex9  $(MAKECMDGOALS)
	@$(MAKE) -C ex10 ex10 $(MAKECMDGOALS)
	@$(MAKE) -C ex11 ex11 $(MAKECMDGOALS)
	@$(MAKE) -C ex12 ex12 $(MAKECMDGOALS)
	@$(MAKE) -C ex13 ex13 $(MAKECMDGOALS)
	@$(MAKE) -C ex14 ex14 $(MAKECMDGOALS)

#
# Run the example programs to see if shared lib's get properly loaded
run: $(wildcard examples/ex*/ex*.C)
	@$(MAKE) -C ex1  $(MAKECMDGOALS)
	@$(MAKE) -C ex2  $(MAKECMDGOALS)
	@$(MAKE) -C ex3  $(MAKECMDGOALS)
	@$(MAKE) -C ex4  $(MAKECMDGOALS)
	@$(MAKE) -C ex5  $(MAKECMDGOALS)
	@$(MAKE) -C ex6  $(MAKECMDGOALS)
	@$(MAKE) -C ex7  $(MAKECMDGOALS)
	@$(MAKE) -C ex8  $(MAKECMDGOALS)
	@$(MAKE) -C ex9  $(MAKECMDGOALS)
	@$(MAKE) -C ex10 $(MAKECMDGOALS)
	@$(MAKE) -C ex11 $(MAKECMDGOALS)
	@$(MAKE) -C ex12 $(MAKECMDGOALS)
	@$(MAKE) -C ex13 $(MAKECMDGOALS)
	@$(MAKE) -C ex14 $(MAKECMDGOALS)
