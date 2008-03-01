# $Id$

test_source := $(shell find ./ -name "*.cc")
tests       := $(patsubst %.cc, %, ${test_source}) 


all: run_tests


build_tests: $(tests)


run_tests: build_tests
	for test in $(tests); do $$test; done


%: %.cc driver.C
	@echo "Building $@ ..."
	@`${LIBMESH_ROOT}/contrib/bin/libmesh-config --cxx` driver.C $< `${LIBMESH_ROOT}/contrib/bin/libmesh-config --cxxflags --include` -o $@ `${LIBMESH_ROOT}/contrib/bin/libmesh-config --ldflags`




echo:
	@echo "test_source = ${test_source}"
	@echo "tests       = ${tests}"



clean:
	rm -f *~ */*~
	rm -f $(tests)
