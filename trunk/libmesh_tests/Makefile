test_source := $(shell find ./ -name "*.cc")
tests       := $(patsubst %.cc, %, ${test_source}) 

build_test := $(echo $(test))
run_test   := $(echo "Running ./$(test)" ; ./$(test))



all: $(tests) $(tests_source) driver.C
	$(foreach test, $(tests), $(build_test))


%: %.cc driver.C
	@echo "Building $@ ..."
	@`${LIBMESH_ROOT}/libmesh/contrib/bin/libmesh-config --cxx` driver.C $< `${LIBMESH_ROOT}/libmesh/contrib/bin/libmesh-config --cxxflags --include` -o $@ `${LIBMESH_ROOT}/libmesh/contrib/bin/libmesh-config --ldflags`




echo:
	@echo "test_source = ${test_source}"
	@echo "tests       = ${tests}"



clean:
	rm -f *~ */*~
	rm -f $(tests)