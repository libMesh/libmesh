LIBMESH_DIR ?= ..

CXX=$(shell ${LIBMESH_DIR}/contrib/bin/libmesh-config --cxx)
CFLAGS=-c -Wall $(shell ${LIBMESH_DIR}/contrib/bin/libmesh-config --cxxflags --include$) -I${CPP_UNIT}/include
LDFLAGS=-L${CPP_UNIT}/lib -lcppunit $(shell ${LIBMESH_DIR}/contrib/bin/libmesh-config --ldflags)
SOURCES=$(shell find ./ -name "*.C")
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=run_tests

all: $(EXECUTABLE)
	$(LIBMESH_RUN) ./$(EXECUTABLE) $(LIBMESH_OPTIONS)
	
$(EXECUTABLE): depend $(SOURCES) $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.C.o:
	$(CXX) $(CFLAGS) $< -o $@

clean:
	rm -rf *~* */*~* *.o */*.o run_tests .depends

depend: $(SOURCES)
	touch .depends
	@$(CXX) $(CFLAGS) -M $(SOURCES) > .depends

-include .depends
# DO NOT DELETE
