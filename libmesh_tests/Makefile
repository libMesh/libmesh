LIBMESH_DIR ?= ..

CC=$(shell ${LIBMESH_DIR}/contrib/bin/libmesh-config --cxx)
CFLAGS=-c -Wall $(shell ${LIBMESH_DIR}/contrib/bin/libmesh-config --cxxflags --include$) -I${CPP_UNIT}/include
LDFLAGS=-L${CPP_UNIT}/lib -lcppunit $(shell ${LIBMESH_DIR}/contrib/bin/libmesh-config --ldflags)
SOURCES=$(shell find ./ -name "*.C")
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=run_tests
DEPEND= makedepend

all: $(EXECUTABLE)
	./$(EXECUTABLE)
	
$(EXECUTABLE): depend $(SOURCES) $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.C.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *~* */*~* *.o */*.o run_tests .depends

depend: $(SOURCES)
	@touch .depends
	@$(DEPEND) -f.depends -- $(CFLAGS) -- $(SOURCES) &>/dev/null

-include .depends
# DO NOT DELETE
