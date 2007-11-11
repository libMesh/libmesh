CC=$(shell ${LIBMESH_ROOT}/contrib/bin/libmesh-config --cxx)
CFLAGS=-g -c -Wall $(shell ${LIBMESH_ROOT}/contrib/bin/libmesh-config --cxxflags --include$)
LDFLAGS=-lcppunit $(shell ${LIBMESH_ROOT}/contrib/bin/libmesh-config --ldflags)
SOURCES=$(shell find ./ -name "*.C")
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=run_tests
DEPEND= makedepend

all: $(EXECUTABLE)
	
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
