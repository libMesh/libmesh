CC=g++
CFLAGS=-c -Wall $(shell ${LIBMESH_ROOT}/contrib/bin/libmesh-config --cxxflags --include$)
LDFLAGS=-lcppunit $(shell ${LIBMESH_ROOT}/contrib/bin/libmesh-config --ldflags)
SOURCES=$(shell find ./ -name "*.C")
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=run_tests
DEPEND= makedepend

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) depend
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.C.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *~* */*~* *.o */*.o run_tests

depend: $(SOURCES)
	@touch .depends
	@$(DEPEND) -f.depends -- $(CFLAGS) -- $(SOURCES) &>/dev/null

-include .depends
