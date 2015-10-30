# Set to appropriate C++ compiler
CPP=g++
CPPFLAGS=-I../../tecsrc ../../libtecio.a
EXECUTABLE=comtest
FILES=$(EXECUTABLE).c

build:
	$(CPP) $(FILES) $(CPPFLAGS) -o $(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE)
