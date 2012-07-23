# Set to appropriate C++ compiler
CPP=g++
CPPFLAGS=-I../../tecsrc ../../libtecio.a
EXECUTABLE=gridsolution
FILES=$(EXECUTABLE).cpp

build:
	$(CPP) $(FILES) $(CPPFLAGS) -o $(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE)
