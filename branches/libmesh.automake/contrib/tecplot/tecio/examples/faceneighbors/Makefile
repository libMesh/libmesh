# Set to appropriate C++ compiler
CPP=g++
CPPFLAGS=-I../../tecsrc ../../tecio.a
EXECUTABLE=faceneighbors
FILES=$(EXECUTABLE).cpp

build:
	$(CPP) $(FILES) $(CPPFLAGS) -o $(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE)