# Set to appropriate C++ compiler
CPP=g++
CPPFLAGS=-I../../tecsrc ../../tecio.a
EXECUTABLE=simtest
FILES=$(EXECUTABLE).c

build:
	$(CPP) $(FILES) $(CPPFLAGS) -o $(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE)
