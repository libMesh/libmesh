// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

// $CXX -o autodiff2 autodiff2.cc `$LIBMESH_DIR/bin/libmesh-config --cppflags --cxxflags --include --libs` -DFUNCTIONPARSER_SUPPORT_DEBUGGING -I $LIBMESH_DIR/include/libmesh

#include "../fparser_ad.hh"

#include <fstream>
#include <iostream>
#include <string>

void write(const char *fname, FunctionParserAD & F)
{
  std::ofstream file;
  file.open(fname);
  for (double x = -1.0; x <= 1.0; x+=0.01)
    file << x << ' ' << F.Eval(&x) << '\n';
  file.close();
}

int main()
{
  std::string function;
  FunctionParserAD fparser;

  //std::string func = "1+2*if(x<0,x,x^2)";
  std::string func = "1+2*x+x+5*x";


  // Parse the input expression into bytecode
  fparser.Parse(func, "x,a");
  std::cout << "Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  fparser.AutoDiff2("x");
  std::cout << std::endl;

  return 0;
}
