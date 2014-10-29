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
  //std::string func = "1+2*x+x*x+5*x^3";
  //std::string func = "A := x+1; sin(A)+cos(A)";
  std::string func = "if(x<0,3*x,sqrt(x))";

  // Parse the input expression into bytecode
  if (fparser.Parse(func, "x,a") != -1)
    return 1;
  std::cout << "Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  fparser.AutoDiff2("x");
  std::cout << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  return 0;
}
