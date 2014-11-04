// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

// $CXX -o autodiff autodiff.cc `$LIBMESH_DIR/bin/libmesh-config --cppflags --cxxflags --include --libs` -I $LIBMESH_DIR/include/libmesh

#include "libmesh/fparser_ad.hh"

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
#ifndef FUNCTIONPARSER_SUPPORT_DEBUGGING
  std::cout << "Please configure libMesh with --enable-fparser-debugging and rebuild to enable bytecode output.\n";
#endif

  std::string function;
  FunctionParserAD fparser;

  std::string func = "2 + 4*x + 8*x^2 + 16*x^3 + 32*x^4";

  // Parse the input expression into bytecode
  fparser.Parse(func, "x,a");
  std::cout << "Input Expression:\n" << func << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  fparser.PrintByteCode(std::cout);
#endif
  std::cout << std::endl;
  write("input_orig.dat", fparser);

  // output optimized version
  fparser.Optimize();
  std::cout << "Optimized Input Expression:\n" << func << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  fparser.PrintByteCode(std::cout);
#endif
  std::cout << std::endl;
  write("input_opt.dat", fparser);

  // Get a copy of the original function
  FunctionParserAD fparser2(fparser);

  // Generate derivative with respect to x
  std::cout << "AutoDiff returned " << fparser.AutoDiff("x") << std::endl;
  std::cout << "Unsimplified derivative:\n";
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  fparser.PrintByteCode(std::cout);
#endif
  std::cout << std::endl;
  write("first_orig.dat", fparser);

  // Run the bytecode optimizer on the derivative
  fparser.Optimize();
  std::cout << "Simplified derivative:\n";
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  fparser.PrintByteCode(std::cout);
#endif
  std::cout << std::endl;
  write("first_opt.dat", fparser);

  // Check if the copied instance is still the original function
  /*std::cout << "Copy of the original:\n";
  fparser2.PrintByteCode(std::cout);
  std::cout << std::endl;*/

  // Generate derivative with respect to x
  std::cout << "AutoDiff 2nd time returned " << fparser.AutoDiff("x") << std::endl;
  std::cout << "Simplified second derivative:\n";
  fparser.Optimize();
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  fparser.PrintByteCode(std::cout);
#endif
  std::cout << std::endl;

  fparser2.JITCompile();
  write("input_jit.dat", fparser2);

  return 0;
}
