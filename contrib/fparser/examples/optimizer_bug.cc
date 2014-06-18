// Simple example file for the the fparser optimization bug
// (c) 2014 by Daniel Schwen
// ========================================================

// $CXX -o optimizer_bug optimizer_bug.cc -I ../../../installed/include/libmesh/ ../../../build/contrib/fparser/.libs/libfparser_la-fp* 

#include "../fparser.hh"

#include <iostream>
#include <fstream>
#include <string>

int main()
{
  std::ofstream myfile;
  double p[3] = {0.0, 2.0, 0.1};
  const double step = 0.01;
  FunctionParser fparser;

  // The White Whale
  std::string func = "c*2*(1-c)^2 - c^2*(1-c)*2";

  // Parse the input expression into bytecode
  fparser.Parse(func, "c");
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  std::cout << "Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;
#endif

  // eval into file
  myfile.open ("pre.dat");
  for (p[0]=0.0; p[0]<=1.0; p[0]+=0.01) myfile << p[0] << ' ' << fparser.Eval(p) << std::endl;
  myfile.close();

  FunctionParser fparser2(fparser);

  // output optimized version
  fparser.Optimize();
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
  std::cout << "Optimized Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;
#endif

  // eval into file
  myfile.open ("post.dat");
  for (p[0]=0.0; p[0]<=1.0; p[0]+=0.01) myfile << p[0] << ' ' << fparser.Eval(p) << std::endl;
  myfile.close();

  // compute difference
  double diff = 0.0;
  for (p[0]=0.0; p[0]<=1.0; p[0]+=step) diff+=step*(fparser.Eval(p)-fparser2.Eval(p));

  std::cout << "Integral of f_optimized()-f_original() on [0,1] = " << diff << "  (should be ~0.0)" << std::endl;
  return 0;
}
