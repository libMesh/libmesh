// Simple example file for the the fparser optimization bug
// (c) 2014 by Daniel Schwen
// ========================================================

// $CXX -o optimizer_bug optimizer_bug.cc -I ../../../installed/include/
// ../../../build/contrib/fparser/.libs/libdbg.a

#include "../fparser.hh"

#include <iostream>
#include <fstream>
#include <string>

int main()
{
  std::ofstream myfile;
  double p[3] = {0.0, 2.0, 0.1};
  const double step = 0.001;
  FunctionParser fparser;

  // // fails:
  // std::string func = "a := 1.81e+16* exp(-1000*abs(c)); 1- 1/(1 + a)";
  // std::string func = "1 - (2e+16 * (3 ^ (-abs(c)) )) ^ -1";
  std::string func = "1 - (5e-17 * (0.3 ^ (-abs(c))))";
  // std::string func = "a := 1.81e+16* exp(-1000*abs(c)); 1- (1/(1 + a))";
  // std::string func = "a := exp(-1000*(abs(c)-0.03742995)); 1- (1/(1 + a))";
  // std::string func = "exp(-100*(abs(c)-8))";

  // // works:
  // std::string func = "a := 1.8e+16* exp(-1000*abs(c)); 1- (1/(1 + a))";
  // std::string func = "a := exp(-1000*(abs(c)-0.03742994)); 1- (1/(1 + a))";

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
  for (p[0] = -0.1; p[0] <= 0.1; p[0] += step)
    myfile << p[0] << ' ' << fparser.Eval(p) << std::endl;
  myfile.close();

  // compute difference
  double diff = 0.0;
  for (p[0] = -0.1; p[0] <= 0.1; p[0] += step)
    diff += step * (fparser.Eval(p) - fparser2.Eval(p));

  std::cout << "Integral of f_optimized()-f_original() over [-0.1,0.1] = "
            << diff << "  (should be ~0.0)" << std::endl;
  return 0;
}
