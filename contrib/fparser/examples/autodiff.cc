// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

#include "../fparser_ad.hh"

#include <iostream>
#include <string>

// clang++ autodiff.cc -o autodiff -L../../../build/contrib/fparser/.libs -lopt -I ../../../installed/include/libmesh/ -DFUNCTIONPARSER_SUPPORT_DEBUGGING

int main()
{
  std::string function;
  FunctionParserAD fparser;

  //std::string func = "x^2+3.345*x+4*a+3";
  //std::string func = "(x^2+3.345*x)*(4*a+3)";
  //std::string func = "3*x^2+4*x+5";
  //std::string func = "4*log(4*x^2+6)-7*x";

  // works:
  //std::string func = "log(2+8*x^2)";
  //std::string func = "(4*x+8*x^2)+(3*x+4*x^2+7)";
  //std::string func = "x*log(x)+(1-x)*log(1-x)";

  //std::string func = "3*a^2 + 7*x^2";
  //std::string func = "3*x + 4*x^2";
  //std::string func = "hypot(3*x,7*sin(x))";
  //std::string func = "if(x-3,x^2,if(x-4,1,7*x))";
  //std::string func = "x^2+x^2*3";

  //std::string func = "sin(x-2)";
  std::string func = "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))+x^3*cos(x)-log(x+2*x^2+3.1*x^3)-x^5*(2+x-cos(x-x^2))";
  //std::string func = "if(x<x*x,x^(-1/2),5)";
  //std::string func = "5+if(x<x*x,1,3)*6+x";


  // Parse the input expression into bytecode
  fparser.Parse(func, "x,a");
  std::cout << "Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  // output optimized version
  fparser.Optimize();
  std::cout << "Optimized Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  // Get a copy of the original function
  FunctionParserAD fparser2(fparser);

  // Generate derivative with respect to x
  std::cout << "AutoDiff returned " << fparser.AutoDiff("x") << std::endl;
  std::cout << "Unsimplified derivative:\n";
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  // Run the bytecode optimizer on the derivative
  fparser.Optimize();
  std::cout << "Simplified derivative:\n";
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  // Check if the copied instance is still the original function
  std::cout << "Copy of the original:\n";
  fparser2.PrintByteCode(std::cout);
  std::cout << std::endl;
  
  double p[1];
  
  for (p[0]=0.0; p[0]<3.141; p[0]+=0.1) 
    std::cout << fparser2.Eval(p) << '\n';

  fparser2.JITCompile();
  std::cout << "\nJIT\n\n";

  for (p[0]=0.0; p[0]<3.141; p[0]+=0.1) 
    std::cout << fparser2.Eval(p) << '\n';

  return 0;
}
