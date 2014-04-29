// Simple example file for teh auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

#include "../fparser_ad.hh"

#include <iostream>
#include <string>

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
  std::string func = "abs(3*x^3)";

  //std::string func = "sin(x^2)";
  //std::string func = "sqrt(3*x)";

  std::cout << "Input Expression:\n" << func << std::endl;

  // Parse the input expression into bytecode
  fparser.Parse(func, "x,a");
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  // Get a copy of the original function
  FunctionParserAD fparser2(fparser);

  // Generate derivative with respect to x
  fparser.AutoDiff("x");
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

  return 0;
}
