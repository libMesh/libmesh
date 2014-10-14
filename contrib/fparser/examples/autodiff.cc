// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

// $CXX -o autodiff autodiff.cc -I $LIBMESH_DIR/include/libmesh -DFUNCTIONPARSER_SUPPORT_DEBUGGING ../../../build_dev/contrib/fparser/.libs/*.o

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
  //std::string func = "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))+x^3*cos(x)-log(x+2*x^2+3.1*x^3)-x^5*(2+x-cos(x-x^2))";
  //std::string func = "if(x<x*x,x^(-1/2),5)";
  //std::string func = "5+if(x<x*x,1,3)*6+x";
  //std::string func = "sin(x^2)";
  //std::string func = "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))";
  //std::string func = "if(x<x*x,x^(-1/2),5)";
  //std::string func = "5+if(x<x*x,1,3)*6+x";
  //std::string func = "1+2*plog(3,4)";
  std::string func = "1+2*plog(x,0.3)";


  // Parse the input expression into bytecode
  fparser.Parse(func, "x,a");
  std::cout << "Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;
  write("input_orig.dat", fparser);

  // output optimized version
  fparser.Optimize();
  std::cout << "Optimized Input Expression:\n" << func << std::endl;
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;
  write("input_opt.dat", fparser);

  // Get a copy of the original function
  FunctionParserAD fparser2(fparser);

  // Generate derivative with respect to x
  std::cout << "AutoDiff returned " << fparser.AutoDiff("x") << std::endl;
  std::cout << "Unsimplified derivative:\n";
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;
  write("first_orig.dat", fparser);

  // Run the bytecode optimizer on the derivative
  fparser.Optimize();
  std::cout << "Simplified derivative:\n";
  fparser.PrintByteCode(std::cout);
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
  fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  fparser2.JITCompile();
  write("input_jit.dat", fparser2);

  return 0;
}
