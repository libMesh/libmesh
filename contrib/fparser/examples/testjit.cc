// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

#include "../fparser_ad.hh"

#include <iostream>
#include <string>
#include <cstdlib>

// clang++ testjit.cc -o testjit -L../../../build/contrib/fparser/.libs -lopt -I ../../../installed/include/libmesh/ -DFUNCTIONPARSER_SUPPORT_DEBUGGING

int main()
{
  const int nfunc = 11;
  const char *func[nfunc] = {
    "if(x>0,x+7,if(x<1.001,x,2*x^2)-9)+100",
    "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))+x^3*cos(-x)-log(2*x^2+3.1*x^4)-x^5*(2+x-cos(x-x^2))",
    "(x<= 1.8) & (x>=0.6)",
    "(x> 1.8) | (x<0.6)",
    "hypot(x+0.77,0.5*x^4)",
    "log10(1/abs(x)+0.1)",
    "x^(1/3)",
    "x % 0.7 + exp(x)",
    "min(x,-0.7) + max(x,1.2)",
    "floor(x)*10 + ceil(x) + 100*int(x) + 1000*trunc(x)",
    "cos(x+3)+sin(x+3)+4*cot(x)+8*sec(x)"
  };

  for (unsigned int i = 0; i < nfunc; ++i)
  {
    FunctionParserAD fparser, fparser2, fparser3;

    // Parse the input expression into bytecode
    if (fparser.Parse(func[i], "x") != -1)
    {
      std::cout << "Parsing failed for " << func[i] << '\n';
      exit(1);
    }
    fparser2.Parse(func[i], "x");
    fparser2.Optimize();
    fparser3.Parse(func[i], "x");

    if (!fparser.JITCompile())
    {
      std::cout << "JIT compile failed for " << func[i] << '\n';
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
      fparser.PrintByteCode(std::cout);
#endif
      exit(1);
    }
    if (!fparser2.JITCompile())
    {
      std::cout << "JIT compile failed for optimized " << func[i] << '\n';
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
      fparser2.PrintByteCode(std::cout);
#endif
      exit(1);
    }

    const double eps = 1e-12;
    double p[1], a, b, c;

    for (p[0]=-2.7; p[0]<2.7; p[0]+=0.1)
    {
      a = fparser.Eval(p);
      b = fparser2.Eval(p);
      c = fparser3.Eval(p);

      if (std::abs(a-c) > eps)
      {
        std::cout << "Eval discrepancy for " << func[i] << '\n';
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
        fparser.PrintByteCode(std::cout);
#endif
        exit(1);
      }
      if (std::abs(b-c) > eps)
      {
        std::cout << "Eval discrepancy for optimized " << func[i] << '\n';
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
        fparser2.PrintByteCode(std::cout);
#endif
        exit(1);
      }
    }
  }

  std::cout << "Success.\n";
  return 0;
}
