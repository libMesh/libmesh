// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

#include "../fparser_ad.hh"

#include <iostream>
#include <string>
#include <sys/time.h>
#include <ctime>
#include <cstdlib>


// clang++ perfjit.cc -o perfjit -L../../../build/contrib/fparser/.libs -lopt -I ../../../installed/include/libmesh/ -DFUNCTIONPARSER_SUPPORT_DEBUGGING

double evalTime(FunctionParserAD & fparser)
{
  double p[1];
  std::clock_t start= std::clock();

  for (int i=0; i<100000; ++i)
    for (p[0]=-2.7; p[0]<2.7; p[0]+=0.1)
      fparser.Eval(p);

  return std::clock() - start;
}

void testFParser(const char *func[], int nfunc, bool ad)
{
  double t=0, to=0, tj=0, toj=0, tcms = 0.0;
  struct timeval tc0, tc1;

  for (unsigned int i = 0; i < nfunc; ++i)
  {
    FunctionParserAD fparser, fparser_o, fparser_j, fparser_oj;

    // Parse
    if (fparser.Parse(func[i], "x") != -1) exit(1);
    fparser_o.Parse(func[i], "x");
    fparser_j.Parse(func[i], "x");
    fparser_oj.Parse(func[i], "x");

    // derivative
    if (ad) {
      if (fparser.AutoDiff("x") != -1) {
        std::cout << func[i] << '\n';
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
        fparser_j.PrintByteCode(std::cout);
#endif
        exit(1);
      }
      fparser_o.AutoDiff("x");
      fparser_j.AutoDiff("x");
      fparser_oj.AutoDiff("x");
    }

    // optimize
    fparser_o.Optimize();
    fparser_oj.Optimize();

    // JIT compile
    gettimeofday(&tc0,NULL);
    if (!fparser_j.JITCompile()) exit(1);
    gettimeofday(&tc1,NULL);
    if (!fparser_oj.JITCompile()) exit(1);

    double tc0ms = (double)tc0.tv_sec * 1000.0 + (double)tc0.tv_usec * .001,
           tc1ms = (double)tc1.tv_sec * 1000.0 + (double)tc1.tv_usec * .001;
    tcms += tc1ms - tc0ms;

    // time
    t   += evalTime(fparser);
    to  += evalTime(fparser_o);
    tj  += evalTime(fparser_j);
    toj += evalTime(fparser_oj);
  }

  std::cout << "FParser                  1 \n";
  std::cout << "FParser optimized        " << (t/to) << " x faster\n";
  std::cout << "FParser JIT              " << (t/tj) << " x faster\n";
  std::cout << "FParser optimized + JIT  " << (t/toj) << " x faster\n";
  std::cout << "   (mean compile time " << (tcms/nfunc) << " ms)\n";
}

void testRegular()
{
  const int nfunc = 7;
  const char *func[nfunc] = {
    "if(x>0,x+7,if(x<1.001,x,2*x^2)-9)+100",
    "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))+x^3*cos(x)-log(2*x^2+3.1*x^4)-x^5*(2+x-cos(x-x^2))",
    "(x<= 1.8) & (x>=0.6)",
    "(x> 1.8) | (x<0.6)",
    "hypot(x+0.77,0.5*x^4)",
    "log10(abs(x))",
    "x^(1/3)"
  };

  testFParser(func, nfunc, false);
}

void testDerivatives()
{
  const int nfunc = 4;
  const char *func[nfunc] = {
    "log(x)*x+log(1-x)*(1-x)",
    "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))+x^3*cos(x)-log(2*x^2+3.1*x^4)-x^5*(2+x-cos(x-x^2))",
    "atan(sin(x)*x^3+4*x-cos(x)+6)*x^2",
    "abs(sin(3+4*x+5*x^2+6*x^3))^(0.03*x^2+0.02*x^4+0.01*x^6)"
  };

  testFParser(func, nfunc, true);
}

int main()
{
  std::cout << "----- regular functions ----\n";
  testRegular();

  std::cout << "----- derivative functions ----\n";
  testDerivatives();

  return 0;
}
