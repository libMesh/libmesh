// Simple example file for the auto-differentiation module
// (c) 2014 by Daniel Schwen
// =======================================================

#include "../fparser_ad.hh"

#include <iostream>
#include <string>
#include <sys/time.h>


// clang++ jittest.cc -o jittest -L../../../build/contrib/fparser/.libs -lopt -I ../../../installed/include/libmesh/ -DFUNCTIONPARSER_SUPPORT_DEBUGGING

int main()
{
  std::string function;
  FunctionParserAD fparser, fparser2;

  std::string func = "sin(3*x)+x*5*sin(3*x)+x^2*(3+sin(3*x))+x^3*cos(x)-log(x+2*x^2+3.1*x^3)-x^5*(2+x-cos(x-x^2))";

  // Parse the input expression into bytecode
  fparser.Parse(func, "x");
  std::cout << "Input Expression:\n" << func << std::endl;
  //fparser.PrintByteCode(std::cout);
  std::cout << std::endl;

  fparser2.Parse(func, "x");
  fparser2.Optimize();

  double p[1];
  timeval time0, time1, time2, time3, time4;
  const int nmax = 1000; //10000;

  gettimeofday(&time0, NULL);
  for (int n=0; n<nmax; ++n)
    for (p[0]=0.1; p[0]<100.0; p[0]+=0.1) 
      fparser.Eval(p);

  gettimeofday(&time1, NULL);

  bool jit = fparser.JITCompile();

  if (jit)
  {
    gettimeofday(&time2, NULL);

    for (int n=0; n<nmax; ++n)
      for (p[0]=0.1; p[0]<100.0; p[0]+=0.1) 
        fparser.Eval(p);
  }

  gettimeofday(&time3, NULL);

  for (int n=0; n<nmax; ++n)
    for (p[0]=0.1; p[0]<100.0; p[0]+=0.1) 
      fparser2.Eval(p);

  gettimeofday(&time4, NULL);


  long t0 = (time0.tv_sec * 1000) + (time0.tv_usec / 1000);
  long t1 = (time1.tv_sec * 1000) + (time1.tv_usec / 1000);
  long t2 = (time2.tv_sec * 1000) + (time2.tv_usec / 1000);
  long t3 = (time3.tv_sec * 1000) + (time3.tv_usec / 1000);
  long t4 = (time4.tv_sec * 1000) + (time4.tv_usec / 1000);

  std::cout << "Unoptimized: " << (t1-t0) << " ms\n";
  if (jit) {
    std::cout << "Compile time:" << (t2-t1) << " ms\n";
    std::cout << "JIT Compiled:" << (t3-t2) << " ms\n";
  }
  std::cout << "FP Optimized:" << (t4-t3) << " ms\n";
  return 0;
}
