#include "libmesh/fparser_ad.hh"

// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

class FParserAutodiffTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE ( FParserAutodiffTest );

  CPPUNIT_TEST ( runTests );
  CPPUNIT_TEST ( registerDerivativeTest );

  CPPUNIT_TEST_SUITE_END ();

private:

  class ADTest {
  public:
    ADTest(const std::string & _func,
           double _min, double _max, double _dx = 1e-6,
           double _reltol = 1e-5, int _steps = 20, double _abstol = 1e-10) :
      func(_func),
      min(_min),
      max(_max),
      dx(_dx),
      reltol(_reltol),
      abstol(_abstol),
      steps(_steps)
    {
      CPPUNIT_ASSERT_MESSAGE ("Failed to parse test function", F.Parse(func, "x") == -1);
      dF.Parse(func, "x");
      dFopt.Parse(func, "x");
      dFaopt.Parse(func, "x");

      CPPUNIT_ASSERT_MESSAGE ("Failed to take derivative of function", dF.AutoDiff("x") == -1);

      dFopt.Optimize();
      CPPUNIT_ASSERT_MESSAGE ("Failed to take derivative of optimized function", dFopt.AutoDiff("x") == -1);

      dFaopt.SetADFlags(FunctionParserAD::ADAutoOptimize, true);
      CPPUNIT_ASSERT_MESSAGE ("Failed to take derivative of auto-optimized function", dFaopt.AutoDiff("x") == -1);
    }

    bool run()
    {
      double x1, x2, vdF, vF1, vF2, fd;
      for (double x = min; x <= max; x += (max-min) / double(steps))
        {
          x1 = x - dx/2.0;
          x2 = x + dx/2.0;

          vF1 = F.Eval(&x1);
          vF2 = F.Eval(&x2);
          fd = (vF2-vF1) / dx;

          // CPPUNIT_ASSERT(std::abs(fd - vdF) > tol)
          // CPPUNIT_ASSERT(std::abs(fd - vdFopt) > tol)
          vdF = dF.Eval(&x);
          if (std::abs(vdF) > abstol && std::abs((fd - vdF)/vdF) > reltol && std::abs(fd - vdF)> abstol)
            {
              std::cout << "Error in " << func << ": " << fd << "!=" << vdF << " at x=" << x << '\n';
              return false;
            }

          vdF = dFopt.Eval(&x);
          if (std::abs(vdF) > abstol && std::abs((fd - vdF)/vdF) > reltol && std::abs(fd - vdF)> abstol)
            {
              std::cout << "Error in opt " << func << ": " << fd << "!=" << vdF << " at x=" << x << '\n';
              return false;
            }

          vdF = dFaopt.Eval(&x);
          if (std::abs(vdF) > abstol && std::abs((fd - vdF)/vdF) > reltol && std::abs(fd - vdF)> abstol)
            {
              std::cout << "Error in auto opt " << func << ": " << fd << "!=" << vdF << " at x=" << x << '\n';
              return false;
            }
        }

      return true;
    }

  private:
    std::string func;
    double min, max, dx, reltol, abstol;
    int steps;

    FunctionParserAD F, dF, dFopt, dFaopt;
  };

  std::vector<ADTest> tests;

public:
  virtual void setUp()
  {
    tests.push_back(ADTest("log(x*x) + log2(2*x) + log10(4*x)", 0.1, 3.0));
    tests.push_back(ADTest("sin(-x) + cos(2*x) + tan(4*x)", -5.0, 5.0, 1e-7, 1e-5, 100));
    tests.push_back(ADTest("sinh(-x) + cosh(x/2) + tanh(x/3)", -4.0, 4.0, 0.0001, 1e-5, 100));
    tests.push_back(ADTest("plog(-x,0.01)", 0.001, 0.05, 0.00001, 1e-5, 100));
    tests.push_back(ADTest("2 + 4*x + 8*x^2 + 16*x^3 + 32*x^4", -5.0, 5.0, 1e-5,1e-4));
    tests.push_back(ADTest("1/x^2", 0.01, 2.0, 1e-8));
    tests.push_back(ADTest("sqrt(x*2)", 0.001, 2.0, 1e-6));
    tests.push_back(ADTest("abs(x*2)", -1.99, 2.0));
    tests.push_back(ADTest("asin(-x)", -0.99, 0.99));
    tests.push_back(ADTest("acos(-x)", -0.99, 0.99));
    tests.push_back(ADTest("atan(-x)", -99, 99));
    tests.push_back(ADTest("x*sin(-x)*log(x)*tanh(x)", 0.001, 5, 1e-8));
    tests.push_back(ADTest("exp(-x) + 2*exp2(x)", -1.0, 2.0));
    tests.push_back(ADTest("hypot(2*x,1) - hypot(1,4*x)", -10, 10.0));
    tests.push_back(ADTest("if(x<0, (-x)^3, x^3)", -1.0, 1.0));
    tests.push_back(ADTest("max(x^2-0.5,0)", -1.5, 1.5));
    tests.push_back(ADTest("min(x^2-0.5,0)", -1.5, 1.5));
    tests.push_back(ADTest("atan2(x,1) + atan2(2,x)", -0.99, 0.99));
    tests.push_back(ADTest("0.767^sin(x)", -1.5, 1.5));
    tests.push_back(ADTest("A := sin(x) + tanh(x); A + sqrt(A) - x", -1.5, 1.5));
    tests.push_back(ADTest("3*sin(2*x)*sin(2*x)", -5.0, 5.0, 1e-7, 1e-5, 100));
  }

  void runTests()
  {
    const unsigned int ntests = tests.size();

    unsigned int passed = 0;
    for (unsigned i = 0; i < ntests; ++i)
      passed += tests[i].run() ? 1 : 0;

    CPPUNIT_ASSERT_EQUAL (passed, ntests);
  }

  void registerDerivativeTest()
  {
    FunctionParserAD R;
    std::string func = "x*a";

    // Parse the input expression into bytecode
    R.Parse(func, "x,a");

    // add a new variable y and map it to the da/dx derivative
    R.AddVariable("y");
    R.RegisterDerivative("a", "x", "y");

    // parameter vector
    double p[3];
    double & x = p[0];
    double & a = p[1];
    double & y = p[2];

    FunctionParserAD dR(R);
    CPPUNIT_ASSERT_EQUAL (dR.AutoDiff("x"), -1);
    dR.Optimize();
    // dR = a+x*y

    FunctionParserAD d2R(dR);
    CPPUNIT_ASSERT_EQUAL (d2R.AutoDiff("x"), -1);
    d2R.Optimize();
    // d2R = 2*y

    // we probe the parsers and check if they agree with the reference solution
    for (x = -1.0; x < 1.0; x+=0.3726)
      for (a = -1.0; a < 1.0; a+=0.2642)
        for (y = -1.0; y < 1.0; y+=0.3156)
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(R.Eval(p), x*a, 1.e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dR.Eval(p), a+x*y, 1.e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(d2R.Eval(p), 2*y, 1.e-12);
          }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION ( FParserAutodiffTest );
