/*==========================================================================
  testbed
  ---------
  Copyright: Juha Nieminen, Joel Yliluoma
  This program (testbed) is distributed under the terms of
  the GNU General Public License (GPL) version 3.
  See gpl.txt for the license text.
============================================================================*/

static const char* const kVersionNumber = "2.3.0.12";

#include "fpconfig.hh"
#include "fparser.hh"

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
#include "fparser_mpfr.hh"
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
#include "fparser_gmpint.hh"
#endif

#include "extrasrc/fpaux.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <cassert>

#define CONST 1.5

#define StringifyHlp(x) #x
#define Stringify(x) StringifyHlp(x)

#ifndef FP_DISABLE_DOUBLE_TYPE
typedef FunctionParser DefaultParser;
#elif defined(FP_SUPPORT_LONG_DOUBLE_TYPE)
typedef FunctionParser_ld DefaultParser;
#elif defined(FP_SUPPORT_FLOAT_TYPE)
typedef FunctionParser_f DefaultParser;
#elif defined(FP_SUPPORT_MPFR_FLOAT_TYPE)
typedef FunctionParser_mpfr DefaultParser;
#else
#error "FunctionParserBase<double> was disabled and no viable floating point alternative has been defined"
#endif

typedef DefaultParser::value_type DefaultValue_t;


namespace
{
    /* Verbosity level:
       0 = No progress output. Error reporting as in verbosity level 1.
       1 = Very brief progress and error output.
       2 = More verbose progress output, full error reporting.
       3 = Very verbose progress output, full error reporting.
    */
    int verbosityLevel = 1;


    const char* getEvalErrorName(int errorCode)
    {
        static const char* const evalErrorNames[6] =
        {
            "no error", "division by zero", "sqrt error", "log error",
            "trigonometric error", "max eval recursion level reached"
        };
        if(errorCode >= 0 && errorCode < 6)
            return evalErrorNames[errorCode];
        return "unknown";
    }

    std::vector<const char*> selectedRegressionTests;

    // Auxiliary functions
    // -------------------
    template<typename Value_t>
    inline Value_t r2d(Value_t x)
    { return x * (Value_t(180) / FUNCTIONPARSERTYPES::fp_const_pi<Value_t>()); }

    template<typename Value_t>
    inline Value_t d2r(Value_t x)
    { return x * (FUNCTIONPARSERTYPES::fp_const_pi<Value_t>() / Value_t(180)); }

    //inline double log10(double x) { return std::log(x) / std::log(10); }

    template<typename Value_t>
    Value_t userDefFuncSqr(const Value_t* p) { return p[0]*p[0]; }

    template<typename Value_t>
    Value_t userDefFuncSub(const Value_t* p) { return p[0]-p[1]; }

    template<typename Value_t>
    Value_t userDefFuncValue(const Value_t*) { return 10; }


    template<typename Value_t>
    class UserDefFuncWrapper:
        public FunctionParserBase<Value_t>::FunctionWrapper
    {
        Value_t (*mFuncPtr)(const Value_t*);
        unsigned mCounter;

     public:
        UserDefFuncWrapper(Value_t (*funcPtr)(const Value_t*)) :
            mFuncPtr(funcPtr), mCounter(0)
        {}

        virtual Value_t callFunction(const Value_t* values)
        {
            ++mCounter;
            return mFuncPtr(values);
        }

        unsigned counter() const { return mCounter; }
    };


    template<typename Value_t>
    inline Value_t Epsilon() { return Value_t(1e-9); }

#ifdef FP_SUPPORT_FLOAT_TYPE
    template<>
    inline float Epsilon<float>() { return 1e-3f; }
#endif

#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
    template<>
    inline long double Epsilon<long double>() { return 1e-10l; }
#endif

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    template<>
    inline MpfrFloat Epsilon<MpfrFloat>()
    {
        static const MpfrFloat eps(2e-20);
        return eps;
    }
#endif

#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
    template<>
    inline std::complex<float> Epsilon<std::complex<float> >()
        { return Epsilon<float>(); }
#endif

#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
    template<>
    inline std::complex<long double> Epsilon<std::complex<long double> >()
        { return Epsilon<long double>(); }
#endif


#ifndef _MSC_VER
    void setAnsiColor(unsigned color)
    {
        static int bold = 0;
        std::cout << "\33[";
        if(color > 7)
        {
            if(!bold) { std::cout << "1;"; bold=1; }
            color -= 7;
        }
        else if(bold) { std::cout << "0;"; bold=0; }
        std::cout << 30+color << "m";
    }

    void setAnsiBold() { std::cout << "\33[1m"; }

    void resetAnsiColor() { std::cout << "\33[0m"; }
#else
    void setAnsiColor(unsigned) {}
    void setAnsiBold() {}
    void resetAnsiColor() {}
#endif
}


//=========================================================================
// Copying testing functions
//=========================================================================
bool TestCopyingNoDeepCopy(DefaultParser p)
{
    DefaultValue_t vars[2] = { 3, 5 };

    if(std::fabs(p.Eval(vars) - 13) > Epsilon<DefaultValue_t>())
    {
        if(verbosityLevel >= 2)
        {
            std::cout
                << "\n - Giving as function parameter (no deep copy) failed."
                << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p.PrintByteCode(std::cout);
#endif
        }
        return false;
    }
    return true;
}

bool TestCopyingDeepCopy(DefaultParser p)
{
    DefaultValue_t vars[2] = { 3, 5 };

    p.Parse("x*y-1", "x,y");

    if(std::fabs(p.Eval(vars) - 14) > Epsilon<DefaultValue_t>())
    {
        if(verbosityLevel >= 2)
        {
            std::cout
                << "\n - Giving as function parameter (deep copy) failed."
                << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p.PrintByteCode(std::cout);
#endif
        }
        return false;
    }
    return true;
}

int TestCopying()
{
    bool retval = true;
    DefaultValue_t vars[2] = { 2, 5 };
    const DefaultValue_t epsilon = Epsilon<DefaultValue_t>();

    DefaultParser p1, p3;
    p1.Parse("x*y-2", "x,y");

    DefaultParser p2(p1);
    if(std::fabs(p2.Eval(vars) - 8) > epsilon)
    {
        retval = false;
        if(verbosityLevel >= 2)
        {
            std::cout << "\n - Copy constructor with no deep copy failed."
                      << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p2.PrintByteCode(std::cout);
#endif
        }
    }

    p2.Parse("x*y-1", "x,y");
    if(std::fabs(p2.Eval(vars) - 9) > epsilon)
    {
        retval = false;
        if(verbosityLevel >= 2)
        {
            std::cout << "\n - Copy constructor with deep copy failed."
                      << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p2.PrintByteCode(std::cout);
#endif
        }
    }

    p3 = p1;
    if(std::fabs(p3.Eval(vars) - 8) > epsilon)
    {
        retval = false;
        if(verbosityLevel >= 2)
        {
            std::cout << "\n - Assignment with no deep copy failed."
                      << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p3.PrintByteCode(std::cout);
#endif
        }
    }

    p3.Parse("x*y-1", "x,y");
    if(std::fabs(p3.Eval(vars) - 9) > epsilon)
    {
        retval = false;
        if(verbosityLevel >= 2)
        {
            std::cout << "\n - Assignment with deep copy failed."
                      << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p3.PrintByteCode(std::cout);
#endif
        }
    }

    if(!TestCopyingNoDeepCopy(p1))
        retval = false;

    // Final test to check that p1 still works:
    if(std::fabs(p1.Eval(vars) - 8) > epsilon)
    {
        retval = false;
        if(verbosityLevel >= 2)
            std::cout << "\n - Failed: p1 was corrupted." << std::endl;
    }

    if(!TestCopyingDeepCopy(p1))
        retval = false;

    // Final test to check that p1 still works:
    if(std::fabs(p1.Eval(vars) - 8) > epsilon)
    {
        retval = false;
        if(verbosityLevel >= 2)
        {
            std::cout << "\n - Failed: p1 was corrupted." << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
            p1.PrintByteCode(std::cout);
#endif
        }
    }

    return retval;
}


//=========================================================================
// Test error situations
//=========================================================================
int TestErrorSituations()
{
    bool retval = true;
    DefaultParser fp, tmpfp;
    fp.AddUnit("unit", 2);
    fp.AddFunction("Value", userDefFuncValue<DefaultValue_t>, 0);
    fp.AddFunction("Sqr", userDefFuncSqr<DefaultValue_t>, 1);
    fp.AddFunctionWrapper("Sub", UserDefFuncWrapper<DefaultValue_t>
                          (userDefFuncSub<DefaultValue_t>), 2);
    tmpfp.Parse("0", "x");

    static const struct
    {
        DefaultParser::ParseErrorType expected_error;
        int                            expected_error_position;
        const char*                    function_string;
    } invalidFuncs[] =
    {
      { DefaultParser::MISSING_PARENTH,     5, "sin(x"},
      { DefaultParser::EXPECT_PARENTH_FUNC, 4, "sin x"},
      { DefaultParser::SYNTAX_ERROR,        2, "x+" },
      { DefaultParser::EXPECT_OPERATOR,     2, "x x"},
      { DefaultParser::UNKNOWN_IDENTIFIER,  4, "sin(y)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   5, "sin(x, 1)" },
      { DefaultParser::EXPECT_OPERATOR,     1, "x, x"},
      { DefaultParser::SYNTAX_ERROR,        2, "x^^2" },
      { DefaultParser::SYNTAX_ERROR,        2, "x**x" },
      { DefaultParser::SYNTAX_ERROR,        2, "x+*x" },
      { DefaultParser::SYNTAX_ERROR,        0, "unit" },
      { DefaultParser::SYNTAX_ERROR,        0, "unit x" },
      { DefaultParser::SYNTAX_ERROR,        2, "x*unit" },
      { DefaultParser::SYNTAX_ERROR,        0, "unit*unit" },
      { DefaultParser::SYNTAX_ERROR,        0, "unit unit" },
      { DefaultParser::EXPECT_OPERATOR,     1, "x(unit)"},
      { DefaultParser::SYNTAX_ERROR,        2, "x+unit" },
      { DefaultParser::SYNTAX_ERROR,        2, "x*unit" },
      { DefaultParser::EMPTY_PARENTH,       1, "()"},
      { DefaultParser::SYNTAX_ERROR,        0, "" },
      { DefaultParser::EXPECT_OPERATOR,     1, "x()"},
      { DefaultParser::EMPTY_PARENTH,       3, "x*()"},
      { DefaultParser::SYNTAX_ERROR,        4, "sin(unit)" },
      { DefaultParser::EXPECT_PARENTH_FUNC, 4, "sin unit"},
      { DefaultParser::EXPECT_OPERATOR,     2, "1..2"},
      { DefaultParser::SYNTAX_ERROR,        1, "(" },
      { DefaultParser::MISM_PARENTH,        0, ")"},
      { DefaultParser::MISSING_PARENTH,     2, "(x"},
      { DefaultParser::EXPECT_OPERATOR,     1, "x)"},
      { DefaultParser::MISM_PARENTH,        0, ")x("},
      { DefaultParser::MISSING_PARENTH,     14,"(((((((x))))))"},
      { DefaultParser::EXPECT_OPERATOR,     15,"(((((((x))))))))"},
      { DefaultParser::EXPECT_OPERATOR,     1, "2x"},
      { DefaultParser::EXPECT_OPERATOR,     3, "(2)x"},
      { DefaultParser::EXPECT_OPERATOR,     3, "(x)2"},
      { DefaultParser::EXPECT_OPERATOR,     1, "2(x)"},
      { DefaultParser::EXPECT_OPERATOR,     1, "x(2)"},
      { DefaultParser::SYNTAX_ERROR,        0, "[x]" },
      { DefaultParser::SYNTAX_ERROR,        0, "@x" },
      { DefaultParser::SYNTAX_ERROR,        0, "$x" },
      { DefaultParser::SYNTAX_ERROR,        0, "{x}" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   5, "max(x)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   8, "max(x, 1, 2)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   6, "if(x,2)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   10,"if(x, 2, 3, 4)" },
      { DefaultParser::MISSING_PARENTH,     6, "Value(x)"},
      { DefaultParser::MISSING_PARENTH,     6, "Value(1+x)"},
      { DefaultParser::MISSING_PARENTH,     6, "Value(1,x)"},
      // Note: ^should these three not return ILL_PARAMS_AMOUNT instead?
      { DefaultParser::ILL_PARAMS_AMOUNT,   4, "Sqr()"},
      { DefaultParser::ILL_PARAMS_AMOUNT,   5, "Sqr(x,1)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   5, "Sqr(1,2,x)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   4, "Sub()" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   5, "Sub(x)" },
      { DefaultParser::ILL_PARAMS_AMOUNT,   7, "Sub(x,1,2)" },
      { DefaultParser::UNKNOWN_IDENTIFIER,  2, "x+Sin(1)" },
      { DefaultParser::UNKNOWN_IDENTIFIER,  0, "sub(1,2)" },
      { DefaultParser::UNKNOWN_IDENTIFIER,  0, "sinx(1)"  },
      { DefaultParser::UNKNOWN_IDENTIFIER,  2, "1+X"      },
      { DefaultParser::UNKNOWN_IDENTIFIER,  0, "eval(x)" }
    };
    const unsigned amnt = sizeof(invalidFuncs)/sizeof(invalidFuncs[0]);
    for(unsigned i = 0; i < amnt; ++i)
    {
        int parse_result = fp.Parse(invalidFuncs[i].function_string, "x");
        if(parse_result < 0)
        {
            retval = false;
            if(verbosityLevel >= 2)
            {
                std::cout << "\n - Parsing the invalid function \""
                          << invalidFuncs[i].function_string
                          << "\" didn't fail\n";
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
                fp.PrintByteCode(std::cout);
#endif
            }
        }
        else if(fp.GetParseErrorType() != invalidFuncs[i].expected_error
             || parse_result != invalidFuncs[i].expected_error_position)
        {
            retval = false;
            if(verbosityLevel >= 2)
            {
                std::cout << "\n - Parsing the invalid function \""
                          << invalidFuncs[i].function_string
                          << "\" produced ";
                if(fp.GetParseErrorType() != invalidFuncs[i].expected_error)
                    std::cout << "wrong error code (" << fp.ErrorMsg() << ")";
                if(parse_result != invalidFuncs[i].expected_error_position)
                    std::cout << "wrong pointer (expected "
                              << invalidFuncs[i].expected_error_position
                              << ", got " << parse_result << ")";
                std::cout << "\n";
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
                fp.PrintByteCode(std::cout);
#endif
            }
        }
    }

    static const char* const invalidNames[] =
    { "s2%", "sin", "(x)", "5x", "2", "\302\240"/*nbsp*/ };
    const unsigned namesAmnt = sizeof(invalidNames)/sizeof(invalidNames[0]);

    for(unsigned i = 0; i < namesAmnt; ++i)
    {
        const char* const n = invalidNames[i];
        if(fp.AddConstant(n, 1))
        {
            retval = false;
            if(verbosityLevel >= 2)
                std::cout << "\n - Adding an invalid name (\"" << n
                          << "\") as constant didn't fail" << std::endl;
        }
        if(fp.AddFunction(n, userDefFuncSqr<DefaultValue_t>, 1))
        {
            retval = false;
            if(verbosityLevel >= 2)
                std::cout << "\n - Adding an invalid name (\"" << n
                          << "\") as funcptr didn't fail" << std::endl;
        }
        if(fp.AddFunction(n, tmpfp))
        {
            retval = false;
            if(verbosityLevel >= 2)
                std::cout << "\n - Adding an invalid name (\"" << n
                          << "\") as funcparser didn't fail" << std::endl;
        }
        if(fp.Parse("0", n) < 0)
        {
            retval = false;
            if(verbosityLevel >= 2)
                std::cout << "\n - Using an invalid name (\"" << n
                          << "\") as variable name didn't fail" << std::endl;
        }
    }

    fp.AddConstant("CONST", 1);
    fp.AddFunction("PTR", userDefFuncSqr<DefaultValue_t>, 1);
    fp.AddFunction("PARSER", tmpfp);

    if(fp.AddConstant("PTR", 1))
    {
        retval = false;
        if(verbosityLevel >= 2)
            std::cout << "\n - Adding a userdef function (\"PTR\") as "
                      << "constant didn't fail" << std::endl;
    }
    if(fp.AddFunction("CONST", userDefFuncSqr<DefaultValue_t>, 1))
    {
        retval = false;
        if(verbosityLevel >= 2)
            std::cout << "\n - Adding a userdef constant (\"CONST\") as "
                      << "funcptr didn't fail" << std::endl;
    }
    if(fp.AddFunction("CONST", tmpfp))
    {
        retval = false;
        if(verbosityLevel >= 2)
            std::cout << "\n - Adding a userdef constant (\"CONST\") as "
                      << "funcparser didn't fail" << std::endl;
    }

    return retval;
}


//=========================================================================
// Thoroughly test whitespaces
//=========================================================================
DefaultValue_t wsFunc(DefaultValue_t x)
{
    return
        x + std::sin((x*-1.5)-(.5*2.0)*(((-x)*1.5+(2-(x)*2.0)*2.0)+(3.0*2.0))+
                     (1.5*2.0))+(cos(x)*2.0);
}

bool testWsFunc(DefaultParser& fp, const std::string& function)
{
    int res = fp.Parse(function, "x");
    if(res > -1)
    {
        if(verbosityLevel >= 2)
            std::cout << "\n - Parsing function:\n\"" << function
                      << "\"\nfailed at char " << res
                      << ": " << fp.ErrorMsg() << std::endl;
        return false;
    }

    DefaultValue_t vars[1];
    for(vars[0] = -2.0; vars[0] <= 2.0; vars[0] += .1)
        if(std::fabs(fp.Eval(vars) - wsFunc(vars[0])) >
           Epsilon<DefaultValue_t>())
        {
            return false;
        }
    return true;
}

int WhiteSpaceTest()
{
    DefaultParser fp;
    fp.AddConstant("const", 1.5);
    fp.AddUnit("unit", 2.0);
    std::string function(" x + sin ( ( x * - 1.5 ) - .5 unit * ( ( ( - x ) * "
                         "const + ( 2 - ( x ) unit ) unit ) + 3 unit ) + "
                         "( const ) unit ) + cos ( x ) unit ");

    if(!testWsFunc(fp, function)) return false;

    static const unsigned char WhiteSpaceTables[][4] =
    {
        { 1, 0x09, 0,0 }, // tab
        { 1, 0x0A, 0,0 }, // linefeed
        { 1, 0x0B, 0,0 }, // vertical tab
        { 1, 0x0D, 0,0 }, // carriage return
        { 1, 0x20, 0,0 }, // space
        { 2, 0xC2,0xA0, 0 }, // U+00A0 (nbsp)
        { 3, 0xE2,0x80,0x80 }, { 3, 0xE2,0x80,0x81 }, // U+2000 to...
        { 3, 0xE2,0x80,0x82 }, { 3, 0xE2,0x80,0x83 }, { 3, 0xE2,0x80,0x84 },
        { 3, 0xE2,0x80,0x85 }, { 3, 0xE2,0x80,0x86 }, { 3, 0xE2,0x80,0x87 },
        { 3, 0xE2,0x80,0x88 }, { 3, 0xE2,0x80,0x89 },
        { 3, 0xE2,0x80,0x8A }, { 3, 0xE2,0x80,0x8B }, // ... U+200B
        { 3, 0xE2,0x80,0xAF }, { 3, 0xE2,0x81,0x9F }, // U+202F and U+205F
        { 3, 0xE3,0x80,0x80 } // U+3000
    };
    const unsigned n_whitespaces =
        sizeof(WhiteSpaceTables)/sizeof(*WhiteSpaceTables);

    for(unsigned i = 0; i < function.size(); ++i)
    {
        if(function[i] == ' ')
        {
            function.erase(i, 1);
            for(std::size_t a = 0; a < n_whitespaces; ++a)
            {
                if(!testWsFunc(fp, function)) return false;
                int length = (int)WhiteSpaceTables[a][0];
                const char* sequence = (const char*)&WhiteSpaceTables[a][1];
                function.insert(i, sequence, length);
                if(!testWsFunc(fp, function)) return false;
                function.erase(i, length);
            }
        }
    }
    return true;
}


//=========================================================================
// Test integer powers
//=========================================================================
bool compareExpValues(DefaultValue_t value,
                      const std::string& funcStr,
                      DefaultValue_t v1,
                      DefaultValue_t v2,
                      bool isOptimized)
{
    const DefaultValue_t epsilon = Epsilon<DefaultValue_t>();
    const DefaultValue_t diff =
        std::fabs(v1) < epsilon ?
        (std::fabs(v2) < epsilon ? std::fabs(v1 - v2) :
         std::fabs((v1 - v2) / v2)) :
        std::fabs((v1 - v2) / v1);
    if(diff > epsilon)
    {
        if(verbosityLevel >= 2)
        {
            std::cout << "\n - For \"" << funcStr << "\" with x=" << value
                      << " the library (";
            if(!isOptimized) std::cout << "not ";
            std::cout << "optimized) returned\n"
                      << std::setprecision(18) << v2
                      << " instead of " << v1 << std::endl;
        }
        return false;
    }
    return true;
}

bool runIntPowTest(DefaultParser& fp, const std::string& funcStr,
                   int exponent, bool isOptimized)
{
    const int absExponent = exponent < 0 ? -exponent : exponent;

    for(int valueOffset = 0; valueOffset <= 5; ++valueOffset)
    {
        const DefaultValue_t value =
            (exponent >= 0 && valueOffset == 0) ? 0.0 :
            1.0+(valueOffset-1)/100.0;
        DefaultValue_t v1 = exponent == 0 ? 1 : value;
        for(int i = 2; i <= absExponent; ++i)
            v1 *= value;
        if(exponent < 0) v1 = 1.0/v1;

        const DefaultValue_t v2 = fp.Eval(&value);

        if(!compareExpValues(value, funcStr, v1, v2, isOptimized))
            return false;
    }

    return true;
}

bool runFractionalPowTest(const std::string& funcStr, double exponent)
{
    DefaultParser fp;
    if(fp.Parse(funcStr, "x") != -1)
    {
        if(verbosityLevel >= 2)
            std::cout << "\n - Parsing \"" << funcStr <<"\" failed: "
                      << fp.ErrorMsg() << "\n";
        return false;
    }

    for(int i = 0; i < 3; ++i)
    {
        for(int valueOffset = 0; valueOffset <= 10; ++valueOffset)
        {
            const DefaultValue_t value =
                (exponent >= 0 && valueOffset == 0) ? 0.0 :
                1.0+(valueOffset-1)/2.0;
            const DefaultValue_t v1 = std::pow(value, exponent);
            const DefaultValue_t v2 = fp.Eval(&value);

            if(!compareExpValues(value, funcStr, v1, v2, i > 0))
                return false;
        }
        fp.Optimize();
    }

    return true;
}

int TestIntPow()
{
    DefaultParser fp;

    for(int exponent = -1300; exponent <= 1300; ++exponent)
    {
        std::ostringstream os;
        os << "x^" << exponent;
        const std::string func = os.str();
        if(fp.Parse(func, "x") != -1)
        {
            if(verbosityLevel >= 2)
                std::cout << "\n - Parsing \"" << func <<"\" failed: "
                          << fp.ErrorMsg() << "\n";
            return false;
        }

        if(!runIntPowTest(fp, func, exponent, false)) return false;
        fp.Optimize();
        if(!runIntPowTest(fp, func, exponent, true)) return false;
    }

    for(int m = -27; m <= 27; ++m)
    {
        for(int n_sqrt=0; n_sqrt<=4; ++n_sqrt)
            for(int n_cbrt=0; n_cbrt<=4; ++n_cbrt)
            {
                if(n_sqrt+n_cbrt == 0) continue;

                std::ostringstream os;
                os << "x^(" << m << "/(1";
                for(int n=0; n<n_sqrt; ++n) os << "*2";
                for(int n=0; n<n_cbrt; ++n) os << "*3";
                os << "))";
                DefaultValue_t exponent = DefaultValue_t(m);
                if(n_sqrt > 0) exponent /= std::pow(2.0, n_sqrt);
                if(n_cbrt > 0) exponent /= std::pow(3.0, n_cbrt);
                if(!runFractionalPowTest(os.str(), exponent)) return false;
            }
    }

    return true;
}


//=========================================================================
// Test UTF-8 parsing
//=========================================================================
namespace
{
    typedef unsigned char UChar;
    struct CharValueRange { const UChar first, last; };

    const CharValueRange validValueRanges[][4] =
    {
        { { 0x30, 0x39 }, { 0, 0 }, { 0, 0 }, { 0, 0 } }, // digits
        { { 0x41, 0x5A }, { 0, 0 }, { 0, 0 }, { 0, 0 } }, // uppercase ascii
        { { 0x5F, 0x5F }, { 0, 0 }, { 0, 0 }, { 0, 0 } }, // underscore
        { { 0x61, 0x7A }, { 0, 0 }, { 0, 0 }, { 0, 0 } }, // lowercase ascii
        // U+0080 through U+009F
        { { 0xC2, 0xC2 }, { 0x80, 0x9F }, { 0, 0 }, { 0, 0 } },
        // U+00A1 through U+00BF
        { { 0xC2, 0xC2 }, { 0xA1, 0xBF }, { 0, 0 }, { 0, 0 } },
        // U+00C0 through U+07FF
        { { 0xC3, 0xDF }, { 0x80, 0xBF }, { 0, 0 }, { 0, 0 } },
        // U+0800 through U+1FFF (skip U+2000..U+200bB, which are whitespaces)
        { { 0xE0, 0xE0 }, { 0xA0, 0xBF }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xE1, 0xE1 }, { 0x80, 0xBF }, { 0x80, 0xBF }, { 0, 0 } },
        // U+200C through U+202E (skip U+202F, which is a whitespace)
        { { 0xE2, 0xE2 }, { 0x80, 0x80 }, { 0x8C, 0xAE }, { 0, 0 } },
        // U+2030 through U+205E (skip U+205F, which is a whitespace)
        { { 0xE2, 0xE2 }, { 0x80, 0x80 }, { 0xB0, 0xBF }, { 0, 0 } },
        { { 0xE2, 0xE2 }, { 0x81, 0x81 }, { 0x80, 0x9E }, { 0, 0 } },
        // U+2060 through U+20FF (skip U+3000, which is a whitespace)
        { { 0xE2, 0xE2 }, { 0x81, 0x81 }, { 0xA0, 0xBF }, { 0, 0 } },
        { { 0xE2, 0xE2 }, { 0x82, 0xBF }, { 0x80, 0xBF }, { 0, 0 } },
        // U+3001 through U+CFFF
        { { 0xE3, 0xE3 }, { 0x80, 0x80 }, { 0x81, 0xBF }, { 0, 0 } },
        { { 0xE3, 0xE3 }, { 0x81, 0xBF }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xE4, 0xEC }, { 0x80, 0xBF }, { 0x80, 0xBF }, { 0, 0 } },
        // U+E000 through U+FFFF
        { { 0xEE, 0xEF }, { 0x80, 0xBF }, { 0x80, 0xBF }, { 0, 0 } },
        // U+10000 through U+FFFFF
        { { 0xF0, 0xF0 }, { 0x90, 0xBF }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF1, 0xF3 }, { 0x80, 0xBF }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        // U+100000 through U+10FFFF
        { { 0xF4, 0xF4 }, { 0x80, 0x8F }, { 0x80, 0xBF }, { 0x80, 0xBF } }
    };
    const unsigned validValueRangesAmount =
        sizeof(validValueRanges)/sizeof(validValueRanges[0]);

    const CharValueRange invalidValueRanges[][4] =
    {
        // spaces:
        { { 0x09, 0x09 }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x0A, 0x0A }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x0B, 0x0B }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x0D, 0x0D }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x20, 0x20 }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0xC2, 0xC2 }, { 0xA0, 0xA0 }, { 0, 0 }, { 0, 0 } },
        { { 0xE2, 0xE2 }, { 0x80, 0x80 }, { 0x80, 0x8B }, { 0, 0 } },
        { { 0xE2, 0xE2 }, { 0x80, 0x80 }, { 0xAF, 0xAF }, { 0, 0 } },
        { { 0xE2, 0xE2 }, { 0x81, 0x81 }, { 0x9F, 0x9F }, { 0, 0 } },
        { { 0xE3, 0xE3 }, { 0x80, 0x80 }, { 0x80, 0x80 }, { 0, 0 } },
        // others:
        { { 0xC0, 0xC1 }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0xED, 0xED }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0xF5, 0xFF }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x21, 0x2F }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x3A, 0x40 }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x5B, 0x5E }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x60, 0x60 }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x7B, 0x7F }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0x80, 0xFF }, { 0, 0 }, { 0, 0 }, { 0, 0 } },
        { { 0xE0, 0xEF }, { 0x80, 0xFF }, { 0, 0 }, { 0, 0 } },
        { { 0xF0, 0xF4 }, { 0x80, 0xFF }, { 0x80, 0xFF }, { 0, 0 } },

        { { 0xC2, 0xDF }, { 0x00, 0x7F }, { 0, 0 }, { 0, 0 } },
        { { 0xC2, 0xDF }, { 0xC0, 0xFF }, { 0, 0 }, { 0, 0 } },

        { { 0xE0, 0xE0 }, { 0x00, 0x9F }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xE0, 0xE0 }, { 0xA0, 0xBF }, { 0x00, 0x7F }, { 0, 0 } },
        { { 0xE0, 0xE0 }, { 0xA0, 0xBF }, { 0xC0, 0xFF }, { 0, 0 } },

        { { 0xE1, 0xEC }, { 0x00, 0x7F }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xE1, 0xEC }, { 0xC0, 0xFF }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xE1, 0xEC }, { 0x80, 0xBF }, { 0x00, 0x7F }, { 0, 0 } },
        { { 0xE1, 0xEC }, { 0x80, 0xBF }, { 0xC0, 0xFF }, { 0, 0 } },

        { { 0xEE, 0xEF }, { 0x00, 0x7F }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xEE, 0xEF }, { 0xC0, 0xFF }, { 0x80, 0xBF }, { 0, 0 } },
        { { 0xEE, 0xEF }, { 0x80, 0xBF }, { 0x00, 0x7F }, { 0, 0 } },
        { { 0xEE, 0xEF }, { 0x80, 0xBF }, { 0xC0, 0xFF }, { 0, 0 } },

        { { 0xF0, 0xF0 }, { 0x00, 0x8F }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF0, 0xF0 }, { 0xC0, 0xFF }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF0, 0xF0 }, { 0x90, 0xBF }, { 0x00, 0x7F }, { 0x80, 0xBF } },
        { { 0xF0, 0xF0 }, { 0x90, 0xBF }, { 0xC0, 0xFF }, { 0x80, 0xBF } },
        { { 0xF0, 0xF0 }, { 0x90, 0xBF }, { 0x80, 0xBF }, { 0x00, 0x7F } },
        { { 0xF0, 0xF0 }, { 0x90, 0xBF }, { 0x80, 0xBF }, { 0xC0, 0xFF } },

        { { 0xF1, 0xF3 }, { 0x00, 0x7F }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF1, 0xF3 }, { 0xC0, 0xFF }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF1, 0xF3 }, { 0x80, 0xBF }, { 0x00, 0x7F }, { 0x80, 0xBF } },
        { { 0xF1, 0xF3 }, { 0x80, 0xBF }, { 0xC0, 0xFF }, { 0x80, 0xBF } },
        { { 0xF1, 0xF3 }, { 0x80, 0xBF }, { 0x80, 0xBF }, { 0x00, 0x7F } },
        { { 0xF1, 0xF3 }, { 0x80, 0xBF }, { 0x80, 0xBF }, { 0xC0, 0xFF } },

        { { 0xF4, 0xF4 }, { 0x00, 0x7F }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF4, 0xF4 }, { 0x90, 0xFF }, { 0x80, 0xBF }, { 0x80, 0xBF } },
        { { 0xF4, 0xF4 }, { 0x80, 0x8F }, { 0x00, 0x7F }, { 0x80, 0xBF } },
        { { 0xF4, 0xF4 }, { 0x80, 0x8F }, { 0xC0, 0xFF }, { 0x80, 0xBF } },
        { { 0xF4, 0xF4 }, { 0x80, 0x8F }, { 0x80, 0xBF }, { 0x00, 0x7F } },
        { { 0xF4, 0xF4 }, { 0x80, 0x8F }, { 0x80, 0xBF }, { 0xC0, 0xFF } }
    };
    const unsigned invalidValueRangesAmount =
        sizeof(invalidValueRanges)/sizeof(invalidValueRanges[0]);

    class CharIter
    {
        const CharValueRange (*valueRanges)[4];
        const unsigned valueRangesAmount;
        UChar charValues[4];
        unsigned rangeIndex, firstRangeIndex, skipIndex;

        void initCharValues()
        {
            for(unsigned i = 0; i < 4; ++i)
                charValues[i] = valueRanges[rangeIndex][i].first;
        }

     public:
        CharIter(bool skipDigits, bool skipLowerCaseAscii):
            valueRanges(validValueRanges),
            valueRangesAmount(validValueRangesAmount),
            rangeIndex(skipDigits ? 1 : 0),
            firstRangeIndex(skipDigits ? 1 : 0),
            skipIndex(skipLowerCaseAscii ? 3 : ~0U)
        {
            initCharValues();
        }

        CharIter():
            valueRanges(invalidValueRanges),
            valueRangesAmount(invalidValueRangesAmount),
            rangeIndex(0), firstRangeIndex(0), skipIndex(~0U)
        {
            initCharValues();
        }

        void appendChar(std::string& dest) const
        {
            for(unsigned i = 0; i < 4; ++i)
            {
                if(charValues[i] == 0) break;
                dest += char(charValues[i]);
            }
        }

        bool next()
        {
            for(unsigned i = 0; i < 4; ++i)
            {
                if(charValues[i] < valueRanges[rangeIndex][i].last)
                {
                    ++charValues[i];
                    return true;
                }
            }
            if(++rangeIndex == skipIndex) ++rangeIndex;
            if(rangeIndex < valueRangesAmount)
            {
                initCharValues();
                return true;
            }
            rangeIndex = firstRangeIndex;
            initCharValues();
            return false;
        }

        void print() const
        {
            std::printf("{");
            for(unsigned i = 0; i < 4; ++i)
            {
                if(charValues[i] == 0) break;
                if(i > 0) std::printf(",");
                std::printf("%02X", unsigned(charValues[i]));
            }
            std::printf("}");
        }
    };

    bool printUTF8TestError(const char* testType,
                            const CharIter* iters, unsigned length,
                            const std::string& identifier)
    {
        if(verbosityLevel >= 2)
        {
            std::printf("\n - %s failed with identifier ", testType);
            for(unsigned i = 0; i < length; ++i)
                iters[i].print();
            std::printf(": \"%s\"\n", identifier.c_str());
        }
        return false;
    }

    bool printUTF8TestError2(const CharIter* iters, unsigned length)
    {
        if(verbosityLevel >= 2)
        {
            std::printf("\n - Parsing didn't fail with invalid identifier ");
            for(unsigned i = 0; i < length; ++i)
                iters[(length-1)-i].print();
            std::printf("\n");
        }
        return false;
    }
}

int UTF8Test()
{
    typedef DefaultParser::value_type Value_t;

    CharIter iters[4] =
        { CharIter(true, false),
          CharIter(false, true),
          CharIter(false, false),
          CharIter(false, false) };
    std::string identifier;
    DefaultParser fp;
    const Value_t value = 0.0;

    for(unsigned length = 1; length <= 4; ++length)
    {
        if(verbosityLevel >= 1)
            std::cout << " " << length << std::flush;
        bool cont = true;
        while(cont)
        {
            identifier.clear();
            for(unsigned i = 0; i < length; ++i)
                iters[i].appendChar(identifier);

            if(fp.Parse(identifier, identifier) >= 0)
                return printUTF8TestError("Parsing", iters, length, identifier);

            if(fp.Eval(&value) != 0.0)
                return printUTF8TestError("Evaluation", iters, length,
                                          identifier);

            cont = false;
            const unsigned step = (length == 1) ? 1 : length-1;
            for(unsigned i = 0; i < length; i += step)
                if(iters[i].next())
                {
                    cont = true;
                    break;
                }
        }
    }

    CharIter invalidIters[3] =
        { CharIter(), CharIter(true, false), CharIter() };
    // test 5: inv
    // test 6: inv + normal
    // test 7: normal + inv

    for(unsigned length = 1; length <= 3; ++length)
    {
        if(verbosityLevel >= 1)
            std::cout << " " << 4+length << std::flush;
        unsigned numchars = length < 3 ? length : 2;
        unsigned firstchar = length < 3 ? 0 : 1;
        bool cont = true;
        while(cont)
        {
            identifier.clear();
            identifier += 'a';
            for(unsigned i = 0; i < numchars; ++i)
                invalidIters[firstchar+i].appendChar(identifier);
            identifier += 'a';

            if(fp.Parse(identifier, identifier) < 0)
                return printUTF8TestError2(invalidIters, length);

            cont = false;
            for(unsigned i = 0; i < numchars; ++i)
                if(invalidIters[firstchar+i].next())
                {
                    cont = true;
                    break;
                }
        }
    }

    return true;
}


//=========================================================================
// Test identifier adding and removal
//=========================================================================
bool AddIdentifier(DefaultParser& fp, const std::string& name, int type)
{
    static DefaultParser anotherParser;
    static bool anotherParserInitialized = false;
    if(!anotherParserInitialized)
    {
        anotherParser.Parse("x", "x");
        anotherParserInitialized = true;
    }

    switch(type)
    {
      case 0: return fp.AddConstant(name, 123);
      case 1: return fp.AddUnit(name, 456);
      case 2: return fp.AddFunction(name, userDefFuncSqr<DefaultValue_t>, 1);
      case 3: return fp.AddFunction(name, anotherParser);
    }
    return false;
}

int TestIdentifiers()
{
    DefaultParser fParser;
    std::vector<std::string> identifierNames(26*26, std::string("AA"));

    unsigned nameInd = 0;
    for(int i1 = 0; i1 < 26; ++i1)
    {
        for(int i2 = 0; i2 < 26; ++i2)
        {
            identifierNames.at(nameInd)[0] = char('A' + i1);
            identifierNames[nameInd][1] = char('A' + i2);

            if(!AddIdentifier(fParser, identifierNames[nameInd], (i1+26*i2)%3))
            {
                if(verbosityLevel >= 2)
                    std::cout << "\n - Failed to add identifier '"
                              << identifierNames[nameInd] << "'\n";
                return false;
            }

            ++nameInd;
        }
    }

    std::random_shuffle(identifierNames.begin(), identifierNames.end());

    for(unsigned nameInd = 0; nameInd <= identifierNames.size(); ++nameInd)
    {
        for(unsigned removedInd = 0; removedInd < nameInd; ++removedInd)
        {
            if(!AddIdentifier(fParser, identifierNames[removedInd], 3))
            {
                if(verbosityLevel >= 2)
                    std::cout
                        << "\n - Failure: Identifier '"
                        << identifierNames[removedInd]
                        << "' was still reserved even after removing it.\n";
                return false;
            }
            if(!fParser.RemoveIdentifier(identifierNames[removedInd]))
            {
                if(verbosityLevel >= 2)
                    std::cout << "\n - Failure: Removing the identifier '"
                              << identifierNames[removedInd]
                              << "' after adding it again failed.\n";
                return false;
            }
        }

        for(unsigned existingInd = nameInd;
            existingInd < identifierNames.size(); ++existingInd)
        {
            if(AddIdentifier(fParser, identifierNames[existingInd], 3))
            {
                if(verbosityLevel >= 2)
                    std::cout << "\n - Failure: Trying to add identifier '"
                              << identifierNames[existingInd]
                              << "' for a second time didn't fail.\n";
                return false;
            }
        }

        if(nameInd < identifierNames.size())
        {
            if(!fParser.RemoveIdentifier(identifierNames[nameInd]))
            {
                if(verbosityLevel >= 2)
                    std::cout << "\n - Failure: Trying to remove identifier '"
                              << identifierNames[nameInd] << "' failed.\n";
                return false;
            }
            if(fParser.RemoveIdentifier(identifierNames[nameInd]))
            {
                if(verbosityLevel >= 2)
                    std::cout << "\n - Failure: Trying to remove identifier '"
                              << identifierNames[nameInd]
                              << "' for a second time didn't fail.\n";
                return false;
            }
        }
    }

    return true;
}


//=========================================================================
// Test user-defined functions
//=========================================================================
namespace
{
    template<int VarsAmount>
    DefaultValue_t userFunction(const DefaultValue_t* p)
    {
        DefaultValue_t result = 1.0;
        for(int i = 0; i < VarsAmount; ++i)
            result += (VarsAmount+i/10.0) * p[i];
        return result;
    }

    DefaultValue_t(*userFunctions[])(const DefaultValue_t*) =
    {
        userFunction<0>, userFunction<1>, userFunction<2>, userFunction<3>,
        userFunction<4>, userFunction<5>, userFunction<6>, userFunction<7>,
        userFunction<8>, userFunction<9>, userFunction<10>, userFunction<11>
    };
    const unsigned userFunctionsAmount =
        sizeof(userFunctions) / sizeof(userFunctions[0]);

    DefaultValue_t nestedFunc1(const DefaultValue_t* p)
    {
        return p[0] + 2.0*p[1] + 3.0*p[2];
    }

    DefaultValue_t nestedFunc2(const DefaultValue_t* p)
    {
        const DefaultValue_t params[3] = { -5.0*p[0], -10.0*p[1], -p[0] };
        return p[0] + 4.0*nestedFunc1(params);
    }

    DefaultValue_t nestedFunc3(const DefaultValue_t* p)
    {
        const DefaultValue_t params1[3] = { 2.5*p[0]+2.0, p[2], p[1]/2.5 };
        const DefaultValue_t params2[2] = { p[1] / 1.5 - 1.0, p[0] - 2.5 };
        return nestedFunc1(params1) + nestedFunc2(params2);
    }
}

int testUserDefinedFunctions()
{
    const DefaultValue_t epsilon = Epsilon<DefaultValue_t>();

    DefaultParser nestedParser1, nestedParser2, nestedParser3;
    nestedParser1.Parse("x + 2.0*y + 3.0*z", "x, y, z");
    nestedParser2.AddFunction("nestedFunc1", nestedParser1);
    nestedParser2.Parse("x + 4.0*nestedFunc1(-5.0*x, -10.0*y, -x)", "x,y");
    nestedParser3.AddFunction("nestedFunc1", nestedParser1);
    nestedParser3.AddFunction("nestedFunc2", nestedParser2);
    nestedParser3.Parse("nestedFunc1(2.5*x+2.0, z, y/2.5) + "
                        "nestedFunc2(y/1.5 - 1.0, x - 2.5)", "x,y,z");

    for(int iteration = 0; iteration < 2; ++iteration)
    {
        DefaultValue_t nestedFuncParams[3];
        for(int i = 0; i < 100; ++i)
        {
            nestedFuncParams[0] = -10.0 + 20.0*i/100.0;
            for(int j = 0; j < 100; ++j)
            {
                nestedFuncParams[1] = -10.0 + 20.0*j/100.0;
                for(int k = 0; k < 100; ++k)
                {
                    nestedFuncParams[2] = -10.0 + 20.0*k/100.0;

                    const DefaultValue_t v1 =
                        nestedParser3.Eval(nestedFuncParams);
                    const DefaultValue_t v2 =
                        nestedFunc3(nestedFuncParams);
                    if(std::fabs(v1-v2) > epsilon)
                    {
                        if(verbosityLevel >= 2)
                            std::cout
                                << "\n - Nested function test failed with "
                                << "parameter values ("
                                << nestedFuncParams[0] << ","
                                << nestedFuncParams[1]
                                << ").\nThe library "
                                << (iteration > 0 ? "(optimized) " : "")
                                << "returned " << v1
                                << " instead of " << v2 << "." << std::endl;
                        return false;
                    }
                }
            }
        }
        nestedParser3.Optimize();
    }

    std::string funcNames[userFunctionsAmount];
    std::string userFunctionParserFunctions[userFunctionsAmount];
    std::string userFunctionParserParameters[userFunctionsAmount];
    DefaultParser userFunctionParsers[userFunctionsAmount];
    DefaultValue_t funcParams[userFunctionsAmount];
    DefaultParser parser1, parser2;

    for(unsigned funcInd = 0; funcInd < userFunctionsAmount; ++funcInd)
    {
        std::ostringstream functionString, paramString;

        functionString << '1';
        for(unsigned paramInd = 0; paramInd < funcInd; ++paramInd)
        {
            functionString << "+" << funcInd+paramInd/10.0
                           << "*p" << paramInd;

            if(paramInd > 0) paramString << ',';
            paramString << "p" << paramInd;
        }

        userFunctionParserFunctions[funcInd] = functionString.str();
        userFunctionParserParameters[funcInd] = paramString.str();

        if(userFunctionParsers[funcInd].Parse
           (userFunctionParserFunctions[funcInd],
            userFunctionParserParameters[funcInd]) >= 0)
        {
            if(verbosityLevel >= 2)
                std::cout << "\n - Failed to parse function\n\""
                          << functionString.str() << "\"\nwith parameters: \""
                          << paramString.str() << "\":\n"
                          << userFunctionParsers[funcInd].ErrorMsg() << "\n";
            return false;
        }

        for(unsigned testInd = 0; testInd < 10; ++testInd)
        {
            for(unsigned paramInd = 0; paramInd < testInd; ++paramInd)
                funcParams[paramInd] = testInd+paramInd;
            const DefaultValue_t result = userFunctions[funcInd](funcParams);
            const DefaultValue_t parserResult =
                userFunctionParsers[funcInd].Eval(funcParams);
            if(std::fabs(result - parserResult) > epsilon)
            {
                if(verbosityLevel >= 2)
                {
                    std::cout << "\n - Function\n\"" << functionString.str()
                              << "\"\nwith parameters (";
                    for(unsigned paramInd = 0; paramInd < testInd; ++paramInd)
                    {
                        if(paramInd > 0) std::cout << ',';
                        std::cout << funcParams[paramInd];
                    }
                    std::cout << ")\nreturned " << parserResult
                              << " instead of " << result << "\n";
                }
                return false;
            }
        }
    }

    for(unsigned funcInd = 0; funcInd < userFunctionsAmount; ++funcInd)
    {
        funcNames[funcInd] = "func00";
        funcNames[funcInd][4] = char('0' + funcInd/10);
        funcNames[funcInd][5] = char('0' + funcInd%10);

        if(!parser1.AddFunction(funcNames[funcInd], userFunctions[funcInd],
                                funcInd))
        {
            if(verbosityLevel >= 2)
                std::cout << "\n - Failed to add user-defined function \""
                          << funcNames[funcInd] << "\".\n";
            return false;
        }
        if(!parser2.AddFunction(funcNames[funcInd],
                                userFunctionParsers[funcInd]))
        {
            if(verbosityLevel >= 2)
                std::cout
                    << "\n - Failed to add user-defined function parser \""
                    << funcNames[funcInd] << "\".\n";
            return false;
        }

        std::ostringstream functionString;
        for(unsigned factorInd = 0; factorInd <= funcInd; ++factorInd)
        {
            if(factorInd > 0) functionString << '+';
            functionString << factorInd+1 << "*"
                           << funcNames[factorInd] << '(';
            for(unsigned paramInd = 0; paramInd < factorInd; ++paramInd)
            {
                if(paramInd > 0) functionString << ',';
                const unsigned value = factorInd*funcInd + paramInd;
                functionString << value << "+x";
            }
            functionString << ')';
        }

        if(parser1.Parse(functionString.str(), "x") >= 0)
        {
            if(verbosityLevel >= 2)
                std::cout << "\n - parser1 failed to parse function\n\""
                          << functionString.str() << "\":\n"
                          << parser1.ErrorMsg() << "\n";
            return false;
        }
        if(parser2.Parse(functionString.str(), "x") >= 0)
        {
            if(verbosityLevel >= 2)
                std::cout << "\n - parser2 failed to parse function\n\""
                          << functionString.str() << "\":\n"
                          << parser2.ErrorMsg() << "\n";
            return false;
        }

        for(unsigned optimizeInd = 0; optimizeInd < 4; ++optimizeInd)
        {
            for(unsigned testInd = 0; testInd < 100; ++testInd)
            {
                const DefaultValue_t x = testInd/10.0;
                DefaultValue_t result = 0.0;
                for(unsigned factorInd = 0; factorInd <= funcInd; ++factorInd)
                {
                    for(unsigned paramInd = 0; paramInd < factorInd; ++paramInd)
                    {
                        const unsigned value = factorInd*funcInd + paramInd;
                        funcParams[paramInd] = value+x;
                    }
                    result +=
                        (factorInd+1) * userFunctions[factorInd](funcParams);
                }

                const DefaultValue_t parser1Result = parser1.Eval(&x);
                const DefaultValue_t parser2Result = parser2.Eval(&x);
                const bool parser1Failed =
                    std::fabs(result - parser1Result) > epsilon;
                const bool parser2Failed =
                    std::fabs(result - parser2Result) > epsilon;

                if(parser1Failed || parser2Failed)
                {
                    if(verbosityLevel >= 2)
                    {
                        std::cout << "\n - For function:\n\""
                                  << functionString.str() << "\"";
                        if(optimizeInd > 0)
                            std::cout << "\n(Optimized " << optimizeInd
                                      << (optimizeInd > 1 ?
                                          " times)" : " time)");
                        std::cout << "\nwith x=" << x
                                  << " parser";
                        if(parser1Failed)
                            std::cout << "1 returned " << parser1Result;
                        else
                            std::cout << "2 returned " << parser2Result;
                        std::cout << " instead of " << result << ".\n";

                        if(parser2Failed)
                        {
                            std::cout << "The user-defined functions are:\n";
                            for(unsigned i = 0; i <= funcInd; ++i)
                                std::cout << funcNames[i] << "=\""
                                          << userFunctionParserFunctions[i]
                                          << "\"\n";
                        }
                    }

                    return false;
                }
            }

            parser1.Optimize();
        }
    }

    return true;
}


//=========================================================================
// Multithreaded test
//=========================================================================
#if defined(FP_USE_THREAD_SAFE_EVAL) || \
    defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
#include <boost/thread.hpp>

class TestingThread
{
    int mThreadNumber;
    DefaultParser* mFp;
    volatile static bool mOk;

    static DefaultValue_t function(const DefaultValue_t* vars)
    {
        const DefaultValue_t x = vars[0], y = vars[1];
        return sin(sqrt(x*x+y*y)) + 2*cos(2*sqrt(2*x*x+2*y*y));
    }

 public:
    TestingThread(int n, DefaultParser* fp):
        mThreadNumber(n), mFp(fp)
    {}

    static bool ok() { return mOk; }

    void operator()()
    {
        const DefaultValue_t epsilon = Epsilon<DefaultValue_t>();
        DefaultValue_t vars[2];
        for(vars[0] = -10.0; vars[0] <= 10.0; vars[0] += 0.02)
        {
            for(vars[1] = -10.0; vars[1] <= 10.0; vars[1] += 0.02)
            {
                if(!mOk) return;

                const DefaultValue_t v1 = function(vars);
                const DefaultValue_t v2 = mFp->Eval(vars);
                /*
                const double scale = pow(10.0, floor(log10(fabs(v1))));
                const double sv1 = fabs(v1) < Epsilon<double>() ? 0 : v1/scale;
                const double sv2 = fabs(v2) < Epsilon<double>() ? 0 : v2/scale;
                const double diff = fabs(sv2-sv1);
                */
                const DefaultValue_t diff =
                    std::fabs(v1) < epsilon ?
                    (std::fabs(v2) < epsilon ?
                     std::fabs(v1 - v2) :
                     std::fabs((v1 - v2) / v2)) :
                    std::fabs((v1 - v2) / v1);

                if(std::fabs(diff) > epsilon)
                {
                    mOk = false;
                    if(verbosityLevel >= 2)
                        std::cout << "\n - Thread " << mThreadNumber
                                  << " failed ([" << vars[0] << "," << vars[1]
                                  << "] -> " << v2 << " vs. " << v1 << ")"
                                  << std::endl;
                    return;
                }
            }
        }
    }
};

volatile bool TestingThread::mOk = true;

int testMultithreadedEvaluation()
{
    DefaultParser fp;
    fp.Parse("sin(sqrt(x*x+y*y)) + 2*cos(2*sqrt(2*x*x+2*y*y))", "x,y");

    if(verbosityLevel >= 1)
        std::cout << " 1" << std::flush;
    boost::thread t1(TestingThread(1, &fp)), t2(TestingThread(2, &fp));
    t1.join();
    t2.join();
    if(!TestingThread::ok()) return false;

    if(verbosityLevel >= 1)
        std::cout << " 2" << std::flush;
    boost::thread
        t3(TestingThread(3, &fp)), t4(TestingThread(4, &fp)),
        t5(TestingThread(5, &fp)), t6(TestingThread(6, &fp));
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    if(!TestingThread::ok()) return false;

    if(verbosityLevel >= 1)
        std::cout << " 3" << std::flush;
    fp.Optimize();
    boost::thread
        t7(TestingThread(7, &fp)), t8(TestingThread(8, &fp)),
        t9(TestingThread(9, &fp));
    t7.join();
    t8.join();
    t9.join();
    if(!TestingThread::ok()) return false;

    return true;
}

#else

int testMultithreadedEvaluation()
{
    return -1;
}

#endif

//=========================================================================
// Test variable deduction
//=========================================================================
template<typename Value_t>
struct TestType
{
    unsigned paramAmount;
    Value_t paramMin, paramMax, paramStep;
    bool useDegrees;

    Value_t (*funcPtr)(const Value_t*);
    double (*doubleFuncPtr)(const double*);
    long (*longFuncPtr)(const long*);

    const char* paramString;
    const char* testName;
    const char* funcString;
};


template<typename Value_t>
bool checkVarString(const char* idString,
                    FunctionParserBase<Value_t> & fp,
                    const TestType<Value_t>& testData,
                    int errorIndex,
                    int variablesAmount, const std::string& variablesString,
                    std::ostream& briefErrorMessages)
{
    const bool stringsMatch =
        (variablesString == testData.paramString);
    if(errorIndex >= 0 ||
       variablesAmount != int(testData.paramAmount) ||
       !stringsMatch)
    {
        if(verbosityLevel >= 2)
        {
            std::cout << "\n" << idString
                      << " ParseAndDeduceVariables() failed with function:\n\""
                      << testData.funcString << "\"\n";
            if(errorIndex >= 0)
                std::cout << "Error index: " << errorIndex
                          << ": " << fp.ErrorMsg() << std::endl;
            else if(!stringsMatch)
                std::cout << "Deduced var string was \"" << variablesString
                          << "\" instead of \""
                          << testData.paramString
                          << "\"." << std::endl;
            else
                std::cout << "Deduced variables amount was "
                          << variablesAmount << " instead of "
                          << testData.paramAmount << "."
                          << std::endl;
        }
        else
        {
            briefErrorMessages << "- " << testData.testName
                               << ": Failed ParseAndDeduceVariables().\n";
        }
        return false;
    }
    return true;
}

template<typename Value_t>
bool testVariableDeduction(FunctionParserBase<Value_t>& fp,
                           const TestType<Value_t>& testData,
                           std::ostream& briefErrorMessages)
{
    static std::string variablesString;
    static std::vector<std::string> variables;

    if(verbosityLevel >= 3)
        std::cout << "(Variable deduction)" << std::flush;

    int variablesAmount = -1;
    int retval = fp.ParseAndDeduceVariables
        (testData.funcString,
         &variablesAmount, testData.useDegrees);
    if(retval >= 0 || variablesAmount !=
       int(testData.paramAmount))
    {
        if(verbosityLevel >= 2)
        {
            std::cout <<
                "\nFirst ParseAndDeduceVariables() failed with function:\n\""
                      << testData.funcString << "\"\n";
            if(retval >= 0)
                std::cout << "Error index: " << retval
                          << ": " << fp.ErrorMsg() << std::endl;
            else
                std::cout << "Deduced variables amount was "
                          << variablesAmount << " instead of "
                          << testData.paramAmount << "."
                          << std::endl;
        }
        else
        {
            briefErrorMessages << "- " << testData.testName
                               << ": Failed ParseAndDeduceVariables().\n";
        }
        return false;
    }

    variablesAmount = -1;
    retval = fp.ParseAndDeduceVariables
        (testData.funcString,
         variablesString,
         &variablesAmount,
         testData.useDegrees);
    if(!checkVarString("Second", fp, testData, retval, variablesAmount,
                       variablesString, briefErrorMessages))
        return false;

    retval = fp.ParseAndDeduceVariables(testData.funcString,
                                        variables,
                                        testData.useDegrees);
    variablesAmount = int(variables.size());
    variablesString.clear();
    for(unsigned i = 0; i < variables.size(); ++i)
    {
        if(i > 0) variablesString += ',';
        variablesString += variables[i];
    }
    return checkVarString("Third", fp, testData, retval, variablesAmount,
                          variablesString, briefErrorMessages);
}


//=========================================================================
// Main test function
//=========================================================================
namespace
{
    template<typename Value_t>
    struct RegressionTests
    {
        static const TestType<Value_t> Tests[];
    };
    template<typename Value_t>
    const TestType<Value_t> RegressionTests<Value_t>::Tests[] = { TestType<Value_t>() };

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    inline MpfrFloat makeMF(const char* str)
    {
        MpfrFloat f;
        f.parseValue(str);
        return f;
    }
#endif

    /* These functions in fparser produce bool values. However,
     * the testing functions require that they produce Value_t's. */
    #define BoolProxy(Fname) \
    template<typename Value_t> \
    Value_t tb_##Fname(const Value_t& a, const Value_t& b) \
        { return Value_t(FUNCTIONPARSERTYPES::Fname(a,b)); }

    BoolProxy(fp_less)
    BoolProxy(fp_lessOrEq)
    BoolProxy(fp_greater)
    BoolProxy(fp_greaterOrEq)
    BoolProxy(fp_equal)
    BoolProxy(fp_nequal)

    template<typename Value_t>
    Value_t fp_truth(const Value_t& a)
    { return Value_t(FUNCTIONPARSERTYPES::fp_truth(a)); }

// Maybe these should be used in the test files instead...
#define fp_less tb_fp_less
#define fp_lessOrEq tb_fp_lessOrEq
#define fp_greater tb_fp_greater
#define fp_greaterOrEq tb_fp_greaterOrEq
#define fp_equal tb_fp_equal
#define fp_nequal tb_fp_nequal

#if defined(FP_SUPPORT_GMP_INT_TYPE) && !defined(FP_SUPPORT_LONG_INT_TYPE)
#define FP_SUPPORT_LONG_INT_TYPE
#include "testbed_tests.inc"
#undef FP_SUPPORT_LONG_INT_TYPE
#else
#include "testbed_tests.inc"
#endif

#undef fp_less
#undef fp_lessOrEq
#undef fp_greater
#undef fp_greaterOrEq
#undef fp_equal
#undef fp_nequal
}

namespace
{
    template<typename Value_t>
    void testAgainstDouble(Value_t*, Value_t, const TestType<Value_t>&,
                           std::ostream&) {}

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    void testAgainstDouble(MpfrFloat* vars, MpfrFloat parserValue,
                           const TestType<MpfrFloat>& testData,
                           std::ostream& error)
    {
        if(!testData.doubleFuncPtr) return;

        double doubleVars[10];
        for(unsigned i = 0; i < 10; ++i) doubleVars[i] = vars[i].toDouble();

        const double Eps = Epsilon<double>();

        const double v1 = testData.doubleFuncPtr(doubleVars);
        const double v2 = parserValue.toDouble();

        /*
        using namespace FUNCTIONPARSERTYPES;
        const double scale = fp_pow(10.0, fp_floor(fp_log10(fp_abs(v1))));
        const double sv1 = fp_abs(v1) < Eps ? 0 : v1/scale;
        const double sv2 = fp_abs(v2) < Eps ? 0 : v2/scale;
        const double diff = fp_abs(sv2-sv1);
        */
        const double diff =
            std::fabs(v1) < Eps ?
            (std::fabs(v2) < Eps ? std::fabs(v1 - v2) :
             std::fabs((v1 - v2) / v2)) :
            std::fabs((v1 - v2) / v1);

        if(diff > Eps)
        {
            using namespace FUNCTIONPARSERTYPES;
            if(verbosityLevel >= 2)
                error << std::setprecision(16) << v2 << " instead of "
                      << std::setprecision(16) << v1
                      << "\n(Difference: "
                      << std::setprecision(16) << v2-v1
                      << ", epsilon: "
                      << std::setprecision(16) << Eps
                      << "; scaled diff "
                      << std::setprecision(16) << diff
                      << ")\nwhen tested against the double function.";
            else
                error << std::setprecision(16) << v2 << " vs "
                      << std::setprecision(16) << v1
                      << " (diff: "
                      << std::setprecision(16) << v2-v1
                      << ", sdiff "
                      << std::setprecision(16) << diff
                      << ") against double.";
        }
    }
#endif

    template<typename Value_t>
    void testAgainstLongInt(Value_t*, Value_t, const TestType<Value_t>&,
                            std::ostream&) {}

#ifdef FP_SUPPORT_GMP_INT_TYPE
    void testAgainstLongInt(GmpInt* vars, GmpInt parserValue,
                            const TestType<GmpInt>& testData,
                            std::ostream& error)
    {
        if(!testData.longFuncPtr) return;

        long longVars[10];
        for(unsigned i = 0; i < 10; ++i) longVars[i] = vars[i].toInt();

        const long longValue = testData.longFuncPtr(longVars);
        if(longValue != parserValue)
        {
            if(verbosityLevel >= 2)
                error << parserValue << " instead of " << longValue
                      << "\nwhen tested against the long int function.";
            else
                error << parserValue << " vs " << longValue
                      << " against long.";
        }
    }
#endif
}

template<typename Value_t>
bool runRegressionTest(FunctionParserBase<Value_t>& fp,
                       const TestType<Value_t>& testData,
                       const std::string& valueType,
                       const Value_t Eps,
                       std::ostream& briefErrorMessages)
{
    Value_t vars[10];
    Value_t fp_vars[10];

    for(unsigned i = 0; i < testData.paramAmount; ++i)
        vars[i] = testData.paramMin;

    while(true)
    {
        unsigned paramInd = 0;
        while(paramInd < testData.paramAmount)
        {
            using namespace FUNCTIONPARSERTYPES;
            /* ^ Import a possible <= operator from that
             *   namespace for this particular comparison only */
            vars[paramInd] += testData.paramStep;
            if(vars[paramInd] <= testData.paramMax) break;
            vars[paramInd++] = testData.paramMin;
        }

        if(paramInd == testData.paramAmount) break;

        for(unsigned i = 0; i < testData.paramAmount; ++i)
            fp_vars[i] = vars[i];

        if(verbosityLevel >= 4)
        {
            std::cout << "Trying (";
            for(unsigned ind = 0; ind < testData.paramAmount; ++ind)
                std::cout << (ind>0 ? ", " : "") << vars[ind];
            std::cout << ")\n" << std::flush;
        }
        const Value_t v1 = testData.funcPtr(vars);
        if(true) /*test Eval() */
        {
            const Value_t v2 = fp.Eval(fp_vars);

            std::ostringstream error;

            if(fp.EvalError() > 0)
            {
                error << "EvalError " << fp.EvalError() << " ("
                      << getEvalErrorName(fp.EvalError()) << ")";
            }
            else if(FUNCTIONPARSERTYPES::IsIntType<Value_t>::result)
            {
                if(v1 != v2)
                {
                    using namespace FUNCTIONPARSERTYPES;
                    if(verbosityLevel >= 2)
                        error << v2 << " instead of " << v1;
                    else
                        error << v2 << " vs " << v1;
                }
                else
                    testAgainstLongInt(vars, v2, testData, error);
            }
            else
            {
                using namespace FUNCTIONPARSERTYPES;
                /*
                const Value_t scale =
                    fp_pow(Value_t(10.0), fp_floor(fp_log10(fp_abs(v1))));
                const Value_t sv1 = fp_abs(v1) < Eps ? 0 : v1/scale;
                const Value_t sv2 = fp_abs(v2) < Eps ? 0 : v2/scale;
                const Value_t diff = fp_abs(sv2-sv1);
                */
                const Value_t diff =
                    fp_abs(v1) < Eps ?
                    (fp_abs(v2) < Eps ? fp_abs(v1 - v2) :
                     fp_abs((v1 - v2) / v2)) :
                    fp_abs((v1 - v2) / v1);
                /*
                const Value_t diff =
                    v1 == Value_t(0) ?
                    (v2 == Value_t(0) ? Value_t(0) :
                     fp_abs((v1 - v2) / v2)) :
                    fp_abs((v1 - v2) / v1);
                */

                if(diff > Eps)
                {
                    using namespace FUNCTIONPARSERTYPES;
                    if(verbosityLevel >= 2)
                        error << std::setprecision(28) << v2 << " instead of "
                              << std::setprecision(28) << v1
                              << "\n(Difference: "
                              << std::setprecision(28) << v2-v1
                              << ", epsilon: "
                              << std::setprecision(28) << Eps
                              << "; scaled diff "
                              << std::setprecision(28) << diff
                              << ")";
                    else
                        error << std::setprecision(16) << v2 << " vs "
                              << std::setprecision(16) << v1
                              << " (diff: "
                              << std::setprecision(16) << v2-v1
                              << ", sdiff "
                              << std::setprecision(16) << diff
                              << ")";
                }
                else
                    testAgainstDouble(vars, v2, testData, error);
            }

            if(!error.str().empty())
            {
                if(verbosityLevel == 2)
                    std::cout << "\n****************************\nTest "
                              << testData.testName
                              << ", function:\n\"" << testData.funcString
                              << "\"\n(" << valueType << ")";

                if(verbosityLevel >= 2)
                {
                    using namespace FUNCTIONPARSERTYPES;
                    // ^ For output of complex numbers according to fparser practice

                    std::cout << std::endl << "Error: For (" << std::setprecision(20);
                    for(unsigned ind = 0; ind < testData.paramAmount; ++ind)
                        std::cout << (ind>0 ? ", " : "") << vars[ind];
                    std::cout << ")\nthe library returned " << error.str()
                              << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
                    fp.PrintByteCode(std::cout);
#endif
                }
                else
                {
                    using namespace FUNCTIONPARSERTYPES;

                    briefErrorMessages << "- " << testData.testName << " (";
                    for(unsigned ind = 0; ind < testData.paramAmount; ++ind)
                        briefErrorMessages << (ind>0 ? "," : "") << vars[ind];
                    briefErrorMessages << "): " << error.str() << "\n";
                }
                return false;
            }
        } /* test Eval() */
    }
    return true;
}

static bool WildMatch(const char *pattern, const char *what)
{
    for(; *what || *pattern; ++what, ++pattern)
        if(*pattern == '*')
        {
            while(*++pattern == '*') {}
            for(; *what; ++what)
                if(WildMatch(pattern, what))
                    return true;
            return !*pattern;
        }
        else if(*pattern != '?' && *pattern != *what)
            return false;
    return true;
}
static bool WildMatch_Dirmask(const char *pattern, const char *what)
{
    std::string testmask = pattern;
    if(testmask.find('/') == testmask.npos) testmask = "*/" + testmask;
    return WildMatch(testmask.c_str(), what);
}
bool IsSelectedTest(const char* testName)
{
    for(std::size_t a=0; a<selectedRegressionTests.size(); ++a)
        if(WildMatch_Dirmask(selectedRegressionTests[a], testName))
            return true;
    return false;
}
/* Asciibetical comparator, with in-string integer values sorted naturally */
bool natcomp(const std::string& a, const std::string& b)
{
    std::size_t ap=0, bp=0;
    while(ap < a.size() && bp < b.size())
    {
        if(a[ap] >= '0' && a[ap] <= '9'
        && b[bp] >= '0' && b[bp] <= '9')
        {
            unsigned long aval = (a[ap++] - '0');
            unsigned long bval = (b[bp++] - '0');
            while(ap < a.size() && a[ap] >= '0' && a[ap] <= '9')
                aval = aval*10ul + (a[ap++] - '0');
            while(bp < b.size() && b[bp] >= '0' && b[bp] <= '9')
                bval = bval*10ul + (b[bp++] - '0');
            if(aval != bval)
                return aval < bval;
        }
        else
        {
            if(a[ap] != b[ap]) return a[ap] < b[ap];
            ++ap; ++bp;
        }
    }
    return (bp < b.size() && ap >= a.size());
}


template<typename Value_t>
bool runRegressionTests(const std::string& valueType)
{
    // Setup the function parser for testing
    // -------------------------------------
    FunctionParserBase<Value_t> fp;

    if(verbosityLevel >= 1)
    {
        setAnsiBold();
        std::cout << "==================== Parser type \"" << valueType
                  << "\" ====================" << std::endl;
        resetAnsiColor();
    }

    bool ret = fp.AddConstant("pi",
                              FUNCTIONPARSERTYPES::fp_const_pi<Value_t>());
    ret = ret && fp.AddConstant("naturalnumber",
                                FUNCTIONPARSERTYPES::fp_const_e<Value_t>());
    ret = ret && fp.AddConstant("logtwo",
                                FUNCTIONPARSERTYPES::fp_const_log2<Value_t>());
    ret = ret && fp.AddConstant("logten",
                                FUNCTIONPARSERTYPES::fp_const_log10<Value_t>());
    ret = ret && fp.AddConstant("CONST", Value_t(CONST));
    if(!ret)
    {
        std::cout << "Ooops! AddConstant() didn't work" << std::endl;
        return false;
    }

    ret = fp.AddUnit("doubled", 2);
    ret = ret && fp.AddUnit("tripled", 3);
    if(!ret)
    {
        std::cout << "Ooops! AddUnit() didn't work" << std::endl;
        return false;
    }

    ret = fp.AddFunctionWrapper
        ("sub", UserDefFuncWrapper<Value_t>(userDefFuncSub<Value_t>), 2);
    ret = ret && fp.AddFunction("sqr", userDefFuncSqr<Value_t>, 1);
    ret = ret && fp.AddFunction("value", userDefFuncValue<Value_t>, 0);
    if(!ret)
    {
        std::cout << "Ooops! AddFunction(ptr) didn't work" << std::endl;
        return false;
    }

    UserDefFuncWrapper<Value_t>* wrapper =
        dynamic_cast<UserDefFuncWrapper<Value_t>*>
        (fp.GetFunctionWrapper("sub"));
    if(!wrapper || wrapper->counter() != 0)
    {
        std::cout << "Ooops! AddFunctionWrapper() didn't work" << std::endl;
        return false;
    }

    FunctionParserBase<Value_t> SqrFun, SubFun, ValueFun;
    if(verbosityLevel >= 3) std::cout << "Parsing SqrFun... ";
    SqrFun.Parse("x*x", "x");
    if(verbosityLevel >= 3) std::cout << "\nParsing SubFun... ";
    SubFun.Parse("x-y", "x,y");
    if(verbosityLevel >= 3) std::cout << "\nParsing ValueFun... ";
    ValueFun.Parse("5", "");
    if(verbosityLevel >= 3) std::cout << std::endl;

    ret = fp.AddFunction("psqr", SqrFun);
    ret = ret && fp.AddFunction("psub", SubFun);
    ret = ret && fp.AddFunction("pvalue", ValueFun);
    if(!ret)
    {
        std::cout << "Ooops! AddFunction(parser) didn't work" << std::endl;
        return false;
    }

    // Test repeated constant addition
    // -------------------------------
   {using namespace FUNCTIONPARSERTYPES; // For a possible custom < operator
    for(Value_t value = 0; value < Value_t(20); value += 1)
    {
        if(!fp.AddConstant("TestConstant", value))
        {
            std::cout << "Ooops2! AddConstant() didn't work" << std::endl;
            return false;
        }

        fp.Parse("TestConstant", "");
        if(fp.Eval(0) != value)
        {
            if(value == Value_t(0)) std::cout << "Usage of 'TestConstant' failed\n";
            else std::cout << "Changing the value of 'TestConstant' failed\n";
            return false;
        }
    }}

    bool allRegressionTestsOk = true;
    std::ostringstream briefErrorMessages;

    std::string prev_test_prefix;
    const unsigned maxtests = ~0U; // unknown
    /*    sizeof(RegressionTests<Value_t>::Tests)
      / sizeof(RegressionTests<Value_t>::Tests[0]); */
    for(unsigned i = 0; i < maxtests; ++i)
    {
        const TestType<Value_t>& testData = RegressionTests<Value_t>::Tests[i];
        if(!testData.testName) break;

        if(!IsSelectedTest(testData.testName)) continue;

        const int retval =
            fp.Parse(testData.funcString, testData.paramString,
                     testData.useDegrees);
        if(retval >= 0)
        {
            std::cout <<
                "With FunctionParserBase<" << valueType << ">"
                "\nin \"" << testData.funcString <<
                "\" (\"" << testData.paramString <<
                "\"), col " << retval <<
                ":\n" << fp.ErrorMsg() << std::endl;
            return false;
        }

        //fp.PrintByteCode(std::cout);
        if(verbosityLevel >= 3)
        {
            std::cout
                << /*std::right <<*/ std::setw(2)
                << testData.testName << ": \""
                << testData.funcString << "\" ("
                << FUNCTIONPARSERTYPES::fp_pow
                ((testData.paramMax - testData.paramMin) /
                 testData.paramStep,
                 Value_t( (int) testData.paramAmount))
                << " param. combinations): " << std::flush;
        }
        else if(verbosityLevel == 2)
        {
            const char* tn = testData.testName;
            const char* p = std::strrchr(tn, '/');
            if(!p)
                { prev_test_prefix = ""; std::cout << tn; }
            else
            {
                std::string path_prefix(tn, p-tn);
                if(path_prefix == prev_test_prefix)
                    std::cout << (p+1);
                else
                    { if(!prev_test_prefix.empty()) std::cout << std::endl;
                      std::cout << tn;
                      prev_test_prefix = path_prefix; }
            }
            std::cout << std::flush << " ";
        }

        bool thisTestOk =
            runRegressionTest(fp, testData, valueType + ", not optimized",
                              Epsilon<Value_t>(), briefErrorMessages);

        if(thisTestOk)
        {
            if(verbosityLevel >= 3) std::cout << "Ok." << std::endl;

            fp.Optimize();
            //fp.PrintByteCode(std::cout);

            if(verbosityLevel >= 3)
                std::cout << "    Optimized: " << std::flush;

            thisTestOk =
                runRegressionTest(fp, testData,
                                  valueType + ", after optimization",
                                  Epsilon<Value_t>(), briefErrorMessages);
            if(thisTestOk)
            {
                if(verbosityLevel >= 3)
                    std::cout << "(Calling Optimize() several times) "
                              << std::flush;

                for(int j = 0; j < 20; ++j)
                    fp.Optimize();

                /* Sometimes literals drift when the optimizer is run many
                   times, which can become significant with floats. The only
                   purpose to test running the optimizer several times is just
                   to see that it doesn't break. It's not intended to be called
                   several times normally. Hence just skip testing with floats,
                   because the drift just causes differences larger than
                   epsilon...
                */
                if(valueType != "float")
                {
                    thisTestOk =
                        runRegressionTest
                        (fp, testData,
                         valueType + ", after several optimization runs",
                         Epsilon<Value_t>(), briefErrorMessages);
                }

                if(thisTestOk)
                {
                    thisTestOk =
                        testVariableDeduction(fp, testData, briefErrorMessages);

                    if(thisTestOk && verbosityLevel >= 3)
                        std::cout << "Ok." << std::endl;
                }
            }
        } // if(thisTestOk)

        if(!thisTestOk) allRegressionTestsOk = false;

        if(verbosityLevel == 1)
            std::cout << (thisTestOk ? "." : "!") << std::flush;
    } // for(unsigned i = 0; i < maxtests; ++i)

    if(allRegressionTestsOk)
    {
        if(verbosityLevel == 1 || verbosityLevel == 2)
            std::cout << std::endl;
    }
    else if(verbosityLevel <= 1)
    {
        if(verbosityLevel == 1) std::cout << "\n";
        std::cout << briefErrorMessages.str() << std::flush;
    }

    if(verbosityLevel >= 2)
        std::cout << "User-defined function \"sub\" was called "
                  << (dynamic_cast<UserDefFuncWrapper<Value_t>*>
                      (fp.GetFunctionWrapper("sub"))->counter())
                  << " times." << std::endl;

    return allRegressionTestsOk;
}

//=========================================================================
// Optimizer tests
//=========================================================================
namespace OptimizerTests
{
    // --------------------------------------------------------------------
    // Optimizer test 1
    // --------------------------------------------------------------------
    /* Tests functions of the form "A(x^B)^C op D(x^E)^F", where:
       - A,D = {sin,cos,tan,sinh,cosh,tanh,exp}
       - B,E = {1,2}
       - C,F = {-2,-1,0,1,2}
       - op = +, *
    */
    struct MathFuncData
    {
        DefaultValue_t (*mathFunc)(DefaultValue_t d);
        const char* funcName;
    };

    const MathFuncData mathFuncs[] =
    {
        { &std::sin, "sin" }, { &std::cos, "cos" }, { &std::tan, "tan" },
        { &std::sinh, "sinh" }, { &std::cosh, "cosh" }, { &std::tanh, "tanh" },
        { &std::exp, "exp" }
    };
    const unsigned mathFuncsAmount = sizeof(mathFuncs) / sizeof(mathFuncs[0]);

    unsigned mathFuncIndexA, mathFuncIndexD;
    int exponent_B, exponent_E;
    int exponent_C, exponent_F;
    unsigned operatorIndex;

    DefaultValue_t evaluateFunction(const DefaultValue_t* params)
    {
        const DefaultValue_t x = params[0];
        const MathFuncData& data1 = mathFuncs[mathFuncIndexA];
        const MathFuncData& data2 = mathFuncs[mathFuncIndexD];

        const DefaultValue_t angle1 =
            (exponent_B == 1 ? x : std::pow(x, exponent_B));
        const DefaultValue_t angle2 =
            (exponent_E == 1 ? x : std::pow(x, exponent_E));
        const DefaultValue_t part1 =
            std::pow(data1.mathFunc(angle1), exponent_C);
        const DefaultValue_t part2 =
            std::pow(data2.mathFunc(angle2), exponent_F);

        if(operatorIndex == 0) return part1 + part2;
        return part1 * part2;
    }

    bool runCurrentTrigCombinationTest()
    {
        const MathFuncData& data1 = mathFuncs[mathFuncIndexA];
        const MathFuncData& data2 = mathFuncs[mathFuncIndexD];

        std::ostringstream os;
        os << data1.funcName << "(x^" << exponent_B << ")^" << exponent_C;
        if(operatorIndex == 0) os << "+";
        else os << "*";
        os << data2.funcName << "(x^" << exponent_E << ")^" << exponent_F;
        const std::string funcString = os.str();

        const TestType<DefaultValue_t> testData =
        {
            1, -4.0, 4.0, 0.49, false, &evaluateFunction, 0, 0, "x",
            "'trig. combo optimizer test'", funcString.c_str()
        };

        DefaultParser parser;
        if(parser.Parse(funcString, "x") >= 0)
        {
            std::cout << "Oops: Function \"" << funcString
                      << "\" was malformed." << std::endl;
            return false;
        }

        std::ostringstream briefErrorMessages;

        if(!runRegressionTest(parser, testData, "DefaultValue_t",
                              Epsilon<DefaultValue_t>(), briefErrorMessages))
        {
            if(verbosityLevel == 1)
                std::cout << "\n - " << briefErrorMessages.str() << std::flush;
            return false;
        }
        return true;
    }

    bool runTrigCombinationTests()
    {
        unsigned testCounter = 0;

        for(mathFuncIndexA = 0;
            mathFuncIndexA < mathFuncsAmount;
            ++mathFuncIndexA)
        {
            for(mathFuncIndexD = 0;
                mathFuncIndexD < mathFuncsAmount;
                ++mathFuncIndexD)
            {
                for(exponent_B = 1; exponent_B <= 2; ++exponent_B)
                {
                    for(exponent_E = 1; exponent_E <= 2; ++exponent_E)
                    {
                        for(exponent_C = -2; exponent_C <= 2; ++exponent_C)
                        {
                            for(exponent_F = -2; exponent_F <= 2; ++exponent_F)
                            {
                                for(operatorIndex = 0;
                                    operatorIndex < 2;
                                    ++operatorIndex)
                                {
                                    ++testCounter;
                                    if(!runCurrentTrigCombinationTest())
                                        return false;
                                }
                            }
                        }
                    }
                }
            }
        }

        if(verbosityLevel >= 1)
            std::cout << " (" << testCounter << ")" << std::flush;
        return true;
    }


    // --------------------------------------------------------------------
    // Optimizer test 2
    // --------------------------------------------------------------------
    /* Tests functions of the form "A op B [op C]", where
       A, B, C = { var, !var, !!var, var comp value }
       var = A -> x, B -> y, C -> z
       comp = { <, <=, =, !=, >, >= }
       value = { -1, -.5, 0, .5, 1 }
       op = { and, or, not and, not or }
    */
    // opIndex = 0-32 for doubles, 0-20 for ints
    const char* getOperandString(char varName, unsigned opIndex)
    {
        if(opIndex <= 2)
        {
            static char operand[] = "!!x";
            operand[2] = varName;
            return operand + (2 - opIndex);
        }

        opIndex -= 3;
        const unsigned compIndex = opIndex % 6, valueIndex = opIndex / 6;
        assert(valueIndex <= 4);

        static const char* const comp[] =
            { "< ", "<=", "= ", "!=", "> ", ">=" };
        static const char* const value[] =
            { "-1 ", "0  ", "1  ", ".5 ", "-.5" };
        static char expression[] = "(x<=-.5)";

        expression[1] = varName;
        expression[2] = comp[compIndex][0];
        expression[3] = comp[compIndex][1];
        expression[4] = value[valueIndex][0];
        expression[5] = value[valueIndex][1];
        expression[6] = value[valueIndex][2];
        return expression;
    }

    template<typename Value_t>
    Value_t getOperandValue(Value_t varValue, unsigned opIndex)
    {
        using namespace FUNCTIONPARSERTYPES;

        switch(opIndex)
        {
          case 0: return varValue;
          case 1: return fp_not(varValue);
          case 2: return fp_notNot(varValue);
        }

        opIndex -= 3;
        const unsigned compIndex = opIndex % 6, valueIndex = opIndex / 6;

        static const Value_t value[] =
            { -1, 0, 1, Value_t(.5), Value_t(-.5) };

        switch(compIndex)
        {
          case 0: return fp_less(varValue, value[valueIndex]);
          case 1: return fp_lessOrEq(varValue, value[valueIndex]);
          case 2: return fp_equal(varValue, value[valueIndex]);
          case 3: return fp_nequal(varValue, value[valueIndex]);
          case 4: return fp_greater(varValue, value[valueIndex]);
          case 5: return fp_greaterOrEq(varValue, value[valueIndex]);
        }
        assert(false);
        return 0;
    }

    // exprIndex = 0-3
    std::string getBooleanExpression(const std::string& operand1,
                                     const std::string& operand2,
                                     unsigned exprIndex)
    {
        switch(exprIndex)
        {
          case 0: return operand1 + "&" + operand2;
          case 1: return operand1 + "|" + operand2;
          case 2: return "!(" + operand1 + "&" + operand2 + ")";
          case 3: return "!(" + operand1 + "|" + operand2 + ")";
        }
        assert(false);
        return "";
    }

    template<typename Value_t>
    Value_t getBooleanValue(Value_t operand1Value, Value_t operand2Value,
                            unsigned exprIndex)
    {
        using namespace FUNCTIONPARSERTYPES;
        switch(exprIndex)
        {
          case 0: return fp_and(operand1Value, operand2Value);
          case 1: return fp_or(operand1Value, operand2Value);
          case 2: return fp_not(fp_and(operand1Value, operand2Value));
          case 3: return fp_not(fp_or(operand1Value, operand2Value));
        }
        assert(false);
        return 0;
    }

    bool updateIndices(unsigned* operandIndices, unsigned* exprIndices,
                       unsigned operands, const unsigned maxOperandIndex)
    {
        for(unsigned oi = 0; oi < operands; ++oi)
        {
            if(++operandIndices[oi] <= maxOperandIndex)
                return true;
            operandIndices[oi] = 0;
        }
        for(unsigned ei = 0; ei < operands-1; ++ei)
        {
            if(++exprIndices[ei] <= 3) return true;
            exprIndices[ei] = 0;
        }
        return false;
    }

    template<typename Value_t, unsigned varsAmount>
    bool runBooleanComparisonEvaluation(const unsigned* operandIndices,
                                        const unsigned* exprIndices,
                                        const unsigned operands,
                                        FunctionParserBase<Value_t>& fparser,
                                        const std::string& functionString,
                                        bool optimized)
    {
        const bool isIntegral = FUNCTIONPARSERTYPES::IsIntType<Value_t>::result;
        const unsigned varValuesToTest = isIntegral ? 3 : 4;

        static const Value_t values[] =
            { -1, 0, 1, Value_t(0.5), Value_t(-0.5) };
        static unsigned valueIndices[varsAmount];
        static Value_t variableValues[varsAmount];

        for(unsigned i = 0; i < operands; ++i) valueIndices[i] = 0;

        bool stop = false;
        while(!stop)
        {
            for(unsigned i = 0; i < operands; ++i)
                variableValues[i] = values[valueIndices[i]];

            const Value_t parserValue = fparser.Eval(variableValues);

            Value_t correctValue = getOperandValue(variableValues[0],
                                                   operandIndices[0]);

            for(unsigned i = 1; i < operands; ++i)
                correctValue =
                    getBooleanValue(correctValue,
                                    getOperandValue(variableValues[i],
                                                    operandIndices[i]),
                                    exprIndices[i-1]);

            if(FUNCTIONPARSERTYPES::fp_nequal(parserValue, correctValue))
            {
                const bool isIntegral =
                    FUNCTIONPARSERTYPES::IsIntType<Value_t>::result;
                if(verbosityLevel >= 2)
                {
                    using namespace FUNCTIONPARSERTYPES;
                    std::cout
                        << "\nFor function \"" << functionString
                        << "\" (";
                    for(unsigned i = 0; i < operands; ++i)
                        std::cout << (i>0 ? "," : "")
                                  << variableValues[i];
                    std::cout
                        << "): Parser<"
                        << (isIntegral ? "long" : "double")
                        << ">"
                        << (optimized ? " (optimized)" : "")
                        << "\nreturned " << parserValue
                        << " instead of " << correctValue
                        << std::endl;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
                    fparser.PrintByteCode(std::cout);
#endif
                }
                else if(verbosityLevel >= 1)
                {
                    using namespace FUNCTIONPARSERTYPES;
                    std::cout << "<" << (isIntegral ? "long" : "double");
                    std::cout << (optimized ? ",optimized" : "");
                    std::cout << ">\"" << functionString
                              << "\"(";
                    for(unsigned i = 0; i < operands; ++i)
                        std::cout << (i>0 ? "," : "")
                                  << variableValues[i];
                    std::cout << "): ";
                    std::cout << parserValue << " vs " << correctValue;
                    std::cout << "\n";
                }
                return false;
            }

            stop = true;
            for(unsigned i = 0; i < operands; ++i)
            {
                if(++valueIndices[i] < varValuesToTest)
                { stop = false; break; }
                valueIndices[i] = 0;
            }
        }

        return true;
    }

    template<typename Value_t>
    bool runBooleanComparisonTestsForType()
    {
        const bool isIntegral = FUNCTIONPARSERTYPES::IsIntType<Value_t>::result;
        const unsigned maxOperandIndex = isIntegral ? 20 : 32;

        const char varNames[] = { 'x', 'y', 'z' };
        const char* const varString = "x,y,z";
        const unsigned varsAmount = sizeof(varNames) / sizeof(varNames[0]);

        unsigned operandIndices[varsAmount];
        unsigned exprIndices[varsAmount - 1];

        unsigned testCounter = 0;
        FunctionParserBase<Value_t> fparser;

        bool errors = false;

        for(unsigned operands = 2; operands <= varsAmount; ++operands)
        {
            for(unsigned i = 0; i < operands; ++i) operandIndices[i] = 0;
            for(unsigned i = 0; i < operands-1; ++i) exprIndices[i] = 0;

            do
            {
                // Generate function string:
                std::string functionString =
                    getOperandString(varNames[0], operandIndices[0]);

                for(unsigned i = 1; i < operands; ++i)
                    functionString =
                        getBooleanExpression
                        (i == 1 ? functionString : "(" + functionString + ")",
                         getOperandString(varNames[i], operandIndices[i]),
                         exprIndices[i-1]);

                //std::cout << '"' << functionString << "\"\n";

                // Parse function string:
                int errorIndex = fparser.Parse(functionString, varString);
                if(errorIndex >= 0)
                {
                    std::cout << "\nOops! Function \"" << functionString
                              << "\" was malformed.\n";
                    return false;
                }

                // Evaluate function and test for correctness:
                if(!runBooleanComparisonEvaluation<Value_t, varsAmount>
                   (operandIndices, exprIndices, operands,
                    fparser, functionString, false))
                {
                    if (verbosityLevel < 1) return false;
                    errors = true;
                }

                fparser.Optimize();

                if(!runBooleanComparisonEvaluation<Value_t, varsAmount>
                   (operandIndices, exprIndices, operands,
                    fparser, functionString, true))
                {
                    if (verbosityLevel < 1) return false;
                    errors = true;
                }

                ++testCounter;
            }
            while(updateIndices(operandIndices, exprIndices,
                                operands, maxOperandIndex));
        }
        if(errors) return false;

        if(verbosityLevel >= 1)
            std::cout << " (" << testCounter << ")" << std::flush;

        return true;
    }
}

int testOptimizer1()
{
    return OptimizerTests::runTrigCombinationTests();
}

int testOptimizer2()
{
    return OptimizerTests::runBooleanComparisonTestsForType<DefaultValue_t>();
}

int testOptimizer3()
{
#ifdef FP_SUPPORT_LONG_INT_TYPE
    return OptimizerTests::runBooleanComparisonTestsForType<long>();
#else
    return -1;
#endif
}


//=========================================================================
// Help output
//=========================================================================
void printAvailableTests(std::vector<std::string>& tests)
{
    std::cout << "Available tests:\n";
    std::size_t column=0;
    std::string prev_test_prefix;

    bool counting_tests = false;
    long last_count     = 0, count_length = 0;

    for(std::size_t a=0; a<tests.size(); ++a)
    {
        std::string tn = tests[a];
        std::size_t p = tn.rfind('/');
        if(p == tn.npos)
            prev_test_prefix = "";
        else
        {
            {
                std::string path_prefix(tn, 0, p);
                if(path_prefix != prev_test_prefix)
                {
                    if(counting_tests && count_length > 1)
                    {
                        std::ostringstream tmp; tmp << "-" << last_count;
                        std::cout << tmp.str(); column += tmp.str().size();
                    }
                    counting_tests = false;
                    if(column) { std::cout << std::endl; column=0; }
                    prev_test_prefix = path_prefix;
                    std::cout << "    " << path_prefix << "/\n";
                }
            }
            tn.erase(0, p+1);
        }
        if(column+tn.size() >= 76) { column=0; std::cout << "\n"; }
        if(column==0) { std::cout << "        "; column+=8; }
        else { std::cout << " "; column+=1; }

        /* TODO: Rewrite this such that backspaces are not needed,
         *       because they don't work with util-linux's "more"
         */
        char* endptr = 0;
        long val = strtol(tn.c_str(), &endptr, 10);
        if(!*endptr)
        {
            if(!counting_tests)
            {
                counting_tests = true; count_length = 1; last_count = val;
            }
            else if(val == last_count+1)
            {
                ++count_length;
                last_count = val; std::cout << "\b"; --column; continue;
            }
            else if(count_length > 1)
            {
                std::ostringstream tmp; tmp << "\b-" << last_count << " ";
                std::cout << tmp.str(); column += tmp.str().size();
                counting_tests = false;
            }
            else counting_tests = false;
        }
        else if(counting_tests && count_length > 1)
        {
            std::ostringstream tmp; tmp << "\b-" << last_count << " ";
            std::cout << tmp.str(); column += tmp.str().size();
            counting_tests = false;
        }
        else counting_tests = false;

        std::cout << tn;
        column += tn.size();
    }
    if(column) std::cout << std::endl;
}

//=========================================================================
// Main
//=========================================================================
int main(int argc, char* argv[])
{
    const char* const optionsHelpText =
        "    -q                Quiet (no progress, brief error reports)\n"
        "    -v                Verbose (progress, full error reports)\n"
        "    -vv               Very verbose\n"
        "    -tests <tests>    Select tests to perform, wildcards ok (implies -noalgo)\n"
        "                      Example: -tests 'cmp*'\n"
        "    -tests help       List available tests\n"
        "    -d                Test double datatype\n"
        "    -f                Test float datatype\n"
        "    -ld               Test long double datatype\n"
        "    -li               Test long int datatype\n"
        "    -mf, -mpfr        Test MpfrFloat datatype\n"
        "    -gi, -gmpint      Test GmpInt datatype\n"
        "    -cd               Test std::complex<double> datatype\n"
        "    -cf               Test std::complex<float> datatype\n"
        "    -cld              Test std::complex<long double> datatype\n"
        "    -algo <n>         Run only algorithmic test <n>\n"
        "    -noalgo           Skip all algorithmic tests\n"
        "    -skipSlowAlgo     Skip slow algorithmic tests\n"
        "    -h, --help        This help\n";

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    MpfrFloat::setDefaultMantissaBits(80);
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
    GmpInt::setDefaultNumberOfBits(80);
#endif

    bool skipSlowAlgo = false;
    bool runAllTypes = true;
    bool runAlgoTests = true;
    bool run_d = false, run_f = false, run_ld = false;
    bool run_li = false, run_mf = false, run_gi = false;
    bool run_cd = false, run_cf = false, run_cld = false;
    unsigned runAlgoTest = 0;

    for(int i = 1; i < argc; ++i)
    {
        if(std::strcmp(argv[i], "-q") == 0) verbosityLevel = 0;
        else if(std::strcmp(argv[i], "-v") == 0) verbosityLevel = 2;
        else if(std::strcmp(argv[i], "-vv") == 0) verbosityLevel = 3;
        else if(std::strcmp(argv[i], "-vvv") == 0) verbosityLevel = 4;
        else if(std::strcmp(argv[i], "-noalgo") == 0) runAlgoTests = false;
        else if(std::strcmp(argv[i], "-skipSlowAlgo") == 0) skipSlowAlgo = true;
        else if(std::strcmp(argv[i], "-algo") == 0)
        {
            if(i+1 < argc) runAlgoTest = std::atoi(argv[++i]);
            runAlgoTests = true;
        }
        else if(std::strcmp(argv[i], "-tests") == 0)
        {
            runAlgoTests = false;

            std::vector<std::string> tests;
#ifndef FP_DISABLE_DOUBLE_TYPE
            for(unsigned a=0; RegressionTests<double>::Tests[a].testName; ++a)
                tests.push_back(RegressionTests<double>::Tests[a].testName);
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
            for(unsigned a=0; RegressionTests<float>::Tests[a].testName; ++a)
                tests.push_back(RegressionTests<float>::Tests[a].testName);
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
            for(unsigned a=0; RegressionTests<long double>::Tests[a].testName; ++a)
                tests.push_back(RegressionTests<long double>::Tests[a].testName);
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
            for(unsigned a=0; RegressionTests<long>::Tests[a].testName; ++a)
                tests.push_back(RegressionTests<long>::Tests[a].testName);
#endif
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
            for(unsigned a=0; RegressionTests<MpfrFloat>::Tests[a].testName; ++a)
                tests.push_back(RegressionTests<MpfrFloat>::Tests[a].testName);
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
            for(unsigned a=0; RegressionTests<GmpInt>::Tests[a].testName; ++a)
                tests.push_back(RegressionTests<GmpInt>::Tests[a].testName);
#endif
            std::sort(tests.begin(), tests.end(), natcomp);
            tests.erase(std::unique(tests.begin(), tests.end()), tests.end());

            if(std::strcmp(argv[i+1], "help") == 0)
            {
                printAvailableTests(tests);
                return 0;
            }
            while(i+1 < argc && argv[i+1][0] != '-')
            {
                const char* t = argv[++i];
                bool ok = false;
                for(std::size_t a=0; a<tests.size(); ++a)
                    if(WildMatch_Dirmask(t, tests[a].c_str()))
                    { ok=true; break; }
                if(!ok)
                {
                    std::cout << "No such test: " << t
                              << "\n\"testbed -tests help\" to list "
                              << "available tests.\n";
                    return -1;
                }
                selectedRegressionTests.push_back(t);
            }
        }
        else if(std::strcmp(argv[i], "-d") == 0
             || std::strcmp(argv[i], "-double") == 0)
            runAllTypes = false, run_d = true;
        else if(std::strcmp(argv[i], "-f") == 0
             || std::strcmp(argv[i], "-float") == 0)
            runAllTypes = false, run_f = true;
        else if(std::strcmp(argv[i], "-ld") == 0
             || std::strcmp(argv[i], "-longdouble") == 0)
            runAllTypes = false, run_ld = true;
        else if(std::strcmp(argv[i], "-li") == 0
             || std::strcmp(argv[i], "-longint") == 0)
            runAllTypes = false, run_li = true;
        else if(std::strcmp(argv[i], "-mf") == 0
             || std::strcmp(argv[i], "-mpfr") == 0)
            runAllTypes = false, run_mf = true;
        else if(std::strcmp(argv[i], "-gi") == 0
             || std::strcmp(argv[i], "-gmpint") == 0)
            runAllTypes = false, run_gi = true;
        else if(std::strcmp(argv[i], "-cd") == 0)
            runAllTypes = false, run_cd = true;
        else if(std::strcmp(argv[i], "-cf") == 0)
            runAllTypes = false, run_cf = true;
        else if(std::strcmp(argv[i], "-cld") == 0)
            runAllTypes = false, run_cld = true;

        else if(std::strcmp(argv[i], "--help") == 0
             || std::strcmp(argv[i], "-help") == 0
             || std::strcmp(argv[i], "-h") == 0
             || std::strcmp(argv[i], "/?") == 0)
        {
            std::cout <<
                "FunctionParser testbed " << kVersionNumber <<
                "\n\nUsage: " << argv[0] << " [<option> ...]\n"
                "\n" << optionsHelpText;
            return 0;
        }
        else if(std::strlen(argv[i]) > 0)
        {
            std::cout << "Unknown option: '" << argv[i] << "'\n";
            return 1;
        }
    }

    if(selectedRegressionTests.empty())
        selectedRegressionTests.push_back("*");

    DefaultParser fp0;

    // Test that the parser doesn't crash if Eval() is called before Parse():
    fp0.Eval(0);

    const char* const delimiterTestFunction = "x+y } ";
    fp0.setDelimiterChar('}');
    int res = fp0.Parse(delimiterTestFunction, "x,y");
    if(fp0.GetParseErrorType() != fp0.FP_NO_ERROR || res != 4)
    {
        std::cout << "Delimiter test \"" << delimiterTestFunction
                  << "\" failed at " << res << ": " << fp0.ErrorMsg()
                  << std::endl;
        return 1;
    }
    fp0.Parse("x+}y", "x,y");
    if(fp0.GetParseErrorType() == fp0.FP_NO_ERROR)
    {
        std::cout << "Erroneous function with delimiter didn't fail"
                  << std::endl;
        return 1;
    }

    bool allTestsOk = true;

    if(!runAllTypes || runAlgoTest == 0)
    {
#ifndef FP_DISABLE_DOUBLE_TYPE
        if(runAllTypes || run_d)
            if(!runRegressionTests<double>("double"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
        if(runAllTypes || run_f)
            if(!runRegressionTests<float>("float"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
        if(runAllTypes || run_ld)
            if(!runRegressionTests<long double>("long double"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
        if(runAllTypes || run_li)
            if(!runRegressionTests<long>("long int"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
        if(runAllTypes || run_mf)
            if(!runRegressionTests<MpfrFloat>("MpfrFloat"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
        if(runAllTypes || run_gi)
            if(!runRegressionTests<GmpInt>("GmpInt"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
        if(runAllTypes || run_cd)
            if(!runRegressionTests<std::complex<double> >("std::complex<double>"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
        if(runAllTypes || run_cf)
            if(!runRegressionTests<std::complex<float> >("std::complex<float>"))
                allTestsOk = false;
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
        if(runAllTypes || run_cld)
            if(!runRegressionTests<std::complex<long double> >("std::complex<long double>"))
                allTestsOk = false;
#endif
    }

////////////////////////////
////////////////////////////
////////////////////////////
////////////////////////////

    // Misc. tests
    // -----------
    const struct
    {
        const char* const testName;
        int(*testFunction)();
    }
    algorithmicTests[] =
    {
        { "Copy constructor and assignment", &TestCopying },
        { "Error situations", &TestErrorSituations },
        { "Whitespaces", &WhiteSpaceTest },
        { "Optimizer test 1 (trig. combinations)", &testOptimizer1 },
        { "Optimizer test 2 (bool combinations, double)",
          (skipSlowAlgo || (!runAllTypes && !run_d)) ? 0 : &testOptimizer2 },
        { "Optimizer test 3 (bool combinations, long)",
          (!runAllTypes && !run_li) ? 0 : &testOptimizer3 },
        { "Integral powers",  &TestIntPow },
        { "UTF8 test", skipSlowAlgo ? 0 : &UTF8Test },
        { "Identifier test", &TestIdentifiers },
        { "Used-defined functions", &testUserDefinedFunctions },
        { "Multithreading", &testMultithreadedEvaluation }
    };

    const unsigned algorithmicTestsAmount =
        sizeof(algorithmicTests) / sizeof(algorithmicTests[0]);

    if(runAlgoTests)
    {
        for(unsigned i = 0; i < algorithmicTestsAmount; ++i)
        {
            if(runAlgoTest >= 1 && runAlgoTest <= algorithmicTestsAmount &&
               runAlgoTest != i+1)
                continue;

            if(verbosityLevel >= 1)
                std::cout << "Algo test " << i+1 << ": "
                          << algorithmicTests[i].testName << std::flush;

            if(!algorithmicTests[i].testFunction)
            {
                if(verbosityLevel >= 1)
                    std::cout << ": Skipped." << std::endl;
                continue;
            }

            int result = algorithmicTests[i].testFunction();

            if(result == 0)
            {
                allTestsOk = false;
                if(verbosityLevel == 0)
                    std::cout << "Algo test " << i+1 << ": "
                              << algorithmicTests[i].testName;
                if(verbosityLevel <= 1)
                    std::cout << ": FAILED." << std::endl;
            }
            else if(verbosityLevel >= 1)
            {
                if(result < 0 )
                    std::cout << ": (No support)" << std::endl;
                else
                    std::cout << ": Ok." << std::endl;
            }
        }
    }

    if(!allTestsOk)
    {
        std::cout << "Some tests failed." << std::endl;
        return 1;
    }
    if(verbosityLevel == 1)
        std::cout << "================= All tests OK =================\n";
    else if(verbosityLevel >= 2)
        std::cout << "================================================\n"
                  << "================= All tests OK =================\n"
                  << "================================================\n";
    return 0;
}
