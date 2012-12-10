#include "fparser.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>

struct Counts { unsigned opcodes, muls; };

namespace
{
    static int Log2up(unsigned long value)
    {
        unsigned long res = 0;
        while(value > (1LU << res)) ++res;
        return res;
    }
    static int PopCount(unsigned long value) // population count
    {
        static const int counts[16] = { 0,1,1,2, 1,2,2,3, 1,2,2,3, 2,3,3,4 };
        int result = 0;
        while(value > 0) { result += counts[value & 15]; value >>= 4; }
        return result;
    }
}

Counts generateOpcodesForExp(unsigned n, bool print, double& value)
{
    Counts retval = { 0, 0 };
    if(n > 1)
    {
        if(n % 2 == 1)
        {
            if(true) // USE_DIVISIONS_AS_WELL
            {
                int next_power_of_2 = Log2up(n);
                int popcount        = PopCount(n);
                //fprintf(stderr, "Considering %u... %u/%u\n", n, popcount,next_power_of_2);
                //fflush(stderr);
                if(n != 3 && popcount >= next_power_of_2 )
                {
                    int morepower = 1 << next_power_of_2;
                    int lesspower = -(n - morepower);
                    const double valueCopy = value;
                    if(print) std::cout << "dup ";
                    retval = generateOpcodesForExp(morepower, print, value);
                    if(print) std::cout << "fetch ";
                    const double high = value;
                    value = valueCopy;
                    Counts r2 = generateOpcodesForExp(lesspower, print, value);
                    if(print) std::cout << "rdiv ";
                    retval.opcodes += r2.opcodes;
                    retval.muls    += r2.muls;
                    retval.muls += 1;
                    retval.opcodes += 3;
                    value = high / value;
                    return retval;
                }
            }
            if(print) std::cout << "dup ";
            const double valueCopy = value;
            retval = generateOpcodesForExp(n-1, print, value);
            retval.opcodes += 2;
            ++retval.muls;
            if(print) std::cout << "mul ";
            value *= valueCopy;
        }
        else
        {
            retval = generateOpcodesForExp(n/2, print, value);
            ++retval.opcodes;
            ++retval.muls;
            if(print) std::cout << "sqr ";
            value *= value;
        }
    }
    return retval;
}

Counts getParserOpcodesAmount(const std::string& func, double& value)
{
    FunctionParser fp;
    std::string line;

    fp.Parse(func, "x");
    fp.Optimize();
    std::ostringstream buf;
    fp.PrintByteCode(buf);
    std::istringstream lines(buf.str());

    Counts counts = { 0, 0 };
    while(std::getline(lines, line).good())
    {
        ++counts.opcodes;
        const std::string end = line.substr(line.size()-3);
        if(end == "mul"
        || end == "div"
        || end == "rdi"
        || end == "sqr") ++counts.muls;
    }
    --counts.opcodes;

    value = fp.Eval(&value);

    return counts;
}

bool compare(double v1, double v2)
{
    const double Epsilon = .000001;
    const double scale = pow(10.0, floor(log10(fabs(v1))));
    double sv1 = fabs(v1) < Epsilon ? 0 : v1/scale;
    double sv2 = fabs(v2) < Epsilon ? 0 : v2/scale;
    double diff = sv2-sv1;
    return std::fabs(diff) < Epsilon;
}

int main()
{
    std::printf
        ("Number of opcodes generated:\n"
         "Func     Naive     Bisq   Func      Naive     Bisq   Func      Naive     Bisq   Func      Naive     Bisq\n"
         "----     -----     ----   ----      -----     ----   ----      -----     ----   ----      -----     ----\n");

    const Counts minimum = { 0, 0 };

    for(unsigned i = 0; i < 100; ++i)
    {
        for(unsigned col = 0; col < 4; ++col)
        {
            const unsigned exponent = i + 100*col;

            const double value = 1.02;
            double result = exponent == 0 ? 1 : value;
            for(unsigned i = 2; i <= exponent; ++i)
                result *= value;

            std::ostringstream funcStream;
            if(exponent < 10) funcStream << " ";
            funcStream << "x^" << exponent;
            const std::string func = funcStream.str();

            double naiveValue = exponent == 0 ? 1 : value;
            Counts naiveOpcodes = exponent < 2 ? minimum :
                generateOpcodesForExp(exponent, false, naiveValue);
            ++naiveOpcodes.opcodes;

            double bisqValue = value;
            const Counts bisqOpcodes = getParserOpcodesAmount(func, bisqValue);

            if(!compare(naiveValue, result))
            {
                std::cerr << "\nFor exponent " << exponent
                          << " naive algorithm returned \n"
                          << std::setprecision(18)
                          << naiveValue << " instead of " << result << "\n";
                return 1;
            }
            if(!compare(bisqValue, result))
            {
                std::cerr << "\nFor exponent " << exponent
                          << " Bisq algorithm returned \n"
                          << bisqValue << " instead of " << result << "\n";
                return 1;
            }

            //fflush(stderr);
            std::printf("%s: %3u (%2u) %3u (%2u)   ", func.c_str(),
                        naiveOpcodes.opcodes, naiveOpcodes.muls,
                        bisqOpcodes.opcodes, bisqOpcodes.muls);
            //fflush(stdout);
        }
        std::printf("\n");
    }
    return 0;

    /*
    for(unsigned i = 2; i < 20; ++i)
    {
        std::cout << "x^" << i << ": ";
        unsigned amount = generateOpcodesForExp(i, true).opcodes;
        std::cout << ": " << amount << "\n";
    }
    return 0;




    std::string function;
    FunctionParser fparser;

    while(true)
    {
        std::cout << "f(x) = ";
        std::getline(std::cin, function);
        int res = fparser.Parse(function, "x");

        if(res >= 0)
            std::cout << std::string(res+7, ' ') << "^\n"
                      << fparser.ErrorMsg() << "\n\n";
        else
        {
            std::cout << "------- Normal: -------\n";
            fparser.PrintByteCode(std::cout);
            std::cout << "------- Optimized: -------\n";
            fparser.Optimize();
            fparser.PrintByteCode(std::cout);
        }
    }
    */
}
