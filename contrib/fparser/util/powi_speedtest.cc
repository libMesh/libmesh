#include "fparser.hh"
#include <cstdio>
#include <sstream>
#include <string>
#include <ctime>
#include <cmath>

std::string getFunction(int exponent)
{
    std::ostringstream os;
    os << "x^" << exponent;
    std::string func = os.str();
    while(func.length() < 6) func = " " + func;
    return func;
}

unsigned getEvalsPerMS(FunctionParser& fp)
{
    const unsigned loops = 3000000;

    std::clock_t iTime = std::clock();
    const double value = 1.02;

    for(unsigned i = 0; i < loops; ++i)
        fp.Eval(&value);

    std::clock_t t = std::clock() - iTime;
    return unsigned(std::floor(double(loops)*CLOCKS_PER_SEC/1000.0/t + .5));
}

int main()
{
    FunctionParser fp;

    std::printf
        ("Evaluations per microsecond:\n"
         "  Func Normal  Optim    Func Normal  Optim    Func Normal  Optim    Func Normal  Optim\n"
         "  ---- ------  -----    ---- ------  -----    ---- ------  -----    ---- ------  -----\n");

    for(int row = 0; row < 100; ++row)
    {
        for(int sign = 1; sign >= -1; sign -= 2)
        {
            for(int column = 0; column < 4; ++column)
            {
                const int exponent = (row + column*100) * sign;
                {
                    const std::string func = getFunction(exponent);
                    fp.Parse(func, "x");

                    const unsigned epms1 = (getEvalsPerMS(fp)+50)/100;
                    fp.Optimize();
                    const unsigned epms2 = (getEvalsPerMS(fp)+50)/100;

                    std::printf("%s %4u.%1u %4u.%1u  ", func.c_str(),
                                epms1/10, epms1%10, epms2/10, epms2%10);
                    std::fflush(stdout);
                }
            }
            std::printf("\n");
        }
    }
}
