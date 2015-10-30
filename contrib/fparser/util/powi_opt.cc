#include <iostream>

#define FP_GENERATING_POWI_TABLE
#include "../fpoptimizer/bytecodesynth.cc"
#include "../fpoptimizer/codetree.hh"
#include "extrasrc/fptypes.hh"

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_ByteCode;

#include <iomanip>
namespace
{
    inline void printHex(std::ostream& dest, unsigned n)
    {
        dest.width(8); dest.fill('0'); std::hex(dest); //uppercase(dest);
        dest << n;
    }
}

static
void PrintByteCode(const std::vector<unsigned>& ByteCode,
                   const std::vector<double>& Immed,
                   std::ostream& dest)
{
    for(unsigned IP = 0, DP = 0; IP < ByteCode.size(); ++IP)
    {
        printHex(dest, IP);
        dest << ": ";

        unsigned opcode = ByteCode[IP];

        switch(opcode)
        {
          case cIf:
              dest << "jz\t";
              printHex(dest, ByteCode[IP+1]+1);
              dest << std::endl;
              IP += 2;
              break;

          case cJump:
              dest << "jump\t";
              printHex(dest, ByteCode[IP+1]+1);
              dest << std::endl;
              IP += 2;
              break;
          case cImmed:
              dest.precision(10);
              dest << "push\t" << Immed[DP++] << std::endl;
              break;

          default:
              if(OPCODE(opcode) < VarBegin)
              {
                  std::string n;
                  unsigned params = 1;
                  switch(opcode)
                  {
                    case cNeg: n = "neg"; break;
                    case cAdd: n = "add"; break;
                    case cSub: n = "sub"; break;
                    case cMul: n = "mul"; break;
                    case cDiv: n = "div"; break;
                    case cMod: n = "mod"; break;
                    case cPow: n = "pow"; break;
                    case cEqual: n = "eq"; break;
                    case cNEqual: n = "neq"; break;
                    case cLess: n = "lt"; break;
                    case cLessOrEq: n = "le"; break;
                    case cGreater: n = "gt"; break;
                    case cGreaterOrEq: n = "ge"; break;
                    case cAnd: n = "and"; break;
                    case cOr: n = "or"; break;
                    case cNot: n = "not"; break;
                    case cDeg: n = "deg"; break;
                    case cRad: n = "rad"; break;

#ifdef FP_SUPPORT_OPTIMIZER
                    case cDup: n = "dup"; break;
                    case cInv: n = "inv"; break;
                    case cSqr: n = "sqr"; break;
                    case cFetch:
                        dest << "cFetch(" << ByteCode[++IP] << ")";
                        break;
                    case cPopNMov:
                    {
                        size_t a = ByteCode[++IP];
                        size_t b = ByteCode[++IP];
                        dest << "cPopNMov(" << a << ", " << b << ")";
                        break;
                    }
#endif

                    default:
                        n      = Functions[opcode-cAbs].name;
                        params = Functions[opcode-cAbs].params;
                  }
                  dest << n;
                  if(params != 1) dest << " (" << params << ")";
                  dest << std::endl;
              }
              else
              {
                  dest << "push\tVar" << opcode-VarBegin << std::endl;
              }
        }
    }
}

static long min(long a,long b) { return a<b?a:b;}

int main()
{
    for(long exponent = 2; exponent < 256; ++exponent)
    {
        CodeTree ct;

        double bestres = 0;
        long bestp = -1;

        /* x^40 / x^5 (rdiv) cannot be used when x=0 */
        for(long p=1; p<256; ++p)
        {
            long factor = p&127;
            if(factor & 64)
            {
                factor = -(factor&63) - 1;
                continue;
            }
            if(p & 128)
            {
                if(factor == 0) continue;
                if(factor==1 || exponent % factor != 0) continue;
                if(factor >= exponent) continue;
            }
            else
            {
                if(factor >= exponent) continue;
                if(factor < 0
                && ((  (-factor&(-factor-1)))
                || -factor <= exponent)
                  ) continue;
            }

            //if(p != powi_table[exponent]) continue;

            powi_table[exponent] = p;

            fprintf(stderr, "For %ld, trying %ld (%ld%s)... ",
                exponent,
                p, factor,
                (p&128) ? ", factor": "");

            ByteCodeSynth synth;
            synth.PushVar(VarBegin);
            AssembleSequence(exponent, MulSequence, synth);

            std::vector<unsigned> byteCode;
            std::vector<double>   immed;
            size_t stacktop_max=1;
            synth.Pull(byteCode, immed, stacktop_max);

            double res = 0;
            for(size_t a=0; a<byteCode.size(); ++a)
            {
                if(byteCode[a] == cMul)
                    res += 7;
                else if(byteCode[a] == cDiv || byteCode[a] == cRDiv)
                    res += 11;
                else if(byteCode[a] == cSqr)
                    res += 6.5;
                else if(byteCode[a] == cDup)
                    res += 1;
                else if(byteCode[a] == cPopNMov)
                    { res += 5; a += 2; }
                else if(byteCode[a] == cFetch)
                    { res += 3.5; a += 1; }
            }

            res += stacktop_max*0.3;

            fprintf(stderr, "gets %g, stackmax %u", res, (unsigned)stacktop_max);
            if(res < bestres
            || bestp == -1
              )
            {
                fprintf(stderr, " -- beats %g (%ld)\n", bestres, bestp);
                bestp   = p;
                bestres = res;
            }
            else
                fprintf(stderr, "\n");
            fflush(stderr);

            PrintByteCode(byteCode, immed, std::cerr);
            std::cerr << std::endl;
        }
        powi_table[exponent] = bestp;
    }
    for(unsigned n=0; n<256; ++n)
    {
        if(n%8 == 0) printf("   ");
        printf("%4d,", powi_table[n]);
        if(n%8 == 7) printf(" /*%4d -%4d */\n", n&~7, (n&~7)|7);
    }
}
