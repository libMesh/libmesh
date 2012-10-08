#include "fpoptimizer/grammar.hh"
#include <algorithm>
#include <sstream>
#include <cmath>
#include <set>

using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

namespace
{
    class TestFunction
    {
    public:
        std::string fparser_test;
        std::string cpp_test;
    public:
        TestFunction()
            : fparser_test(),
              cpp_test()
        {
        }

        TestFunction(const std::string& s)
            : fparser_test(s),
              cpp_test(s)
        {
        }

        TestFunction(const std::string& fp, const std::string& cp)
            : fparser_test(fp),
              cpp_test(cp)
        {
        }

        struct TestInfo
        {
            const char* fp_pref, *cp_pref;
            const char* fp_sep,  *cp_sep;
            const char* fp_suff, *cp_suff;

            explicit TestInfo(OPCODE opcode)
            {
                fp_pref = cp_pref = "((";
                fp_sep = cp_sep = ",";
                fp_suff = cp_suff = "))";
                switch(opcode)
                {
                    case cMul: fp_sep = cp_sep = ")*("; return;
                    case cAdd: fp_sep = cp_sep = ")+("; return;
                    case cSub: fp_sep = cp_sep = ")-("; return;
                    case cDiv: fp_sep = cp_sep = ")/("; return;
                    case cAnd: fp_sep = cp_sep = ")&("; return;
                    case cOr:  fp_sep = cp_sep = ")|("; return;
                    case cEqual:  fp_sep = ")=("; cp_sep = ")==("; return;
                    case cNEqual:  fp_sep = cp_sep = ")!=("; return;
                    case cLess:  fp_sep = cp_sep = ")<("; return;
                    case cLessOrEq:  fp_sep = cp_sep = ")<=("; return;
                    case cGreater:  fp_sep = cp_sep = ")>("; return;
                    case cGreaterOrEq:  fp_sep = cp_sep = ")>=("; return;
                    case cMod: fp_sep = ")%("; cp_pref = "(fmod("; cp_sep = ","; return;
                    case cPow: fp_sep = ")^("; cp_pref =  "(pow("; cp_sep = ","; return;
                    case cNeg: fp_pref = cp_pref =  "(-("; return;
                    case cNot: fp_pref = cp_pref =  "(!("; return;
                    case cNotNot: fp_pref = cp_pref =  "(!!("; return;
                    case cInv: fp_pref = cp_pref = "(1/("; return;
                    case cCot: fp_pref = "cot"; cp_pref = "(1/tan("; return;
                    case cCsc: fp_pref = "csc"; cp_pref = "(1/sin("; return;
                    case cSec: fp_pref = "sec"; cp_pref = "(1/cos("; return;
                    case cInt: fp_pref = "int"; cp_pref = "(floor(0.5+("; return;
                    #define op(opcode, fpcode, cpcode) \
                        case opcode: \
                            fp_pref = #fpcode "("; \
                            cp_pref = #cpcode "("; \
                            fp_suff = cp_suff = ")"; return
                    op(cAbs,abs,abs);
                    op(cAcos,acos,acos);
                    op(cAcosh,acosh,fp_acosh);
                    op(cAsin,asin,asin);
                    op(cAsinh,asinh,fp_asinh);
                    op(cAtan,atan,atan);
                    op(cAtan2,atan2,atan2);
                    op(cAtanh,atanh,fp_atanh);
                    op(cCeil,ceil,ceil);
                    op(cCos,cos,cos);
                    op(cCosh,cosh,cosh);
                    case cIf:
                        fp_pref = "if((";
                        return;
                    op(cExp,exp,exp);
                    op(cExp2,exp2,exp2);
                    op(cFloor,floor,floor);
                    op(cLog,log,log);
                    op(cLog10,log10,log10);
                    op(cLog2,log2,log2);
                    op(cMax,max,max);
                    op(cMin,min,min);
                    op(cSin,sin,sin);
                    op(cSinh,sinh,sinh);
                    op(cSqrt,sqrt,sqrt);
                    op(cTan,tan,tan);
                    op(cTanh,tanh,tanh);
                    op(cTrunc,trunc,trunc);
                    #undef op
                    case cImmed: return; // does not occur
                    case cJump: return; // does not occur
                    case cDup: return; // does not occur
                    case cFetch: return; // does not occur
                    case cPopNMov: return; // does not occur
                    case cSqr: return; // does not occur
                    case cRDiv: return; // does not occur
                    case cRSub: return; // does not occur
                    case cRSqrt: return; // does not occur (?)
                    case cNop: return; // does not occur
                    case VarBegin: return; // does not occur
                }
            }
        };
    };

    class TestGenerator
    {
        std::map<unsigned, std::string> holders;
    public:
        TestFunction CreateTest(const ParamSpec_ParamHolder& holder)
        {
            std::map<unsigned, std::string>::const_iterator i = holders.find(holder.index);
            if(i != holders.end()) return TestFunction( i->second );
            std::string result;
            if(holder.constraints & Constness_Const)
            {
            get_const_instead:;
                double value = 4.91;
                switch(ImmedConstraint_Value( holder.constraints & ValueMask) )
                {
                    case Value_IsInteger: case Value_EvenInt:
                        value = std::floor(value);
                        break;
                    case Value_OddInt:
                        value = 3.0;
                        break;
                    case Value_Logical:
                        value = 1.0;
                    default: break;
                }
                if(holder.constraints & Oneness_One)
                    value = 1.0;
                if(holder.constraints & Sign_Negative)
                    value = -value;
                std::ostringstream r;
                r.precision(15);
                r << value;
                result = r.str();
            }
            else
            {
                if(holder.constraints & OnenessMask) goto get_const_instead;
                static const char* const predefined[] =
                { "x", "y", "z", "sinh(x)" };
                result = predefined[(holder.index+4-2) % 4];
                // any expression evaluating to the given constraints
                switch(ImmedConstraint_Value( holder.constraints & ValueMask) )
                {
                    case Value_IsInteger:
                        result = "floor(" + result + ")";
                        break;
                    case Value_OddInt:
                    case Value_EvenInt:
                    case Value_NonInteger:
                        goto get_const_instead;
                    case Value_Logical:
                        result = "(" + result + " < 3)";
                        break;
                }
                switch(ImmedConstraint_Sign( holder.constraints & SignMask) )
                {
                    case Sign_Positive:
                        result = "abs(" + result + ")";
                        break;
                    case Sign_Negative:
                        result = "(-cosh(" + result + "))";
                        break;
                }
            }
            return TestFunction( holders[holder.index] = result );
        }

        TestFunction CreateTest(const ParamSpec_NumConstant& n)
        {
            std::ostringstream s;
            s.precision(15);
            s << n.constvalue;
            return TestFunction( s.str() );
        }

        TestFunction CreateTest(const ParamSpec_SubFunctionData& tree, unsigned constraints)
        {
            TestFunction::TestInfo info ( tree.subfunc_opcode );
            std::vector<ParamSpec> params;
            for(unsigned p=0; p<tree.param_count; ++p)
                params.push_back( ParamSpec_Extract(tree.param_list, p) );

            std::vector<TestFunction> result_params;
            for(unsigned p=0; p<tree.param_count; ++p)
            {
                const ParamSpec& param = params[p];
                switch(param.first)
                {
                    case NumConstant:
                    {
                        const ParamSpec_NumConstant& n = *(const ParamSpec_NumConstant*)param.second;
                        result_params.push_back( CreateTest( n ) );
                        break;
                    }
                    case ParamHolder:
                    {
                        const ParamSpec_ParamHolder& n = *(const ParamSpec_ParamHolder*)param.second;
                        result_params.push_back( CreateTest( n ) );
                        break;
                    }
                    case SubFunction:
                    {
                        const ParamSpec_SubFunction& n = *(const ParamSpec_SubFunction*)param.second;
                        result_params.push_back( CreateTest( n.data, n.constraints ) );
                        break;
                    }
                }
            }
            if(tree.restholder_index != 0)
            {
                // add two random params
                if(constraints & Sign_Positive)
                    result_params.push_back( TestFunction( "abs(q)" ) );
                else
                {
                    result_params.push_back( TestFunction( "q" ) );
                    result_params.push_back( TestFunction( "w" ) );
                }
            }

            if(tree.match_type != PositionalParams)
                std::random_shuffle( result_params.begin(), result_params.end() );
            TestFunction result;
            result.fparser_test += info.fp_pref;
            result.cpp_test     += info.cp_pref;
            for(unsigned p=0; p<result_params.size(); ++p)
            {
                if(tree.subfunc_opcode == cIf)
                {
                    static const char* const fp_list[3] = { "", "),(", "),(" };
                    static const char* const cp_list[3] = { "", ")?(", "):(" };
                    info.fp_sep = fp_list[p];
                    info.cp_sep = cp_list[p];
                }
                if(p > 0)
                {
                    result.fparser_test += info.fp_sep;
                    result.cpp_test     += info.cp_sep;
                }
                result.fparser_test += result_params[p].fparser_test;
                result.cpp_test     += result_params[p].cpp_test;
            }
            result.fparser_test += info.fp_suff;
            result.cpp_test     += info.cp_suff;
            return result;
        }
    };
}

static std::set<const Rule*> Rules;

static void FindRules(const Grammar& g)
{
    for(unsigned a=0; a<g.rule_count; ++a)
        Rules.insert(&grammar_rules[g.rule_list[a]]);
}

int main()
{
    FindRules(grammar_optimize_round1);
    FindRules(grammar_optimize_round2);
    FindRules(grammar_optimize_round3);

    for(std::set<const Rule*>::const_iterator
        i = Rules.begin();
        i != Rules.end();
        ++i)
    {
        const Rule& r = **i;

        if(r.logical_context)
        {
            // Skipping logical context rule
            // FIXME: Instead, devise a test that utilizes logical context
            continue;
        }

        ParamSpec_SubFunctionData in_func = r.match_tree;
        if(r.ruletype == ReplaceParams
        && (in_func.subfunc_opcode == cAdd
        || in_func.subfunc_opcode == cMul
        || in_func.subfunc_opcode == cAnd
        || in_func.subfunc_opcode == cOr))
        {
            in_func.restholder_index = 7;
        }
        TestGenerator gen;
        TestFunction test = gen.CreateTest(in_func, 0);

        TestFunction repl;
        if(r.ruletype == ReplaceParams)
        {
            ParamSpec_SubFunctionData repl_func =
            { r.repl_param_count, r.repl_param_list,
              r.match_tree.subfunc_opcode,
              PositionalParams, in_func.restholder_index };
            repl = gen.CreateTest(repl_func, 0);
        }
        else
        {
            ParamSpec p = ParamSpec_Extract(r.repl_param_list, 0);
            if(p.first == SubFunction)
                repl = gen.CreateTest(*(const ParamSpec_SubFunctionData*)p.second, 0);
            else if(p.first == ParamHolder)
                repl = gen.CreateTest(*(const ParamSpec_ParamHolder*)p.second);
            else
                repl = gen.CreateTest(*(const ParamSpec_NumConstant*)p.second);
        }
        ParamSpec_SubFunction tmp = {r.match_tree,0,0};

        std::cout << "echo '---------NEW TEST-----------'\n";
        std::cout << "echo 'Rule: ";
        FPoptimizer_Grammar::DumpParam( ParamSpec(SubFunction, (const void*) &tmp) );
        std::cout << "'\n";
        if(r.ruletype == ProduceNewTree)
        {
            std::cout << "echo '  ->  ";
            FPoptimizer_Grammar::DumpParam(
                                        ParamSpec_Extract(r.repl_param_list, 0) );
            std::cout << "'\n";
        }
        else
        {
            std::cout << "echo '  :   ";
            FPoptimizer_Grammar::DumpParams(
                                        r.repl_param_list, r.repl_param_count );
            std::cout << "'\n";
        }

        std::cout << "./functioninfo '" << test.fparser_test
                  << "' '" << repl.fparser_test << "'\n";
        /*
        std::cout << test.fparser_test <<
             "\n" << test.cpp_test <<
             "\n\n";
        */
    }
}
