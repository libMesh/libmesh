#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <algorithm>
#include <assert.h>

#include "optimize.hh"

using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace
{
    /* Synthesize the given grammatic parameter into the codetree */
    template<typename Value_t>
    CodeTree<Value_t> SynthesizeParam(
        const ParamSpec& parampair,
        MatchInfo<Value_t> & info,
        bool inner = true)
    {
        switch( parampair.first )
        {
            case NumConstant:
              { const ParamSpec_NumConstant<Value_t>& param = *(const ParamSpec_NumConstant<Value_t>*) parampair.second;
                return CodeTreeImmed( param.constvalue );
              }
            case ParamHolder:
              { const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                return info.GetParamHolderValue( param.index );
              }
            case SubFunction:
              { const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                CodeTree<Value_t> tree;
                tree.SetOpcode( param.data.subfunc_opcode );
                for(unsigned a=0; a < param.data.param_count; ++a)
                {
                    CodeTree<Value_t> nparam =
                        SynthesizeParam( ParamSpec_Extract<Value_t>(param.data.param_list, a),
                                         info, true );
                    tree.AddParamMove(nparam);
                }
                if(param.data.restholder_index != 0)
                {
                    std::vector<CodeTree<Value_t> > trees
                        ( info.GetRestHolderValues( param.data.restholder_index ) );
                    tree.AddParamsMove(trees);
                    // ^note: this fails if the same restholder is synth'd twice
                    if(tree.GetParamCount() == 1)
                    {
                        /* Convert cMul <1> into <1> when <1> only contains one operand.
                         * This is redundant code; it is also done in ConstantFolding(),
                         * but it might be better for performance to do it here, too.
                         */
                        assert(tree.GetOpcode() == cAdd || tree.GetOpcode() == cMul
                            || tree.GetOpcode() == cMin || tree.GetOpcode() == cMax
                            || tree.GetOpcode() == cAnd || tree.GetOpcode() == cOr
                            || tree.GetOpcode() == cAbsAnd || tree.GetOpcode() == cAbsOr);
                        tree.Become(tree.GetParam(0));
                    }
                    else if(tree.GetParamCount() == 0)
                    {
                        switch(tree.GetOpcode())
                        {
                            case cAdd: case cOr:
                                tree = CodeTreeImmed(Value_t(0));
                                break;
                            case cMul: case cAnd:
                                tree = CodeTreeImmed(Value_t(1));
                            default: break;
                        }
                    }
                }
                if(inner) tree.Rehash();
                return tree;
              }
        }
        return CodeTree<Value_t> ();
    }
}

namespace FPoptimizer_Optimize
{
    template<typename Value_t>
    void SynthesizeRule(
        const Rule& rule,
        CodeTree<Value_t>& tree,
        MatchInfo<Value_t>& info)
    {
        switch(rule.ruletype)
        {
            case ProduceNewTree:
            {
                tree.Become(
                    SynthesizeParam( ParamSpec_Extract<Value_t>(rule.repl_param_list, 0),
                                     info, false ) );
                break;
            }
            case ReplaceParams:
            default:
            {
                /* Delete the matched parameters from the source tree */
                std::vector<unsigned> list = info.GetMatchedParamIndexes();
                std::sort(list.begin(), list.end());
                for(size_t a=list.size(); a-->0; )
                    tree.DelParam( list[a] );

                /* Synthesize the replacement params */
                for(unsigned a=0; a < rule.repl_param_count; ++a)
                {
                    CodeTree<Value_t> nparam
                     = SynthesizeParam( ParamSpec_Extract<Value_t>(rule.repl_param_list, a),
                                        info, true );
                    tree.AddParamMove(nparam);
                }
                break;
            }
        }
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_Optimize
{
#define FP_INSTANTIATE(type) \
    template \
    void SynthesizeRule( \
        const Rule& rule, \
        CodeTree<type>& tree, \
        MatchInfo<type>& info);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
