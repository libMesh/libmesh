#include <cmath>
#include <list>
#include <cassert>

#include "codetree.hh"
#include "extrasrc/fptypes.hh"
#include "consts.hh"
#include "optimize.hh"
#include "bytecodesynth.hh"

//#include "grammar.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    using namespace FPoptimizer_CodeTree;

    template<typename Value_t>
    bool AssembleSequence(
                  const CodeTree<Value_t>& tree, long count,
                  const FPoptimizer_ByteCode::SequenceOpCode<Value_t>& sequencing,
                  FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
                  size_t max_bytecode_grow_length);

    /*
    Trigonomic operations are expensive.
    If we need to synthesize a tan(x), we could
    first check if we can synthesize it by utilizing
    some expressions that are already in the stack.
    Such as: 1 / cot(x)
             sin(x) / cos(x)
             sin(x) * sec(x)
             (1/csc(x)) / cos(x)
             sec(x) / csc(x)
    This table lists the instructions for building
    each of these functions from components that
    might already exist in the stack.

    It is up to CSE to make it possible for this
    opportunity to arise as often as possible.
    */
    static const struct SinCosTanDataType
    {
        OPCODE whichopcode;
        OPCODE inverse_opcode;
        enum { nominator,
               denominator,
               inverse_nominator,
               inverse_denominator };
        OPCODE codes[4];
    } SinCosTanData[12] =
    {
        { cTan, cCot, { cSin,cCos, cCsc,cSec } },
        { cCot, cCot, { cCos,cSin, cSec,cCsc } },
        // tan(x) = 1/cot(x)
        //        = sin(x) / cos(x)
        //        = sin(x) * sec(x)
        //        = (1/csc(x)) / cos(x)
        //        = sec(x) / csc(x)

        { cCos, cSec, { cSin,cTan, cCsc,cCot } },
        { cSec, cCos, { cTan,cSin, cCot,cCsc } },
        // cos(x) = 1/sec(x) = sin(x) / tan(x)

        { cSin, cCsc, { cCos,cCot, cSec,cTan } },
        { cCsc, cSin, { cCot,cCos, cTan,cSec } },
        // sin(x) = 1/csc(x) = cos(x)/cot(x)

        { cTanh, cNop, { cSinh,cCosh, cNop,cNop } },
        { cSinh, cNop, { cTanh,cNop,  cNop,cCosh } },
        { cCosh, cNop, { cSinh,cTanh, cNop,cNop } },
        { cNop, cTanh, { cCosh,cSinh, cNop,cNop } },
        { cNop, cSinh, { cNop, cTanh, cCosh,cNop } },
        { cNop, cCosh, { cTanh,cSinh, cNop,cNop } }
    };
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void CodeTree<Value_t>::SynthesizeByteCode(
        std::vector<unsigned>& ByteCode,
        std::vector<Value_t>&   Immed,
        size_t& stacktop_max)
    {
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Making bytecode for:\n";
        DumpTreeWithIndent(*this);
    #endif
        while(RecreateInversionsAndNegations())
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "One change issued, produced:\n";
            DumpTreeWithIndent(*this);
        #endif
            FixIncompleteHashes();

            using namespace FPoptimizer_Optimize;
            using namespace FPoptimizer_Grammar;
            const void* g = (const void*)&grammar_optimize_recreate;
            while(ApplyGrammar(*(const Grammar*)g, *this))
                {   FixIncompleteHashes();
                }
        }
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Actually synthesizing, after recreating inv/neg:\n";
        DumpTreeWithIndent(*this);
    #endif

        FPoptimizer_ByteCode::ByteCodeSynth<Value_t> synth;

        /* Then synthesize the actual expression */
        SynthesizeByteCode(synth, false);
        /* The "false" parameters tells SynthesizeByteCode
         * that at the outermost synthesizing level, it does
         * not matter if leftover temps are left in the stack.
         */
        synth.Pull(ByteCode, Immed, stacktop_max);
    }

    template<typename Value_t>
    void CodeTree<Value_t>::SynthesizeByteCode(
        FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
        bool MustPopTemps) const
    {
        // If the synth can already locate our operand in the stack,
        // never mind synthesizing it again, just dup it.
        if(synth.FindAndDup(*this))
        {
            return;
        }

        for(size_t a=0; a<12; ++a)
        {
            const SinCosTanDataType& data = SinCosTanData[a];
            if(data.whichopcode != cNop)
            {
                if(GetOpcode() != data.whichopcode) continue;

                CodeTree invtree;
                invtree.SetParams(GetParams());
                invtree.SetOpcode( data.inverse_opcode );
                invtree.Rehash(false);
                if(synth.FindAndDup(invtree))
                {
                    synth.AddOperation(cInv,1,1);
                    synth.StackTopIs(*this);
                    return;
                }
            }
            else
            {
                // cNop indicates that there's no dedicated
                // opcode that indicates an inverted function.
                // For example, we have inv(cosh(x)).
                if(GetOpcode() != cInv) continue;
                if(GetParam(0).GetOpcode() != data.inverse_opcode) continue;
                if(synth.FindAndDup(GetParam(0)))
                {
                    synth.AddOperation(cInv,1,1);
                    synth.StackTopIs(*this);
                    return;
                }
            }

            // Check which trees we can find
            size_t   found[4];
            for(size_t b=0; b<4; ++b)
            {
                CodeTree tree;
                if(data.codes[b] == cNop)
                {
                    tree.SetOpcode(cInv);
                    CodeTree subtree;
                    subtree.SetParams(GetParams());
                    subtree.SetOpcode(data.codes[b ^ 2]);
                    subtree.Rehash(false);
                    tree.AddParamMove(subtree);
                }
                else
                {
                    tree.SetParams(GetParams());
                    tree.SetOpcode(data.codes[b]);
                }
                tree.Rehash(false);
                found[b] = synth.FindPos(tree);
            }

            if(found[ data.nominator ]   != ~size_t(0)
            && found[ data.denominator ] != ~size_t(0))
            {
                // tan(x) = sin(x) / cos(x)
                synth.DoDup( found[data.nominator] );
                synth.DoDup( found[data.denominator] );
                synth.AddOperation(cDiv,2,1);
                synth.StackTopIs(*this);
                return;
            }

            if(found[ data.nominator ]           != ~size_t(0)
            && found[ data.inverse_denominator ] != ~size_t(0))
            {
                // tan(x) = sin(x) * sec(x)
                synth.DoDup( found[data.nominator] );
                synth.DoDup( found[data.inverse_denominator] );
                synth.AddOperation(cMul,2,1);
                synth.StackTopIs(*this);
                return;
            }

            if(found[ data.inverse_nominator ]   != ~size_t(0)
            && found[ data.inverse_denominator ] != ~size_t(0))
            {
                // tan(x) = 1 / (csc(x) / sec(x)) = sec(x) / csc(x)
                synth.DoDup( found[data.inverse_nominator] );
                synth.DoDup( found[data.inverse_denominator] );
                synth.AddOperation(cRDiv,2,1);
                synth.StackTopIs(*this);
                return;
            }

            if(found[ data.inverse_nominator ]   != ~size_t(0)
            && found[ data.denominator ]         != ~size_t(0))
            {
                // tan(x) = 1 / csc(x) / cos(x) = 1 / (csc(x) * cos(x))
                synth.DoDup( found[data.inverse_nominator] );
                synth.DoDup( found[data.denominator] );
                synth.AddOperation(cMul,2,1);
                synth.AddOperation(cInv,1,1);
                synth.StackTopIs(*this);
                return;
            }
        }

        size_t n_subexpressions_synthesized = SynthCommonSubExpressions(synth);

        switch(GetOpcode())
        {
            case VarBegin:
                synth.PushVar(GetVar());
                break;
            case cImmed:
                synth.PushImmed(GetImmed());
                break;
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd:
            case cOr:
            case cAbsAnd:
            case cAbsOr:
            {
                if(GetOpcode() == cMul) // Special treatment for cMul sequences
                {
                    // If the paramlist contains an Immed, and that Immed
                    // fits in a long-integer, try to synthesize it
                    // as add-sequences instead.
                    bool did_muli = false;
                    for(size_t a=0; a<GetParamCount(); ++a)
                    {
                        if(GetParam(a).IsImmed() && isLongInteger(GetParam(a).GetImmed()))
                        {
                            long value = makeLongInteger(GetParam(a).GetImmed());

                            CodeTree tmp(*this, typename CodeTree::CloneTag());
                            tmp.DelParam(a);
                            tmp.Rehash();
                            if(AssembleSequence(
                                tmp, value,
                                FPoptimizer_ByteCode::SequenceOpcodes<Value_t>::AddSequence,
                                synth,
                                MAX_MULI_BYTECODE_LENGTH))
                            {
                                did_muli = true;
                                break;
                            }
                        }
                    }
                    if(did_muli)
                        break; // done
                }

                // If any of the params is currently a copy of
                // the stack topmost item, treat it first.
                int n_stacked = 0;
                std::vector<bool> done( GetParamCount() , false );
                CodeTree synthed_tree;
                synthed_tree.SetOpcode(GetOpcode());
                for(;;)
                {
                    bool found = false;
                    for(size_t a=0; a<GetParamCount(); ++a)
                    {
                        if(done[a]) continue;
                        if(synth.IsStackTop(GetParam(a)))
                        {
                            found = true;
                            done[a] = true;
                            GetParam(a).SynthesizeByteCode(synth);
                            synthed_tree.AddParam(GetParam(a));
                            if(++n_stacked > 1)
                            {
                                // Cumulate at the earliest opportunity.
                                synth.AddOperation(GetOpcode(), 2); // stack state: -2+1 = -1
                                synthed_tree.Rehash(false);
                                synth.StackTopIs(synthed_tree);
                                n_stacked = n_stacked - 2 + 1;
                            }
                        }
                    }
                    if(!found) break;
                }

                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    if(done[a]) continue;
                    GetParam(a).SynthesizeByteCode(synth);
                    synthed_tree.AddParam(GetParam(a));
                    if(++n_stacked > 1)
                    {
                        // Cumulate at the earliest opportunity.
                        synth.AddOperation(GetOpcode(), 2); // stack state: -2+1 = -1
                        synthed_tree.Rehash(false);
                        synth.StackTopIs(synthed_tree);
                        n_stacked = n_stacked - 2 + 1;
                    }
                }
                if(n_stacked == 0)
                {
                    // Uh, we got an empty cAdd/cMul/whatever...
                    // Synthesize a default value.
                    // This should never happen.
                    switch(GetOpcode())
                    {
                        case cAdd:
                        case cOr:
                        case cAbsOr:
                            synth.PushImmed(0);
                            break;
                        case cMul:
                        case cAnd:
                        case cAbsAnd:
                            synth.PushImmed(1);
                            break;
                        case cMin:
                        case cMax:
                            //synth.PushImmed(NaN);
                            synth.PushImmed(0);
                            break;
                        default:
                            break;
                    }
                    ++n_stacked;
                }
                assert(n_stacked == 1);
                break;
            }
            case cPow:
            {
                const CodeTree& p0 = GetParam(0);
                const CodeTree& p1 = GetParam(1);

                if(!p1.IsImmed()
                || !isLongInteger(p1.GetImmed())
                || !AssembleSequence( /* Optimize integer exponents */
                        p0, makeLongInteger(p1.GetImmed()),
                        FPoptimizer_ByteCode::SequenceOpcodes<Value_t>::MulSequence,
                        synth,
                        MAX_POWI_BYTECODE_LENGTH)
                  )
                {
                    p0.SynthesizeByteCode(synth);
                    p1.SynthesizeByteCode(synth);
                    synth.AddOperation(GetOpcode(), 2); // Create a vanilla cPow.
                }
                break;
            }
            case cIf:
            case cAbsIf:
            {
                // Assume that the parameter count is 3 as it should.
                typename FPoptimizer_ByteCode::ByteCodeSynth<Value_t>::IfData ifdata;

                GetParam(0).SynthesizeByteCode(synth); // expression

                synth.SynthIfStep1(ifdata, GetOpcode());

                GetParam(1).SynthesizeByteCode(synth); // true branch

                synth.SynthIfStep2(ifdata);

                GetParam(2).SynthesizeByteCode(synth); // false branch

                synth.SynthIfStep3(ifdata);
                break;
            }
            case cFCall:
            case cPCall:
            {
                // Assume that the parameter count is as it should.
                for(size_t a=0; a<GetParamCount(); ++a)
                    GetParam(a).SynthesizeByteCode(synth);
                synth.AddOperation(GetOpcode(), (unsigned) GetParamCount());
                synth.AddOperation(0x80000000u | GetFuncNo(), 0, 0);
                break;
            }
            default:
            {
                // Assume that the parameter count is as it should.
                for(size_t a=0; a<GetParamCount(); ++a)
                    GetParam(a).SynthesizeByteCode(synth);
                synth.AddOperation(GetOpcode(), (unsigned) GetParamCount());
                break;
            }
        }

        // Tell the synthesizer which tree was just produced in the stack
        synth.StackTopIs(*this);

        // If we added subexpressions, peel them off the stack now
        if(MustPopTemps && n_subexpressions_synthesized > 0)
        {
            size_t top = synth.GetStackTop();
            synth.DoPopNMov(top-1-n_subexpressions_synthesized, top-1);
        }
    }
}

namespace
{
    template<typename Value_t>
    bool AssembleSequence(
        const CodeTree<Value_t>& tree, long count,
        const FPoptimizer_ByteCode::SequenceOpCode<Value_t>& sequencing,
        FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
        size_t max_bytecode_grow_length)
    {
        if(count != 0)
        {
            FPoptimizer_ByteCode::ByteCodeSynth<Value_t> backup = synth;

            tree.SynthesizeByteCode(synth);

            // Ignore the size generated by subtree
            size_t bytecodesize_backup = synth.GetByteCodeSize();

            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);

            size_t bytecode_grow_amount = synth.GetByteCodeSize() - bytecodesize_backup;
            if(bytecode_grow_amount > max_bytecode_grow_length)
            {
                synth = backup;
                return false;
            }
            return true;
        }
        else
        {
            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);
            return true;
        }
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template void CodeTree<type>::SynthesizeByteCode( \
        std::vector<unsigned>& ByteCode, \
        std::vector<type>&   Immed, \
        size_t& stacktop_max);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
