#include <cmath>
#include <cassert>

#include "codetree.hh"
#include "optimize.hh"
#include "opcodename.hh"
#include "grammar.hh"
#include "extrasrc/fptypes.hh"

#include "consts.hh"
#include "fparser.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    using namespace FPoptimizer_CodeTree;

    #define FactorStack std::vector /* typedef*/

    const struct PowiMuliType
    {
        unsigned opcode_square;
        unsigned opcode_cumulate;
        unsigned opcode_invert;
        unsigned opcode_half;
        unsigned opcode_invhalf;
    } iseq_powi = {cSqr,cMul,cInv,cSqrt,cRSqrt},
      iseq_muli = {~unsigned(0), cAdd,cNeg, ~unsigned(0),~unsigned(0) };

    template<typename Value_t>
    Value_t ParsePowiMuli(
        const PowiMuliType& opcodes,
        const std::vector<unsigned>& ByteCode, size_t& IP,
        size_t limit,
        size_t factor_stack_base,
        FactorStack<Value_t>& stack)
    {
        Value_t result(1);
        while(IP < limit)
        {
            if(ByteCode[IP] == opcodes.opcode_square)
            {
                if(!isInteger(result)) break;
                result *= 2;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_invert)
            {
                result = -result;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_half)
            {
                if(result > Value_t(0) && isEvenInteger(result))
                    break;
                result *= Value_t(0.5);
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_invhalf)
            {
                if(result > Value_t(0) && isEvenInteger(result))
                    break;
                result *= Value_t(-0.5);
                ++IP;
                continue;
            }

            size_t dup_fetch_pos = IP;
            Value_t lhs(1);

            if(ByteCode[IP] == cFetch)
            {
                unsigned index = ByteCode[++IP];
                if(index < factor_stack_base
                || size_t(index-factor_stack_base) >= stack.size())
                {
                    // It wasn't a powi-fetch after all
                    IP = dup_fetch_pos;
                    break;
                }
                lhs = stack[index - factor_stack_base];
                // Note: ^This assumes that cFetch of recentmost
                //        is always converted into cDup.
                goto dup_or_fetch;
            }
            if(ByteCode[IP] == cDup)
            {
                lhs = result;
                goto dup_or_fetch;

            dup_or_fetch:
                stack.push_back(result);
                ++IP;
                Value_t subexponent = ParsePowiMuli
                    (opcodes,
                     ByteCode, IP, limit,
                     factor_stack_base, stack);
                if(IP >= limit || ByteCode[IP] != opcodes.opcode_cumulate)
                {
                    // It wasn't a powi-dup after all
                    IP = dup_fetch_pos;
                    break;
                }
                ++IP; // skip opcode_cumulate
                stack.pop_back();
                result += lhs*subexponent;
                continue;
            }
            break;
        }
        return result;
    }

    template<typename Value_t>
    Value_t ParsePowiSequence
        (const std::vector<unsigned>& ByteCode, size_t& IP,
         size_t limit, size_t factor_stack_base)
    {
        FactorStack<Value_t> stack;
        stack.push_back( Value_t(1) );
        return ParsePowiMuli(iseq_powi, ByteCode, IP, limit, factor_stack_base, stack);
    }

    template<typename Value_t>
    Value_t ParseMuliSequence
        (const std::vector<unsigned>& ByteCode, size_t& IP,
         size_t limit, size_t factor_stack_base)
    {
        FactorStack<Value_t> stack;
        stack.push_back( Value_t(1) );
        return ParsePowiMuli(iseq_muli, ByteCode, IP, limit, factor_stack_base, stack);
    }

    template<typename Value_t>
    class CodeTreeParserData
    {
    public:
        explicit CodeTreeParserData(bool k_powi)
            : stack(), clones(), keep_powi(k_powi) { }

        void Eat(size_t nparams, OPCODE opcode)
        {
            CodeTree<Value_t> newnode;
            newnode.SetOpcode(opcode);

            std::vector<CodeTree<Value_t> > params = Pop(nparams);
            newnode.SetParamsMove(params);

            if(!keep_powi)
            switch(opcode)
            {
                //        asinh: log(x + sqrt(x*x + 1))
                //cAsinh [x] -> cLog (cAdd x (cPow (cAdd (cPow x 2) 1) 0.5))
                // Note: ^ Replacement function refers to x twice

                //        acosh: log(x + sqrt(x*x - 1))
                //cAcosh [x] -> cLog (cAdd x (cPow (cAdd (cPow x 2) -1) 0.5))

                //        atanh: log( (1+x) / (1-x)) / 2
                //cAtanh [x] -> cMul (cLog (cMul (cAdd 1 x) (cPow (cAdd 1 (cMul -1 x)) -1))) 0.5

                //        asin: atan2(x, sqrt(1-x*x))
                //cAsin[x] -> cAtan2 [x (cPow [(cAdd 1 (cMul (cPow [x 2] -1)) 0.5])]

                //        acos: atan2(sqrt(1-x*x), x)
                //cAcos[x] -> cAtan2 [(cPow [(cAdd 1 (cMul (cPow [x 2] -1)) 0.5]) x]

                //     The hyperbolic functions themselves are:
                //        sinh: (exp(x)-exp(-x)) / 2  = exp(-x) * (exp(2*x)-1) / 2
                //cSinh [x] -> cMul 0.5 (cPow [CONSTANT_EI x]) (cAdd [-1 (cPow [CONSTANT_2E x])])

                //        cosh: (exp(x)+exp(-x)) / 2  = exp(-x) * (exp(2*x)+1) / 2
                //        cosh(-x) = cosh(x)
                //cCosh [x] -> cMul 0.5 (cPow [CONSTANT_EI x]) (cAdd [ 1 (cPow [CONSTANT_2E x])])

                //        tanh: sinh/cosh = (exp(2*x)-1) / (exp(2*x)+1)
                //cTanh [x] -> (cMul (cAdd {(cPow [CONSTANT_2E x]) -1}) (cPow [(cAdd {(cPow [CONSTANT_2E x]) 1}) -1]))
                case cTanh:
                {
                    CodeTree<Value_t> sinh, cosh;
                    sinh.SetOpcode(cSinh); sinh.AddParam(newnode.GetParam(0)); sinh.Rehash();
                    cosh.SetOpcode(cCosh); cosh.AddParamMove(newnode.GetParam(0)); cosh.Rehash();
                    CodeTree<Value_t> pow;
                    pow.SetOpcode(cPow);
                    pow.AddParamMove(cosh);
                    pow.AddParam(CodeTreeImmed(Value_t(-1)));
                    pow.Rehash();
                    newnode.SetOpcode(cMul);
                    newnode.SetParamMove(0, sinh);
                    newnode.AddParamMove(pow);
                    break;
                }

                //        tan: sin/cos
                //cTan [x] -> (cMul (cSin [x]) (cPow [(cCos [x]) -1]))
                case cTan:
                {
                    CodeTree<Value_t> sin, cos;
                    sin.SetOpcode(cSin); sin.AddParam(newnode.GetParam(0)); sin.Rehash();
                    cos.SetOpcode(cCos); cos.AddParamMove(newnode.GetParam(0)); cos.Rehash();
                    CodeTree<Value_t> pow;
                    pow.SetOpcode(cPow);
                    pow.AddParamMove(cos);
                    pow.AddParam(CodeTreeImmed(Value_t(-1)));
                    pow.Rehash();
                    newnode.SetOpcode(cMul);
                    newnode.SetParamMove(0, sin);
                    newnode.AddParamMove(pow);
                    break;
                }

                case cPow:
                {
                    const CodeTree<Value_t>& p0 = newnode.GetParam(0);
                    const CodeTree<Value_t>& p1 = newnode.GetParam(1);
                    if(p1.GetOpcode() == cAdd)
                    {
                        // convert x^(a + b) into x^a * x^b just so that
                        // some optimizations can be run on it.
                        // For instance, exp(log(x)*-61.1 + log(z)*-59.1)
                        // won't be changed into exp(log(x*z)*-61.1)*z^2
                        // unless we do this.
                        std::vector<CodeTree<Value_t> > mulgroup(p1.GetParamCount());
                        for(size_t a=0; a<p1.GetParamCount(); ++a)
                        {
                            CodeTree<Value_t> pow;
                            pow.SetOpcode(cPow);
                            pow.AddParam(p0);
                            pow.AddParam(p1.GetParam(a));
                            pow.Rehash();
                            mulgroup[a].swap(pow);
                        }
                        newnode.SetOpcode(cMul);
                        newnode.SetParamsMove(mulgroup);
                    }
                    break;
                }

                // Should we change sin(x) into cos(pi/2-x)
                //               or cos(x) into sin(pi/2-x)?
                //                        note: cos(x-pi/2) = cos(pi/2-x) = sin(x)
                //                        note: sin(x-pi/2) = -sin(pi/2-x) = -cos(x)
                default: break;
            }

            newnode.Rehash(!keep_powi);
        /*
            using namespace FPoptimizer_Grammar;
            bool recurse = false;
            while(ApplyGrammar(pack.glist[0], newnode, recurse)) // intermediate
            { //std::cout << "Rerunning 1\n";
                newnode.FixIncompleteHashes();
                recurse = true;
            }
        */
            FindClone(newnode, false);
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "POP " << nparams << ", " << FP_GetOpcodeName(opcode)
                      << "->" << FP_GetOpcodeName(newnode.GetOpcode())
                      << ": PUSH ";
            DumpTree(newnode);
            std::cout <<std::endl;
            DumpHashes(newnode);
        #endif
            stack.push_back(newnode);
        }

        void EatFunc(size_t nparams, OPCODE opcode, unsigned funcno)
        {
            CodeTree<Value_t> newnode = CodeTreeFuncOp<Value_t> (opcode, funcno);
            std::vector<CodeTree<Value_t> > params = Pop(nparams);
            newnode.SetParamsMove(params);
            newnode.Rehash(false);
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "POP " << nparams << ", PUSH ";
            DumpTree(newnode);
            std::cout << std::endl;
            DumpHashes(newnode);
        #endif
            FindClone(newnode);
            stack.push_back(newnode);
        }

        void AddConst(const Value_t& value)
        {
            CodeTree<Value_t> newnode = CodeTreeImmed(value);
            FindClone(newnode);
            Push(newnode);
        }

        void AddVar(unsigned varno)
        {
            CodeTree<Value_t> newnode = CodeTreeVar<Value_t>(varno);
            FindClone(newnode);
            Push(newnode);
        }

        void SwapLastTwoInStack()
        {
            stack[stack.size()-1].swap( stack[stack.size()-2] );
        }

        void Dup()
        {
            Fetch(stack.size()-1);
        }

        void Fetch(size_t which)
        {
            Push(stack[which]);
        }

        template<typename T>
        void Push(T tree)
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "PUSH ";
            DumpTree(tree);
            std::cout << std::endl;
            DumpHashes(tree);
        #endif
            stack.push_back(tree);
        }

        void PopNMov(size_t target, size_t source)
        {
            stack[target] = stack[source];
            stack.resize(target+1);
        }

        CodeTree<Value_t> PullResult()
        {
            clones.clear();
            CodeTree<Value_t> result(stack.back());
            stack.resize(stack.size()-1);
            return result;
        }
        std::vector<CodeTree<Value_t> > Pop(size_t n_pop)
        {
            std::vector<CodeTree<Value_t> > result(n_pop);
            for(unsigned n=0; n<n_pop; ++n)
                result[n].swap(stack[stack.size()-n_pop+n]);
        #ifdef DEBUG_SUBSTITUTIONS
            for(size_t n=n_pop; n-- > 0; )
            {
                std::cout << "POP ";
                DumpTree(result[n]);
                std::cout << std::endl;
                DumpHashes(result[n]);
            }
        #endif
            stack.resize(stack.size()-n_pop);
            return result;
        }

        size_t GetStackTop() const { return stack.size(); }
    private:
        void FindClone(CodeTree<Value_t> & /*tree*/, bool /*recurse*/ = true)
        {
            // Disabled: Causes problems in optimization when
            // the same subtree is included in logical and non-logical
            // contexts: optimizations applied to the logical one will
            // mess up the non-logical one.
            return;
            /*
            std::multimap<fphash_t, CodeTree>::const_iterator
                i = clones.lower_bound(tree.GetHash());
            for(; i != clones.end() && i->first == tree.GetHash(); ++i)
            {
                if(i->second.IsIdenticalTo(tree))
                    tree.Become(i->second);
            }
            if(recurse)
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                    FindClone(tree.GetParam(a));
            clones.insert(std::make_pair(tree.GetHash(), tree));
            */
        }
    private:
        std::vector<CodeTree<Value_t> > stack;
        std::multimap<fphash_t, CodeTree<Value_t> > clones;

        bool keep_powi;

    private:
        CodeTreeParserData(const CodeTreeParserData&);
        CodeTreeParserData& operator=(const CodeTreeParserData&);
    };

    template<typename Value_t>
    struct IfInfo
    {
        CodeTree<Value_t> condition;
        CodeTree<Value_t> thenbranch;
        size_t endif_location;

        IfInfo(): condition(), thenbranch(), endif_location() { }
    };
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void CodeTree<Value_t>::GenerateFrom(
        const typename FunctionParserBase<Value_t>::Data& fpdata,
        bool keep_powi)
    {
        std::vector<CodeTree<Value_t> > var_trees;
        var_trees.reserve(fpdata.mVariablesAmount);
        for(unsigned n=0; n<fpdata.mVariablesAmount; ++n)
        {
            var_trees.push_back( CodeTreeVar<Value_t> (n+VarBegin) );
        }
        GenerateFrom(fpdata,var_trees,keep_powi);
    }

    template<typename Value_t>
    void CodeTree<Value_t>::GenerateFrom(
        const typename FunctionParserBase<Value_t>::Data& fpdata,
        const std::vector<CodeTree>& var_trees,
        bool keep_powi)
    {
        const std::vector<unsigned>& ByteCode = fpdata.mByteCode;
        const std::vector<Value_t>&  Immed    = fpdata.mImmed;

        /*for(unsigned i=0; i<ByteCode.size(); ++i)
            fprintf(stderr, "by[%u/%u]=%u\n", i, (unsigned)ByteCode.size(), (unsigned) ByteCode[i]);*/

    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "ENTERS GenerateFrom()\n";
    #endif
        CodeTreeParserData<Value_t> sim(keep_powi);
        std::vector<IfInfo<Value_t> > if_stack;

        for(size_t IP=0, DP=0; ; ++IP)
        {
        after_powi:
            while(!if_stack.empty() &&
              (   // Normal If termination rule:
                  if_stack.back().endif_location == IP
                  // This rule matches when cJumps are threaded:
               || (IP < ByteCode.size() && ByteCode[IP] == cJump
                   && if_stack.back().thenbranch.IsDefined())
              ))
            {
                // The "else" of an "if" ends here
                CodeTree elsebranch = sim.PullResult();
                sim.Push(if_stack.back().condition);
                sim.Push(if_stack.back().thenbranch);
                sim.Push(elsebranch);
                sim.Eat(3, cIf);
                if_stack.pop_back();
            }
            if(IP >= ByteCode.size()) break;

            unsigned opcode = ByteCode[IP];
            if((opcode == cSqr || opcode == cDup
             || (opcode == cInv && !IsIntType<Value_t>::result)
             || opcode == cNeg
             || opcode == cSqrt || opcode == cRSqrt
             || opcode == cFetch))
            {
                // Parse a powi sequence
                size_t was_ip = IP;
                Value_t exponent = ParsePowiSequence<Value_t>(
                    ByteCode, IP, if_stack.empty() ? ByteCode.size() : if_stack.back().endif_location,
                    sim.GetStackTop()-1);
                if(exponent != Value_t(1.0))
                {
                    //std::cout << "Found exponent at " << was_ip << ": " << exponent << "\n";
                    sim.AddConst(exponent);
                    sim.Eat(2, cPow);
                    goto after_powi;
                }
                if(opcode == cDup
                || opcode == cFetch
                || opcode == cNeg)
                {
                    Value_t factor = ParseMuliSequence<Value_t>(
                        ByteCode, IP, if_stack.empty() ? ByteCode.size() : if_stack.back().endif_location,
                        sim.GetStackTop()-1);
                    if(factor != Value_t(1.0))
                    {
                        //std::cout << "Found factor at " << was_ip << ": " << factor << "\n";
                        sim.AddConst(factor);
                        sim.Eat(2, cMul);
                        goto after_powi;
                    }
                }
                IP = was_ip;
            }
            if(OPCODE(opcode) >= VarBegin)
            {
                unsigned index = opcode-VarBegin;
                /*std::fprintf(stderr, "IP=%u, opcode=%u, VarBegin=%u, ind=%u\n",
                    (unsigned)IP,(unsigned)opcode,(unsigned)VarBegin, (unsigned)index);*/
                sim.Push(var_trees[index]);
            }
            else
            {
                switch( OPCODE(opcode) )
                {
                    // Specials
                    case cIf:
                    case cAbsIf:
                    {
                        if_stack.resize(if_stack.size() + 1);
                        CodeTree res( sim.PullResult() );
                        if_stack.back().condition.swap( res );
                        if_stack.back().endif_location = ByteCode.size();
                        IP += 2; // dp,sp for elsebranch are irrelevant.
                        continue;
                    }
                    case cJump:
                    {
                        CodeTree res( sim.PullResult() );
                        if_stack.back().thenbranch.swap( res );
                        if_stack.back().endif_location = ByteCode[IP+1]+1;
                        IP += 2;
                        continue;
                    }
                    case cImmed:
                        sim.AddConst(Immed[DP++]);
                        break;
                    case cDup:
                        sim.Dup();
                        break;
                    case cNop:
                        break;
                    case cFCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        assert(funcno < fpdata.mFuncPtrs.size());
                        unsigned params = fpdata.mFuncPtrs[funcno].mParams;
                        sim.EatFunc(params, OPCODE(opcode), funcno);
                        break;
                    }
                    case cPCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        assert(funcno < fpdata.mFuncParsers.size());
                        const FunctionParserBase<Value_t>& p =
                            *fpdata.mFuncParsers[funcno].mParserPtr;
                        unsigned params = fpdata.mFuncParsers[funcno].mParams;

                        /* Inline the procedure call */
                        /* Works because cPCalls can never recurse */
                        std::vector<CodeTree> paramlist = sim.Pop(params);
                        CodeTree pcall_tree;
                        pcall_tree.GenerateFrom(*p.mData, paramlist);
                        sim.Push(pcall_tree);
                        break;
                    }
                    // Unary operators requiring special attention
                    case cInv:  // already handled by powi_opt
                        //sim.Eat(1, cInv);
                        //break;
                        sim.AddConst(1);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cDiv);
                        break;
                    case cNeg: // already handled by powi_opt
                        sim.Eat(1, cNeg);
                        break;
                        sim.AddConst(0);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cSub);
                        break;
                    case cSqr: // already handled by powi_opt
                        //sim.Eat(1, cSqr);
                        //break;
                        sim.AddConst(2);
                        sim.Eat(2, cPow);
                        break;
                    // Unary functions requiring special attention
                    case cSqrt: // already handled by powi_opt
                        sim.AddConst( Value_t(0.5) );
                        sim.Eat(2, cPow);
                        break;
                    case cRSqrt: // already handled by powi_opt
                        sim.AddConst( Value_t(-0.5) );
                        sim.Eat(2, cPow);
                        break;
                    case cCbrt:
                        sim.AddConst(Value_t(1) / Value_t(3));
                        sim.Eat(2, cPow);
                        break;
                    case cDeg:
                        sim.AddConst(fp_const_rad_to_deg<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cRad:
                        sim.AddConst(fp_const_deg_to_rad<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cExp:
                        if(keep_powi) goto default_function_handling;
                        sim.AddConst(fp_const_e<Value_t>());
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cPow);
                        break;
                    case cExp2: // from fpoptimizer
                        if(keep_powi) goto default_function_handling;
                        sim.AddConst(2.0);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cPow);
                        break;
                    case cCot:
                        sim.Eat(1, cTan);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cCsc:
                        sim.Eat(1, cSin);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cSec:
                        sim.Eat(1, cCos);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cInt: // int(x) = floor(x + 0.5)
                    #ifndef __x86_64
                        if(keep_powi) { sim.Eat(1, cInt); break; }
                    #endif
                        sim.AddConst( Value_t(0.5) );
                        sim.Eat(2, cAdd);
                        sim.Eat(1, cFloor);
                        break;
                    case cLog10:
                        sim.Eat(1, cLog);
                        sim.AddConst(fp_const_log10inv<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cLog2:
                        sim.Eat(1, cLog);
                        sim.AddConst(fp_const_log2inv<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cLog2by: // x y     -> log(x)*CONSTANT_L2I*y
                        sim.SwapLastTwoInStack();   // y x
                        sim.Eat(1, cLog);           // y log(x)
                        sim.AddConst(fp_const_log2inv<Value_t>()); // y log(x) CONSTANT_L2I
                        sim.Eat(3, cMul);           // y*log(x)*CONSTANT_L2I
                        break;
                    case cHypot: // x y -> sqrt(x*x + y*y)
                        sim.AddConst(2);
                        sim.Eat(2, cPow); // x y^2
                        sim.SwapLastTwoInStack();
                        sim.AddConst(2);
                        sim.Eat(2, cPow); // y^2 x^2
                        sim.Eat(2, cAdd); // y^2 + x^2
                        sim.AddConst( Value_t(0.5) );
                        sim.Eat(2, cPow); // (y^2 + x^2)^0.5
                        break;
                    case cSinCos:
                        sim.Dup();
                        sim.Eat(1, cSin);
                        sim.SwapLastTwoInStack();
                        sim.Eat(1, cCos);
                        break;
                    case cSinhCosh:
                        sim.Dup();
                        sim.Eat(1, cSinh);
                        sim.SwapLastTwoInStack();
                        sim.Eat(1, cCosh);
                        break;
                    //case cLog:
                    //    sim.Eat(1, cLog2);
                    //    sim.AddConst(fp_const_log2<Value_t>());
                    //    sim.Eat(2, cMul);
                    //    break;
                    // Binary operators requiring special attention
                    case cRSub: // from fpoptimizer
                        sim.SwapLastTwoInStack();
                        // Passthru to cSub
                    case cSub:
                        if(keep_powi) { sim.Eat(2, cSub); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cMul); // -x is x*-1
                        sim.Eat(2, cAdd); // Minus is negative adding
                        break;
                    case cRDiv: // from fpoptimizer
                        sim.SwapLastTwoInStack();
                        // Passthru to cDiv
                    case cDiv:
                        if(keep_powi || IsIntType<Value_t>::result)
                        {
                            sim.Eat(2, cDiv);
                            break;
                        }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow); // 1/x is x^-1
                        sim.Eat(2, cMul); // Divide is inverse multiply
                        break;
                    // Binary operators not requiring special attention
                    case cAdd: case cMul:
                    case cMod: case cPow:
                    case cEqual: case cLess: case cGreater:
                    case cNEqual: case cLessOrEq: case cGreaterOrEq:
                    case cAnd: case cOr:
                    case cAbsAnd: case cAbsOr:
                        sim.Eat(2, OPCODE(opcode));
                        break;
                    // Unary operators not requiring special attention
                    case cNot:
                    case cNotNot: // from fpoptimizer
                    case cAbsNot:
                    case cAbsNotNot:
                        sim.Eat(1, OPCODE(opcode));
                        break;
                    // Special opcodes generated by fpoptimizer itself
                    case cFetch:
                        sim.Fetch(ByteCode[++IP]);
                        break;
                    case cPopNMov:
                    {
                        unsigned stackOffs_target = ByteCode[++IP];
                        unsigned stackOffs_source = ByteCode[++IP];
                        sim.PopNMov(stackOffs_target, stackOffs_source);
                        break;
                    }

                    default:
                    default_function_handling:;
                        unsigned funcno = opcode-cAbs;
                        assert(funcno < FUNC_AMOUNT);
                        const FuncDefinition& func = Functions[funcno];
                        sim.Eat(func.params, OPCODE(opcode));
                        break;
                }
            }
        }
        Become(sim.PullResult());
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Produced tree:\n";
        DumpTreeWithIndent(*this);
    #endif
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template \
    void CodeTree<type>::GenerateFrom( \
        const FunctionParserBase<type>::Data& fpdata, \
        bool keep_powi);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
