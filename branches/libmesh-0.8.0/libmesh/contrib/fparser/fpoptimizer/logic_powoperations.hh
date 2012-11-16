#include "codetree.hh"

#include <limits>

/****
#ifdef _MSC_VER
#include <float.h>
#define isinf(x) (!_finite(x))
#endif
*/

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /**************************************/
    /* OPERATIONS DONE TO POW()           */
    /**************************************/

    template<typename Value_t>
    int maxFPExponent()
    {
        return std::numeric_limits<Value_t>::max_exponent;
    }

    template<typename Value_t>
    bool fPExponentIsTooLarge(Value_t base, Value_t exponent)
    {
        if(base < Value_t(0)) return true;
        if(fp_equal(base, Value_t(0)) || fp_equal(base, Value_t(1)))
            return false;
        return exponent >= Value_t(maxFPExponent<Value_t>()) / fp_log2(base);
    }

    template<typename Value_t>
    bool ConstantFolding_PowOperations(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cPow);

        if(tree.GetParam(0).IsImmed()
        && tree.GetParam(1).IsImmed())
        {
            Value_t const_value = fp_pow(tree.GetParam(0).GetImmed(),
                                        tree.GetParam(1).GetImmed());
            tree.ReplaceWithImmed(const_value);
            return false;
        }
        if(tree.GetParam(1).IsImmed()
        && fp_equal(tree.GetParam(1).GetImmed(), Value_t(1)))
        {
            // Used to be: float(getimmed()) == 1.0
            // Conversion through a float type value gets rid of
            // awkward abs(x)^1 generated from exp(log(x^6)/6),
            // without sacrificing as much precision as fp_equal() does.
            // x^1 = x
            tree.Become(tree.GetParam(0));
            return true; // rerun optimization (opcode changed)
        }
        if(tree.GetParam(0).IsImmed()
        && fp_equal(tree.GetParam(0).GetImmed(), Value_t(1)))
        {
            // 1^x = 1
            tree.ReplaceWithImmed(1);
            return false;
        }

        // 5^(20*x) = (5^20)^x
        if(tree.GetParam(0).IsImmed()
        && tree.GetParam(1).GetOpcode() == cMul)
        {
            bool changes = false;
            Value_t base_immed = tree.GetParam(0).GetImmed();
            CodeTree<Value_t> mulgroup = tree.GetParam(1);
            for(size_t a=mulgroup.GetParamCount(); a-->0; )
                if(mulgroup.GetParam(a).IsImmed())
                {
                    Value_t imm = mulgroup.GetParam(a).GetImmed();
                    //if(imm >= 0.0)
                    {
                        /****
                        Value_t new_base_immed = fp_pow(base_immed, imm);
                        if(isinf(new_base_immed)
                        || fp_equal(new_base_immed, Value_t(0)))
                        {
                            // It produced an infinity. Do not change.
                            break;
                        }
                        */
                        if(fPExponentIsTooLarge(base_immed, imm))
                            break;

                        Value_t new_base_immed = fp_pow(base_immed, imm);
                        if(fp_equal(new_base_immed, Value_t(0)))
                            break;

                        if(!changes)
                        {
                            changes = true;
                            mulgroup.CopyOnWrite();
                        }
                        base_immed = new_base_immed;
                        mulgroup.DelParam(a);
                        break; //
                    }
                }
            if(changes)
            {
                mulgroup.Rehash();
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before pow-mul change: "; DumpTree(tree);
                std::cout << "\n";
            #endif
                tree.GetParam(0).Become(CodeTreeImmed<Value_t> (base_immed));
                tree.GetParam(1).Become(mulgroup);
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After pow-mul change: "; DumpTree(tree);
                std::cout << "\n";
            #endif
            }
        }
        // (x*20)^2 = x^2 * 20^2
        if(tree.GetParam(1).IsImmed()
        && tree.GetParam(0).GetOpcode() == cMul)
        {
            Value_t exponent_immed = tree.GetParam(1).GetImmed();
            Value_t factor_immed   = 1.0;
            bool changes = false;
            CodeTree<Value_t>& mulgroup = tree.GetParam(0);
            for(size_t a=mulgroup.GetParamCount(); a-->0; )
                if(mulgroup.GetParam(a).IsImmed())
                {
                    Value_t imm = mulgroup.GetParam(a).GetImmed();
                    //if(imm >= 0.0)
                    {
                        /****
                        Value_t new_factor_immed = fp_pow(imm, exponent_immed);
                        if(isinf(new_factor_immed)
                        || fp_equal(new_factor_immed, Value_t(0)))
                        {
                            // It produced an infinity. Do not change.
                            break;
                        }
                        */
                        if(fPExponentIsTooLarge(imm, exponent_immed))
                            break;

                        Value_t new_factor_immed = fp_pow(imm, exponent_immed);
                        if(fp_equal(new_factor_immed, Value_t(0)))
                            break;

                        if(!changes)
                        {
                            changes = true;
                            mulgroup.CopyOnWrite();
                        }
                        factor_immed *= new_factor_immed;
                        mulgroup.DelParam(a);
                        break; //
                    }
                }
            if(changes)
            {
                mulgroup.Rehash();
                CodeTree<Value_t> newpow;
                newpow.SetOpcode(cPow);
                newpow.SetParamsMove(tree.GetParams());
                newpow.Rehash(false);
                tree.SetOpcode(cMul);
                tree.AddParamMove(newpow);
                tree.AddParam( CodeTreeImmed<Value_t>(factor_immed) );
                return true; // rerun optimization (opcode changed)
            }
        }

        // (x^3)^2 = x^6
        // NOTE: If 3 is even and 3*2 is not, x must be changed to abs(x).
        if(tree.GetParam(0).GetOpcode() == cPow
        && tree.GetParam(1).IsImmed()
        && tree.GetParam(0).GetParam(1).IsImmed())
        {
            Value_t a = tree.GetParam(0).GetParam(1).GetImmed();
            Value_t b = tree.GetParam(1).GetImmed();
            Value_t c = a * b; // new exponent
            if(isEvenInteger(a) // a is an even int?
            && !isEvenInteger(c)) // c is not?
            {
                CodeTree<Value_t> newbase;
                newbase.SetOpcode(cAbs);
                newbase.AddParam(tree.GetParam(0).GetParam(0));
                newbase.Rehash();
                tree.SetParamMove(0, newbase);
            }
            else
                tree.SetParam(0, tree.GetParam(0).GetParam(0));
            tree.SetParam(1, CodeTreeImmed<Value_t>(c));
        }
        return false; // No changes that require a rerun
    }
}
