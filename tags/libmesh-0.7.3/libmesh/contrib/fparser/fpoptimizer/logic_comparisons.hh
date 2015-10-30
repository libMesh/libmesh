#include "codetree.hh"

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /*****************************************/
    /* RANGE-BASED TREATMENTS TO COMPARISONS */
    /*****************************************/

    struct RangeComparisonData
    {
        enum Decision
        {
            MakeFalse=0,
            MakeTrue=1,
            MakeNEqual=2,
            MakeEqual=3,
            MakeNotNotP0=4,
            MakeNotNotP1=5,
            MakeNotP0=6,
            MakeNotP1=7,
            Unchanged=8
        };
        enum WhatDoWhenCase
        {
            Never =0,
            Eq0   =1, // val==0
            Eq1   =2, // val==1
            Gt0Le1=3, // val>0 && val<=1
            Ge0Lt1=4  // val>=0 && val<1
        };

        Decision if_identical; // What to do when operands are identical
        Decision if_always[4]; // What to do if Always <, <=, >, >=
        struct { Decision what : 4; WhatDoWhenCase when : 4; }
            p0_logical_a, p1_logical_a,
            p0_logical_b, p1_logical_b;

        template<typename Value_t>
        Decision Analyze(const CodeTree<Value_t>& a, const CodeTree<Value_t>& b) const
        {
            if(a.IsIdenticalTo(b))
                return if_identical;

            range<Value_t> p0 = CalculateResultBoundaries(a);
            range<Value_t> p1 = CalculateResultBoundaries(b);
            if(p0.max.known && p1.min.known)
            {
                if(p0.max.val <  p1.min.val && if_always[0] != Unchanged)
                    return if_always[0]; // p0 < p1
                if(p0.max.val <= p1.min.val && if_always[1] != Unchanged)
                    return if_always[1]; // p0 <= p1
            }
            if(p0.min.known && p1.max.known)
            {
                if(p0.min.val >  p1.max.val && if_always[2] != Unchanged)
                    return if_always[2]; // p0 > p1
                if(p0.min.val >= p1.max.val && if_always[3] != Unchanged)
                    return if_always[3]; // p0 >= p1
            }

            if(IsLogicalValue(a))
            {
                if(p0_logical_a.what != Unchanged)
                    if(TestCase(p0_logical_a.when, p1)) return p0_logical_a.what;
                if(p0_logical_b.what != Unchanged)
                    if(TestCase(p0_logical_b.when, p1)) return p0_logical_b.what;
            }
            if(IsLogicalValue(b))
            {
                if(p1_logical_a.what != Unchanged)
                    if(TestCase(p1_logical_a.when, p0)) return p1_logical_a.what;
                if(p1_logical_b.what != Unchanged)
                    if(TestCase(p1_logical_b.when, p0)) return p1_logical_b.what;
            }
            return Unchanged;
        }

        template<typename Value_t>
        static bool TestCase(WhatDoWhenCase when, const range<Value_t>& p)
        {
            if(!p.min.known || !p.max.known) return false;
            switch(when)
            {
                case Eq0: return p.min.val==Value_t(0.0) && p.max.val==p.min.val;
                case Eq1: return p.min.val==Value_t(1.0) && p.max.val==p.max.val;
                case Gt0Le1: return p.min.val>Value_t(0) && p.max.val<=Value_t(1);
                case Ge0Lt1: return p.min.val>=Value_t(0) && p.max.val<Value_t(1);
                default:;
            }
            return false;
        }
    };

    namespace RangeComparisonsData
    {
        static const RangeComparisonData Data[6] =
        {
            // cEqual:
            // Case:      p0 == p1  Antonym: p0 != p1
            // Synonym:   p1 == p0  Antonym: p1 != p0
            { RangeComparisonData::MakeTrue,  // If identical: always true
              {RangeComparisonData::MakeFalse,  // If Always p0 < p1: always false
               RangeComparisonData::Unchanged,
               RangeComparisonData::MakeFalse,  // If Always p0 > p1: always false
               RangeComparisonData::Unchanged},
             // NotNot(p0) if p1==1    NotNot(p1) if p0==1
             //    Not(p0) if p1==0       Not(p1) if p0==0
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotP0, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotP1, RangeComparisonData::Eq0}
            },
            // cNEqual:
            // Case:      p0 != p1  Antonym: p0 == p1
            // Synonym:   p1 != p0  Antonym: p1 == p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeTrue,  // If Always p0 < p1: always true
               RangeComparisonData::Unchanged,
               RangeComparisonData::MakeTrue,  // If Always p0 > p1: always true
               RangeComparisonData::Unchanged},
             // NotNot(p0) if p1==0    NotNot(p1) if p0==0
             //    Not(p0) if p1==1       Not(p1) if p0==1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotP0, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotP1, RangeComparisonData::Eq1}
            },
            // cLess:
            // Case:      p0 < p1   Antonym: p0 >= p1
            // Synonym:   p1 > p0   Antonym: p1 <= p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeTrue,  // If Always p0  < p1: always true
               RangeComparisonData::MakeNEqual,
               RangeComparisonData::MakeFalse, // If Always p0 > p1: always false
               RangeComparisonData::MakeFalse},// If Always p0 >= p1: always false
             // Not(p0)   if p1>0 & p1<=1    --   NotNot(p1) if p0>=0 & p0<1
              {RangeComparisonData::MakeNotP0,    RangeComparisonData::Gt0Le1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cLessOrEq:
            // Case:      p0 <= p1  Antonym: p0 > p1
            // Synonym:   p1 >= p0  Antonym: p1 < p0
            { RangeComparisonData::MakeTrue,   // If identical: always true
              {RangeComparisonData::Unchanged, // If Always p0  < p1: ?
               RangeComparisonData::MakeTrue,  // If Always p0 <= p1: always true
               RangeComparisonData::MakeFalse, // If Always p0  > p1: always false
               RangeComparisonData::MakeEqual},// If Never  p0  < p1:  use cEqual
             // Not(p0)    if p1>=0 & p1<1   --   NotNot(p1) if p0>0 & p0<=1
              {RangeComparisonData::MakeNotP0,    RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Gt0Le1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cGreater:
            // Case:      p0 >  p1  Antonym: p0 <= p1
            // Synonym:   p1 <  p0  Antonym: p1 >= p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeFalse, // If Always p0  < p1: always false
               RangeComparisonData::MakeFalse, // If Always p0 <= p1: always false
               RangeComparisonData::MakeTrue,  // If Always p0  > p1: always true
               RangeComparisonData::MakeNEqual},
             // NotNot(p0) if p1>=0 & p1<1   --   Not(p1)   if p0>0 & p0<=1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::MakeNotP1,    RangeComparisonData::Gt0Le1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cGreaterOrEq:
            // Case:      p0 >= p1  Antonym: p0 < p1
            // Synonym:   p1 <= p0  Antonym: p1 > p0
            { RangeComparisonData::MakeTrue,   // If identical: always true
              {RangeComparisonData::MakeFalse, // If Always p0  < p1: always false
               RangeComparisonData::MakeEqual, // If Always p0 >= p1: always true
               RangeComparisonData::Unchanged, // If always p0  > p1: ?
               RangeComparisonData::MakeTrue}, // If Never  p0  > p1:  use cEqual
             // NotNot(p0) if p1>0 & p1<=1   --   Not(p1)    if p0>=0 & p0<1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Gt0Le1},
              {RangeComparisonData::MakeNotP1,    RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            }
        };
    }

    template<typename Value_t>
    bool ConstantFolding_Comparison(CodeTree<Value_t>& tree)
    {
        using namespace RangeComparisonsData;

        assert(tree.GetOpcode() >= cEqual && tree.GetOpcode() <= cGreaterOrEq);

        switch(Data[tree.GetOpcode()-cEqual].
            Analyze(tree.GetParam(0), tree.GetParam(1)))
        {
            case RangeComparisonData::MakeFalse:
                tree.ReplaceWithImmed(0); return true;
            case RangeComparisonData::MakeTrue:
                tree.ReplaceWithImmed(1); return true;
            case RangeComparisonData::MakeEqual:  tree.SetOpcode(cEqual); return true;
            case RangeComparisonData::MakeNEqual: tree.SetOpcode(cNEqual); return true;
            case RangeComparisonData::MakeNotNotP0: tree.SetOpcode(cNotNot); tree.DelParam(1); return true;
            case RangeComparisonData::MakeNotNotP1: tree.SetOpcode(cNotNot); tree.DelParam(0); return true;
            case RangeComparisonData::MakeNotP0: tree.SetOpcode(cNot); tree.DelParam(1); return true;
            case RangeComparisonData::MakeNotP1: tree.SetOpcode(cNot); tree.DelParam(0); return true;
            case RangeComparisonData::Unchanged:;
        }

        // Any reversible functions:
        //   sin(x)  -> ASIN: Not doable, x can be cyclic
        //   asin(x) -> SIN: doable.
        //                   Invalid combinations are caught by
        //                   range-estimation. Threshold is at |pi/2|.
        //   acos(x) -> COS: doable.
        //                   Invalid combinations are caught by
        //                   range-estimation. Note that though
        //                   the range is contiguous, it is direction-flipped.
        //    log(x) -> EXP: no problem
        //   exp2, exp10: Converted to cPow, done by grammar.
        //   atan(x) -> TAN: doable.
        //                   Invalid combinations are caught by
        //                   range-estimation. Threshold is at |pi/2|.
        //   sinh(x) -> ASINH: no problem
        //   tanh(x) -> ATANH: no problem, but atanh is limited to -1..1
        //                     Invalid combinations are caught by
        //                     range-estimation, but the exact value
        //                     of 1.0 still needs checking, because
        //                     it involves infinity.
        if(tree.GetParam(1).IsImmed())
            switch(tree.GetParam(0).GetOpcode())
            {
                case cAsin:
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_sin(tree.GetParam(1).GetImmed())));
                    return true;
                case cAcos:
                    // -1..+1 --> pi..0 (polarity-flipping)
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_cos(tree.GetParam(1).GetImmed())));
                    tree.SetOpcode( tree.GetOpcode()==cLess ? cGreater
                                  : tree.GetOpcode()==cLessOrEq ? cGreaterOrEq
                                  : tree.GetOpcode()==cGreater ? cLess
                                  : tree.GetOpcode()==cGreaterOrEq ? cLessOrEq
                                  : tree.GetOpcode() );
                    return true;
                case cAtan:
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_tan(tree.GetParam(1).GetImmed())));
                    return true;
                case cLog:
                    // Different logarithms have a constant-multiplication,
                    // which is no problem.
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_exp(tree.GetParam(1).GetImmed())));
                    return true;
                case cSinh:
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_asinh(tree.GetParam(1).GetImmed())));
                    return true;
                case cTanh:
                    if(fp_less(fp_abs(tree.GetParam(1).GetImmed()), Value_t(1)))
                    {
                        tree.SetParam(0, tree.GetParam(0).GetParam(0));
                        tree.SetParam(1, CodeTreeImmed(fp_atanh(tree.GetParam(1).GetImmed())));
                        return true;
                    }
                    break;
                default: break;
            }

        return false;
    }
}
