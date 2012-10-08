#include <algorithm>

#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include "codetree.hh"
#include "optimize.hh"
#include "consts.hh"

#include <assert.h>

#include "rangeestimation.hh"
#include "constantfolding.hh"


#include "logic_boolgroups.hh"
/* ^For ConstantFolding_AndLogic()
 *      ConstantFolding_OrLogic()
 *      ConstantFolding_MulLogicItems()
 *      ConstantFolding_AddLogicItems()
 */

#include "logic_collections.hh"
/* ^For ConstantFolding_MulGrouping()
 *      ConstantFolding_AddGrouping()
 */

#include "logic_ifoperations.hh"
/* ^For ConstantFolding_IfOperations()
 */

#include "logic_powoperations.hh"
/* ^For ConstantFolding_PowOperations()
 */

#include "logic_comparisons.hh"
/* ^For ConstantFolding_Comparison()
 */

#define FP_MUL_COMBINE_EXPONENTS

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /*************************************/
    /* ADOPTING SAME-TYPE CHILDREN       */
    /*************************************/

    template<typename Value_t>
    static void AdoptChildrenWithSameOpcode(CodeTree<Value_t>& tree)
    {
        /* If the list contains another list of the same kind, assimilate it */
      #ifdef DEBUG_SUBSTITUTIONS
        bool assimilated = false;
      #endif
        for(size_t a=tree.GetParamCount(); a-- > 0; )
            if(tree.GetParam(a).GetOpcode() == tree.GetOpcode())
            {
              #ifdef DEBUG_SUBSTITUTIONS
                if(!assimilated)
                {
                    std::cout << "Before assimilation: "; DumpTree(tree);
                    std::cout << "\n";
                    assimilated = true;
                }
              #endif
                // Assimilate its children and remove it
                tree.AddParamsMove(tree.GetParam(a).GetUniqueRef().GetParams(), a);
            }
      #ifdef DEBUG_SUBSTITUTIONS
        if(assimilated)
        {
            std::cout << "After assimilation:   "; DumpTree(tree);
            std::cout << "\n";
        }
      #endif
    }
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void ConstantFolding(CodeTree<Value_t>& tree)
    {
        tree.Sort(); // Otherwise "0 <= acos(x)" does not get properly optimized
    #ifdef DEBUG_SUBSTITUTIONS
        void* stackptr=0;
        std::cout << "[" << (&stackptr) << "]Runs ConstantFolding for: ";
        DumpTree(tree);
        std::cout << "\n";
        DumpHashes(tree); std::cout << std::flush;
    #endif
        if(false)
        {
    redo:;
            tree.Sort();
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "[" << (&stackptr) << "]Re-runs ConstantFolding: ";
            DumpTree(tree);
            std::cout << "\n";
            DumpHashes(tree);
        #endif
        }

        // Insert here any hardcoded constant-folding optimizations
        // that you want to be done whenever a new subtree is generated.
        /* Not recursive. */

        if(tree.GetOpcode() != cImmed)
        {
            range<Value_t> p = CalculateResultBoundaries(tree);
            if(p.min.known && p.max.known && p.min.val == p.max.val)
            {
                // Replace us with this immed
                tree.ReplaceWithImmed(p.min.val);
                goto do_return;
            }
        }

        if(false)
        {
            ReplaceTreeWithOne:  tree.ReplaceWithImmed(Value_t(1)); goto do_return;
            ReplaceTreeWithZero: tree.ReplaceWithImmed(Value_t(0)); goto do_return;
            ReplaceTreeWithParam0:
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before replace: ";
                std::cout << std::hex
                          << '[' << tree.GetHash().hash1
                          << ',' << tree.GetHash().hash2
                          << ']' << std::dec;
                DumpTree(tree);
                std::cout << "\n";
              #endif
                tree.Become(tree.GetParam(0));
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After replace: ";
                std::cout << std::hex
                          << '[' << tree.GetHash().hash1
                          << ',' << tree.GetHash().hash2
                          << ']' << std::dec;
                DumpTree(tree);
                std::cout << "\n";
              #endif
                goto redo;
        }

        /* Constant folding */
        switch(tree.GetOpcode())
        {
            case cImmed:
                break; // nothing to do
            case VarBegin:
                break; // nothing to do

            case cAnd:
            case cAbsAnd:
            {
                AdoptChildrenWithSameOpcode(tree);
                bool has_nonlogical_values = false;
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                {
                    if(!IsLogicalValue(tree.GetParam(a))) has_nonlogical_values = true;
                    switch(GetLogicalValue(tree.GetParam(a), tree.GetOpcode()==cAbsAnd))
                    {
                        case IsNever: goto ReplaceTreeWithZero;
                        case IsAlways: tree.DelParam(a); break; // x & y & 1 = x & y;  x & 1 = !!x
                        case Unknown: default: ;
                    }
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithOne;
                    case 1: tree.SetOpcode(tree.GetOpcode()==cAnd ? cNotNot : cAbsNotNot); goto redo; // Replace self with the single operand
                    default:
                        if(tree.GetOpcode()==cAnd || !has_nonlogical_values)
                            if(ConstantFolding_AndLogic(tree)) goto redo;
                }
                break;
            }
            case cOr:
            case cAbsOr:
            {
                AdoptChildrenWithSameOpcode(tree);
                bool has_nonlogical_values = false;
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                {
                    if(!IsLogicalValue(tree.GetParam(a))) has_nonlogical_values = true;
                    switch(GetLogicalValue(tree.GetParam(a), tree.GetOpcode()==cAbsOr))
                    {
                        case IsAlways: goto ReplaceTreeWithOne;
                        case IsNever: tree.DelParam(a); break;
                        case Unknown: default: ;
                    }
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: tree.SetOpcode(tree.GetOpcode()==cOr ? cNotNot : cAbsNotNot); goto redo; // Replace self with the single operand
                    default:
                        if(tree.GetOpcode()==cOr || !has_nonlogical_values)
                            if(ConstantFolding_OrLogic(tree)) goto redo;
                }
                break;
            }
            case cNot:
            case cAbsNot:
            {
                unsigned opposite = 0;
                switch(tree.GetParam(0).GetOpcode())
                {
                    case cEqual:       opposite = cNEqual; break;
                    case cNEqual:      opposite = cEqual; break;
                    case cLess:        opposite = cGreaterOrEq; break;
                    case cGreater:     opposite = cLessOrEq; break;
                    case cLessOrEq:    opposite = cGreater; break;
                    case cGreaterOrEq: opposite = cLess; break;
                    case cNotNot:      opposite = cNot; break;
                    case cNot:         opposite = cNotNot; break;
                    case cAbsNot:      opposite = cAbsNotNot; break;
                    case cAbsNotNot:   opposite = cAbsNot; break;
                    default: break;
                }
                if(opposite)
                {
                    tree.SetOpcode(OPCODE(opposite));
                    tree.SetParamsMove(tree.GetParam(0).GetUniqueRef().GetParams());
                    goto redo;
                }

                // If the sub-expression evaluates to approx. zero, yield one.
                // If the sub-expression evaluates to approx. nonzero, yield zero.
                switch(GetLogicalValue(tree.GetParam(0), tree.GetOpcode()==cAbsNot))
                {
                    case IsAlways: goto ReplaceTreeWithZero;
                    case IsNever: goto ReplaceTreeWithOne;
                    case Unknown: default: ;
                }
                if(tree.GetOpcode() == cNot && GetPositivityInfo(tree.GetParam(0)) == IsAlways)
                    tree.SetOpcode(cAbsNot);

                if(tree.GetParam(0).GetOpcode() == cIf
                || tree.GetParam(0).GetOpcode() == cAbsIf)
                {
                    CodeTree<Value_t> iftree = tree.GetParam(0);
                    const CodeTree<Value_t>& ifp1 = iftree.GetParam(1);
                    const CodeTree<Value_t>& ifp2 = iftree.GetParam(2);
                    if(ifp1.GetOpcode() == cNot
                    || ifp1.GetOpcode() == cAbsNot)
                    {
                        // cNot [(cIf [x (cNot[y]) z])] -> cIf [x (cNotNot[y]) (cNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree<Value_t> p1;
                        p1.SetOpcode(ifp1.GetOpcode()==cNot ? cNotNot : cAbsNotNot);
                        p1.AddParam(ifp1.GetParam(0));
                        p1.Rehash();
                        tree.AddParamMove(p1);
                        CodeTree<Value_t> p2;
                        p2.SetOpcode(tree.GetOpcode());
                        p2.AddParam(ifp2);
                        p2.Rehash();
                        tree.AddParamMove(p2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                    if(ifp2.GetOpcode() == cNot
                    || ifp2.GetOpcode() == cAbsNot)
                    {
                        // cNot [(cIf [x y (cNot[z])])] -> cIf [x (cNot[y]) (cNotNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree<Value_t> p1;
                        p1.SetOpcode(tree.GetOpcode());
                        p1.AddParam(ifp1);
                        p1.Rehash();
                        tree.AddParamMove(p1);
                        CodeTree<Value_t> p2;
                        p2.SetOpcode(ifp2.GetOpcode()==cNot ? cNotNot : cAbsNotNot);
                        p2.AddParam(ifp2.GetParam(0));
                        p2.Rehash();
                        tree.AddParamMove(p2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                }
                break;
            }
            case cNotNot:
            case cAbsNotNot:
            {
                // The function of cNotNot is to protect a logical value from
                // changing. If the parameter is already a logical value,
                // then the cNotNot opcode is redundant.
                if(IsLogicalValue(tree.GetParam(0)))
                    goto ReplaceTreeWithParam0;

                // If the sub-expression evaluates to approx. zero, yield zero.
                // If the sub-expression evaluates to approx. nonzero, yield one.
                switch(GetLogicalValue(tree.GetParam(0), tree.GetOpcode()==cAbsNotNot))
                {
                    case IsNever: goto ReplaceTreeWithZero;
                    case IsAlways: goto ReplaceTreeWithOne;
                    case Unknown: default: ;
                }
                if(tree.GetOpcode() == cNotNot && GetPositivityInfo(tree.GetParam(0)) == IsAlways)
                    tree.SetOpcode(cAbsNotNot);

                if(tree.GetParam(0).GetOpcode() == cIf
                || tree.GetParam(0).GetOpcode() == cAbsIf)
                {
                    CodeTree<Value_t> iftree = tree.GetParam(0);
                    const CodeTree<Value_t>& ifp1 = iftree.GetParam(1);
                    const CodeTree<Value_t>& ifp2 = iftree.GetParam(2);
                    if(ifp1.GetOpcode() == cNot
                    || ifp1.GetOpcode() == cAbsNot)
                    {
                        // cNotNot [(cIf [x (cNot[y]) z])] -> cIf [x (cNot[y]) (cNotNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        tree.AddParam(ifp1);
                        CodeTree<Value_t> p2;
                        p2.SetOpcode(tree.GetOpcode());
                        p2.AddParam(ifp2);
                        p2.Rehash();
                        tree.AddParamMove(p2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                    if(ifp2.GetOpcode() == cNot
                    || ifp2.GetOpcode() == cAbsNot)
                    {
                        // cNotNot [(cIf [x y (cNot[z])])] -> cIf [x (cNotNot[y]) (cNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree<Value_t> p1;
                        p1.SetOpcode(tree.GetOpcode());
                        p1.AddParam(ifp1);
                        p1.Rehash();
                        tree.AddParamMove(p1);
                        tree.AddParam(ifp2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                }
                break;
            }
            case cIf:
            case cAbsIf:
            {
                if(ConstantFolding_IfOperations(tree))
                    goto redo;
                break;
            }
            case cMul:
            {
            NowWeAreMulGroup: ;
                AdoptChildrenWithSameOpcode(tree);
                // If one sub-expression evalutes to exact zero, yield zero.
                Value_t immed_product = Value_t(1);
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    if(!tree.GetParam(a).IsImmed()) continue;
                    // ^ Only check constant values
                    Value_t immed = tree.GetParam(a).GetImmed();
                    if(immed == Value_t(0) ) goto ReplaceTreeWithZero;
                    immed_product *= immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1 || (n_immeds == 1 && fp_equal(immed_product, Value_t(1))))
                    needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "cMul: Will add new immed " << immed_product << "\n";
                #endif
                    for(size_t a=tree.GetParamCount(); a-->0; )
                        if(tree.GetParam(a).IsImmed())
                        {
                        #ifdef DEBUG_SUBSTITUTIONS
                            std::cout << " - For that, deleting immed " << tree.GetParam(a).GetImmed();
                            std::cout << "\n";
                        #endif
                            tree.DelParam(a);
                        }
                    if(!fp_equal(immed_product, Value_t(1)))
                        tree.AddParam( CodeTreeImmed<Value_t> (immed_product) );
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithOne;
                    case 1: goto ReplaceTreeWithParam0; // Replace self with the single operand
                    default:
                        if(ConstantFolding_MulGrouping(tree)) goto redo;
                        if(ConstantFolding_MulLogicItems(tree)) goto redo;
                }
                break;
            }
            case cAdd:
            {
                AdoptChildrenWithSameOpcode(tree);
                Value_t immed_sum = 0.0;
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    if(!tree.GetParam(a).IsImmed()) continue;
                    // ^ Only check constant values
                    Value_t immed = tree.GetParam(a).GetImmed();
                    immed_sum += immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1 || (n_immeds == 1 && immed_sum == Value_t(0)))
                    needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "cAdd: Will add new immed " << immed_sum << "\n";
                    std::cout << "In: "; DumpTree(tree);
                    std::cout << "\n";
                #endif
                    for(size_t a=tree.GetParamCount(); a-->0; )
                        if(tree.GetParam(a).IsImmed())
                        {
                        #ifdef DEBUG_SUBSTITUTIONS
                            std::cout << " - For that, deleting immed " << tree.GetParam(a).GetImmed();
                            std::cout << "\n";
                        #endif
                            tree.DelParam(a);
                        }
                    if(!(immed_sum == Value_t(0.0)))
                        tree.AddParam( CodeTreeImmed<Value_t> (immed_sum) );
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: goto ReplaceTreeWithParam0; // Replace self with the single operand
                    default:
                        if(ConstantFolding_AddGrouping(tree)) goto redo;
                        if(ConstantFolding_AddLogicItems(tree)) goto redo;
                }
                break;
            }
            case cMin:
            {
                AdoptChildrenWithSameOpcode(tree);
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the larger one.
                 *
                 * Algorithm: 1. figure out the smallest maximum of all operands.
                 *            2. eliminate all operands where their minimum is
                 *               larger than the selected maximum.
                 */
                size_t preserve=0;
                range<Value_t> smallest_maximum;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    while(a+1 < tree.GetParamCount() && tree.GetParam(a).IsIdenticalTo(tree.GetParam(a+1)))
                        tree.DelParam(a+1);
                    range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                    if(p.max.known && (!smallest_maximum.max.known || (p.max.val) < smallest_maximum.max.val))
                    {
                        smallest_maximum.max.val = p.max.val;
                        smallest_maximum.max.known = true;
                        preserve=a;
                }   }
                if(smallest_maximum.max.known)
                    for(size_t a=tree.GetParamCount(); a-- > 0; )
                    {
                        range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                        if(p.min.known && a != preserve && p.min.val >= smallest_maximum.max.val)
                            tree.DelParam(a);
                    }
                //fprintf(stderr, "Remains: %u\n", (unsigned)tree.GetParamCount());
                if(tree.GetParamCount() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }
            case cMax:
            {
                AdoptChildrenWithSameOpcode(tree);
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the smaller one.
                 *
                 * Algorithm: 1. figure out the biggest minimum of all operands.
                 *            2. eliminate all operands where their maximum is
                 *               smaller than the selected minimum.
                 */
                size_t preserve=0;
                range<Value_t> biggest_minimum;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    while(a+1 < tree.GetParamCount() && tree.GetParam(a).IsIdenticalTo(tree.GetParam(a+1)))
                        tree.DelParam(a+1);
                    range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                    if(p.min.known && (!biggest_minimum.min.known || p.min.val > biggest_minimum.min.val))
                    {
                        biggest_minimum.min.val = p.min.val;
                        biggest_minimum.min.known = true;
                        preserve=a;
                }   }
                if(biggest_minimum.min.known)
                {
                    //fprintf(stderr, "Removing all where max < %g\n", biggest_minimum.min.val);
                    for(size_t a=tree.GetParamCount(); a-- > 0; )
                    {
                        range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                        if(p.max.known && a != preserve && (p.max.val) < biggest_minimum.min.val)
                        {
                            //fprintf(stderr, "Removing %g\n", p.max.val);
                            tree.DelParam(a);
                        }
                    }
                }
                //fprintf(stderr, "Remains: %u\n", (unsigned)tree.GetParamCount());
                if(tree.GetParamCount() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }

            case cEqual:
            case cNEqual:
            case cLess:
            case cGreater:
            case cLessOrEq:
            case cGreaterOrEq:
                if(ConstantFolding_Comparison(tree)) goto redo;
                break;

            case cAbs:
            {
                /* If we know the operand is always positive, cAbs is redundant.
                 * If we know the operand is always negative, use actual negation.
                 */
                range<Value_t> p0 = CalculateResultBoundaries( tree.GetParam(0) );
                if(p0.min.known && p0.min.val >= Value_t(0.0))
                    goto ReplaceTreeWithParam0;
                if(p0.max.known && p0.max.val <= fp_const_negativezero<Value_t>())
                {
                    /* abs(negative) = negative*-1 */
                    tree.SetOpcode(cMul);
                    tree.AddParam( CodeTreeImmed(Value_t(1)) );
                    /* The caller of ConstantFolding() will do Sort() and Rehash() next.
                     * Thus, no need to do it here. */
                    /* We were changed into a cMul group. Do cMul folding. */
                    goto NowWeAreMulGroup;
                }
                /* If the operand is a cMul group, find elements
                 * that are always positive and always negative,
                 * and move them out, e.g. abs(p*n*x*y) = p*(-n)*abs(x*y)
                 */
                if(tree.GetParam(0).GetOpcode() == cMul)
                {
                    const CodeTree<Value_t>& p = tree.GetParam(0);
                    std::vector<CodeTree<Value_t> > pos_set;
                    std::vector<CodeTree<Value_t> > neg_set;
                    for(size_t a=0; a<p.GetParamCount(); ++a)
                    {
                        p0 = CalculateResultBoundaries( p.GetParam(a) );
                        if(p0.min.known && p0.min.val >= Value_t(0.0))
                            { pos_set.push_back(p.GetParam(a)); }
                        if(p0.max.known && p0.max.val <= fp_const_negativezero<Value_t>())
                            { neg_set.push_back(p.GetParam(a)); }
                    }
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "Abs: mul group has " << pos_set.size()
                              << " pos, " << neg_set.size() << "neg\n";
                #endif
                    if(!pos_set.empty() || !neg_set.empty())
                    {
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-Before: ";
                        DumpTree(tree);
                        std::cout << "\n" << std::flush;
                        DumpHashes(tree, std::cout);
                #endif
                        CodeTree<Value_t> pclone;
                        pclone.SetOpcode(cMul);
                        for(size_t a=0; a<p.GetParamCount(); ++a)
                        {
                            p0 = CalculateResultBoundaries( p.GetParam(a) );
                            if((p0.min.known && p0.min.val >= Value_t(0.0))
                            || (p0.max.known && p0.max.val <= fp_const_negativezero<Value_t>()))
                                {/*pclone.DelParam(a);*/}
                            else
                                pclone.AddParam( p.GetParam(a) );
                            /* Here, p*n*x*y -> x*y.
                             * p is saved in pos_set[]
                             * n is saved in neg_set[]
                             */
                        }
                        pclone.Rehash();
                        CodeTree<Value_t> abs_mul;
                        abs_mul.SetOpcode(cAbs);
                        abs_mul.AddParamMove(pclone);
                        abs_mul.Rehash();
                        CodeTree<Value_t> mulgroup;
                        mulgroup.SetOpcode(cMul);
                        mulgroup.AddParamMove(abs_mul); // cAbs[whatever remains in p]
                        mulgroup.AddParamsMove(pos_set);
                        /* Now:
                         * mulgroup  = p * Abs(x*y)
                         */
                        if(!neg_set.empty())
                        {
                            if(neg_set.size() % 2)
                                mulgroup.AddParam( CodeTreeImmed(Value_t(-1)) );
                            mulgroup.AddParamsMove(neg_set);
                            /* Now:
                             * mulgroup = p * n * -1 * Abs(x*y)
                             */
                        }
                        tree.Become(mulgroup);
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-After: ";
                        DumpTree(tree, std::cout);
                        std::cout << "\n" << std::flush;
                        DumpHashes(tree, std::cout);
                #endif
                        /* We were changed into a cMul group. Do cMul folding. */
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            }

            #define HANDLE_UNARY_CONST_FUNC(funcname) \
                if(tree.GetParam(0).IsImmed()) \
                    { tree.ReplaceWithImmed( funcname(tree.GetParam(0).GetImmed()) ); \
                      goto do_return; }

            case cLog:
                HANDLE_UNARY_CONST_FUNC(fp_log);
                if(tree.GetParam(0).GetOpcode() == cPow)
                {
                    CodeTree<Value_t> pow = tree.GetParam(0);
                    if(GetPositivityInfo(pow.GetParam(0)) == IsAlways)  // log(posi ^ y) = y*log(posi)
                    {
                        pow.CopyOnWrite();
                        pow.SetOpcode(cLog);
                        tree.SetOpcode(cMul);
                        tree.AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        tree.SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                    if(GetEvennessInfo(pow.GetParam(1)) == IsAlways) // log(x ^ even) = even*log(abs(x))
                    {
                        pow.CopyOnWrite();
                        CodeTree<Value_t> abs;
                        abs.SetOpcode(cAbs);
                        abs.AddParamMove(pow.GetParam(0));
                        abs.Rehash();
                        pow.SetOpcode(cLog);
                        tree.SetOpcode(cMul);
                        pow.SetParamMove(0, abs);
                        tree.AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        tree.SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                }
                else if(tree.GetParam(0).GetOpcode() == cAbs)
                {
                    // log(abs(x^y)) = y*log(abs(x))
                    CodeTree<Value_t> pow = tree.GetParam(0).GetParam(0);
                    if(pow.GetOpcode() == cPow)
                    {
                        pow.CopyOnWrite();
                        CodeTree<Value_t> abs;
                        abs.SetOpcode(cAbs);
                        abs.AddParamMove(pow.GetParam(0));
                        abs.Rehash();
                        pow.SetOpcode(cLog);
                        tree.SetOpcode(cMul);
                        pow.SetParamMove(0, abs);
                        tree.AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        tree.SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            case cAcosh: HANDLE_UNARY_CONST_FUNC(fp_acosh); break;
            case cAsinh: HANDLE_UNARY_CONST_FUNC(fp_asinh); break;
            case cAtanh: HANDLE_UNARY_CONST_FUNC(fp_atanh); break;
            case cAcos: HANDLE_UNARY_CONST_FUNC(fp_acos); break;
            case cAsin: HANDLE_UNARY_CONST_FUNC(fp_asin); break;
            case cAtan: HANDLE_UNARY_CONST_FUNC(fp_atan); break;
            case cCosh: HANDLE_UNARY_CONST_FUNC(fp_cosh); break;
            case cSinh: HANDLE_UNARY_CONST_FUNC(fp_sinh); break;
            case cTanh: HANDLE_UNARY_CONST_FUNC(fp_tanh); break;
            case cSin: HANDLE_UNARY_CONST_FUNC(fp_sin); break;
            case cCos: HANDLE_UNARY_CONST_FUNC(fp_cos); break;
            case cTan: HANDLE_UNARY_CONST_FUNC(fp_tan); break;
            case cCeil:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_ceil); break;
            case cTrunc:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_trunc); break;
            case cFloor:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_floor); break;
            case cInt:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_int); break;
            case cCbrt: HANDLE_UNARY_CONST_FUNC(fp_cbrt); break; // converted into cPow x 0.33333
            case cSqrt: HANDLE_UNARY_CONST_FUNC(fp_sqrt); break; // converted into cPow x 0.5
            case cExp: HANDLE_UNARY_CONST_FUNC(fp_exp); break; // convered into cPow CONSTANT_E x
            case cLog2: HANDLE_UNARY_CONST_FUNC(fp_log2); break;
            case cLog10: HANDLE_UNARY_CONST_FUNC(fp_log10); break;

            case cLog2by:
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_log2(tree.GetParam(0).GetImmed()) * tree.GetParam(1).GetImmed() );
                      goto do_return; }
                break;

            case cArg: HANDLE_UNARY_CONST_FUNC(fp_arg); break;
            case cConj: HANDLE_UNARY_CONST_FUNC(fp_conj); break;
            case cImag: HANDLE_UNARY_CONST_FUNC(fp_imag); break;
            case cReal: HANDLE_UNARY_CONST_FUNC(fp_real); break;
            case cPolar:
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_polar(tree.GetParam(0).GetImmed(), tree.GetParam(1).GetImmed() ) );
                      goto do_return; }
                break;

            case cMod: /* Can more be done than this? */
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_mod(tree.GetParam(0).GetImmed(), tree.GetParam(1).GetImmed()) );
                      goto do_return; }
                break;

            case cAtan2:
            {
                /* Range based optimizations for (y,x):
                 * If y is +0 and x <= -0, +pi is returned
                 * If y is -0 and x <= -0, -pi is returned (assumed never happening)
                 * If y is +0 and x >= +0, +0 is returned
                 * If y is -0 and x >= +0, -0 is returned  (assumed never happening)
                 * If x is +-0 and y < 0, -pi/2 is returned
                 * If x is +-0 and y > 0, +pi/2 is returned
                 * Otherwise, perform constant folding when available
                 * If we know x <> 0, convert into atan(y / x)
                 *   TODO: Figure out whether the above step is wise
                 *         It allows e.g. atan2(6*x, 3*y) -> atan(2*x/y)
                 *         when we know y != 0
                 */
                range<Value_t> p0 = CalculateResultBoundaries( tree.GetParam(0) );
                range<Value_t> p1 = CalculateResultBoundaries( tree.GetParam(1) );
                if(tree.GetParam(0).IsImmed()
                && fp_equal(tree.GetParam(0).GetImmed(), Value_t(0)))   // y == 0
                {
                    if(p1.max.known && (p1.max.val) < Value_t(0))    // y == 0 && x < 0
                        { tree.ReplaceWithImmed( fp_const_pi<Value_t>() ); goto do_return; }
                    if(p1.min.known && p1.min.val >= Value_t(0.0))  // y == 0 && x >= 0.0
                        { tree.ReplaceWithImmed( Value_t(0) ); goto do_return; }
                }
                if(tree.GetParam(1).IsImmed()
                && fp_equal(tree.GetParam(1).GetImmed(), Value_t(0)))   // x == 0
                {
                    if(p0.max.known && (p0.max.val) < Value_t(0))   // y < 0 && x == 0
                        { tree.ReplaceWithImmed( -fp_const_pihalf<Value_t>() ); goto do_return; }
                    if(p0.min.known && p0.min.val > Value_t(0))     // y > 0 && x == 0
                        { tree.ReplaceWithImmed(  fp_const_pihalf<Value_t>() ); goto do_return; }
                }
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_atan2(tree.GetParam(0).GetImmed(),
                                             tree.GetParam(1).GetImmed()) );
                      goto do_return; }
                if((p1.min.known && p1.min.val > Value_t(0))            // p1 != 0.0
                || (p1.max.known && (p1.max.val) < fp_const_negativezero<Value_t>())) // become atan(p0 / p1)
                {
                    CodeTree<Value_t> pow_tree;
                    pow_tree.SetOpcode(cPow);
                    pow_tree.AddParamMove(tree.GetParam(1));
                    pow_tree.AddParam(CodeTreeImmed(Value_t(-1)));
                    pow_tree.Rehash();
                    CodeTree<Value_t> div_tree;
                    div_tree.SetOpcode(cMul);
                    div_tree.AddParamMove(tree.GetParam(0));
                    div_tree.AddParamMove(pow_tree);
                    div_tree.Rehash();
                    tree.SetOpcode(cAtan);
                    tree.SetParamMove(0, div_tree);
                    tree.DelParam(1);
                }
                break;
            }

            case cPow:
            {
                if(ConstantFolding_PowOperations(tree)) goto redo;
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.cc and thus
             * they will never occur in the calling context for the
             * most of the parsing context. They may however occur
             * at the late phase, so we deal with them.
             */
            case cDiv: // converted into cPow y -1
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed()
                && tree.GetParam(1).GetImmed() != Value_t(0.0))
                    { tree.ReplaceWithImmed( tree.GetParam(0).GetImmed() / tree.GetParam(1).GetImmed() );
                      goto do_return; }
                break;
            case cInv: // converted into cPow y -1
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(0).GetImmed() != Value_t(0.0))
                    { tree.ReplaceWithImmed( Value_t(1) / tree.GetParam(0).GetImmed() );
                      goto do_return; }
                // Note: Could use (mulgroup)^immed optimization from cPow
                break;
            case cSub: // converted into cMul y -1
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( tree.GetParam(0).GetImmed() - tree.GetParam(1).GetImmed() );
                      goto do_return; }
                break;
            case cNeg: // converted into cMul x -1
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( -tree.GetParam(0).GetImmed() );
                      goto do_return; }
                break;
            case cRad: // converted into cMul x CONSTANT_RD
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( RadiansToDegrees( tree.GetParam(0).GetImmed() ) );
                      goto do_return; }
                break;
            case cDeg: // converted into cMul x CONSTANT_DR
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( DegreesToRadians( tree.GetParam(0).GetImmed() ) );
                      goto do_return; }
                break;
            case cSqr: // converted into cMul x x
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( tree.GetParam(0).GetImmed() * tree.GetParam(0).GetImmed() );
                      goto do_return; }
                break;
            case cExp2: // converted into cPow 2.0 x
                HANDLE_UNARY_CONST_FUNC(fp_exp2); break;
            case cRSqrt: // converted into cPow x -0.5
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( Value_t(1) / fp_sqrt(tree.GetParam(0).GetImmed()) );
                      goto do_return; }
                break;
            case cCot: // converted into cMul (cPow (cTan x) -1)
                if(tree.GetParam(0).IsImmed())
                    { Value_t tmp = fp_tan(tree.GetParam(0).GetImmed());
                      if(fp_nequal(tmp, Value_t(0)))
                      { tree.ReplaceWithImmed( Value_t(1) / tmp );
                        goto do_return; } }
                break;
            case cSec: // converted into cMul (cPow (cCos x) -1)
                if(tree.GetParam(0).IsImmed())
                    { Value_t tmp = fp_cos(tree.GetParam(0).GetImmed());
                      if(fp_nequal(tmp, Value_t(0)))
                      { tree.ReplaceWithImmed( Value_t(1) / tmp );
                        goto do_return; } }
                break;
            case cCsc: // converted into cMul (cPow (cSin x) -1)
                if(tree.GetParam(0).IsImmed())
                    { Value_t tmp = fp_sin(tree.GetParam(0).GetImmed());
                      if(fp_nequal(tmp, Value_t(0)))
                      { tree.ReplaceWithImmed( Value_t(1) / tmp );
                        goto do_return; } }
                break;
            case cHypot: // converted into cSqrt(cAdd(cMul(x x), cMul(y y)))
                if(tree.GetParam(0).IsImmed() && tree.GetParam(1).IsImmed())
                {
                    tree.ReplaceWithImmed( fp_hypot(tree.GetParam(0).GetImmed(),
                                           tree.GetParam(1).GetImmed()) );
                    goto do_return;
                }
                break;

            /* Opcodes that do not occur in the tree for other reasons */
            case cRDiv: // version of cDiv
            case cRSub: // version of cSub
            case cDup:
            case cFetch:
            case cPopNMov:
            case cSinCos:
            case cSinhCosh:
            case cNop:
            case cJump:
                break; /* Should never occur */

            /* Opcodes that we can't do anything about */
            case cPCall:
            case cFCall:
                break;
        }
    do_return:;
#ifdef DEBUG_SUBSTITUTIONS
        std::cout << "[" << (&stackptr) << "]Done ConstantFolding, result: ";
        DumpTree(tree);
        std::cout << "\n";
        DumpHashes(tree);
#endif
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template void ConstantFolding(CodeTree<type>& );
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
