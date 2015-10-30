#include "codetree.hh"

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /***************************/
    /* LOGIC (AND, OR, NOT)    */
    /***************************/

    struct ComparisonSetBase
    {
        enum { Lt_Mask = 0x1,   // 1=less
               Eq_Mask = 0x2,   // 2=equal
               Le_Mask = 0x3,   // 1+2 = Less or Equal
               Gt_Mask = 0x4,   // 4=greater
               Ne_Mask = 0x5,   // 4+1 = Greater or Less, i.e. Not equal
               Ge_Mask = 0x6 }; // 4+2 = Greater or Equal
        static int Swap_Mask(int m) { return (m&Eq_Mask)
                                  | ((m&Lt_Mask) ? Gt_Mask : 0)
                                  | ((m&Gt_Mask) ? Lt_Mask : 0); }
        enum RelationshipResult
        {
            Ok,
            BecomeZero,
            BecomeOne,
            Suboptimal
        };
        enum ConditionType
        {
            cond_or,
            cond_and,
            cond_mul,
            cond_add
        };
    };

    template<typename Value_t>
    struct ComparisonSet: public ComparisonSetBase /* For optimizing And, Or */
    {
        struct Comparison
        {
            CodeTree<Value_t> a;
            CodeTree<Value_t> b;
            int relationship;

            Comparison() : a(),b(), relationship() {}
        };
        std::vector<Comparison> relationships;
        struct Item
        {
            CodeTree<Value_t> value;
            bool negated;

            Item() : value(), negated(false) {}
        };
        std::vector<Item> plain_set;
        int const_offset;

        ComparisonSet():
            relationships(),
            plain_set(),
            const_offset(0)
        {
        }

        RelationshipResult AddItem(
            const CodeTree<Value_t>& a,
            bool negated,
            ConditionType type)
        {
            for(size_t c=0; c<plain_set.size(); ++c)
                if(plain_set[c].value.IsIdenticalTo(a))
                {
                    if(negated != plain_set[c].negated)
                    {
                        switch(type)
                        {
                            case cond_or:
                                return BecomeOne;
                            case cond_add:
                                plain_set.erase(plain_set.begin() + c);
                                const_offset += 1;
                                return Suboptimal;
                            case cond_and:
                            case cond_mul:
                                return BecomeZero;
                        }
                    }
                    return Suboptimal;
                }
            Item pole;
            pole.value   = a;
            pole.negated = negated;
            plain_set.push_back(pole);
            return Ok;
        }

        /* Note: Trees are passed by-value so we can use swap() on them safely. */
        RelationshipResult AddRelationship
            (CodeTree<Value_t> a,
             CodeTree<Value_t> b,
             int reltype,
             ConditionType type)
        {
            switch(type)
            {
                case cond_or:
                    if(reltype == 7) return BecomeOne;
                    break;
                case cond_add:
                    if(reltype == 7) { const_offset += 1; return Suboptimal; }
                    break;
                case cond_and:
                case cond_mul:
                    if(reltype == 0) return BecomeZero;
                    break;
            }

            if(!(a.GetHash() < b.GetHash()))
            {
                a.swap(b);
                reltype = Swap_Mask(reltype);
            }

            for(size_t c=0; c<relationships.size(); ++c)
            {
                if(relationships[c].a.IsIdenticalTo(a)
                && relationships[c].b.IsIdenticalTo(b))
                {
                    switch(type)
                    {
                        case cond_or:
                        {
                            int newrel = relationships[c].relationship | reltype;
                            if(newrel == 7) return BecomeOne;
                            relationships[c].relationship = newrel;
                            break;
                        }
                        case cond_and:
                        case cond_mul:
                        {
                            int newrel = relationships[c].relationship & reltype;
                            if(newrel == 0) return BecomeZero;
                            relationships[c].relationship = newrel;
                            break;
                        }
                        case cond_add:
                        {
                            int newrel_or  = relationships[c].relationship | reltype;
                            int newrel_and = relationships[c].relationship & reltype;
                            if(newrel_or  == 5 // < + >
                            && newrel_and == 0)
                            {
                                // (x<y) + (x>y) = x!=y
                                relationships[c].relationship = Ne_Mask;
                                return Suboptimal;
                            }
                            if(newrel_or  == 7
                            && newrel_and == 0)
                            {
                                // (x<y) + (x>=y) = 1
                                // (x<=y) + (x>y) = 1
                                // (x=y) + (x!=y) = 1
                                const_offset += 1;
                                relationships.erase(relationships.begin()+c);
                                return Suboptimal;
                            }
                            if(newrel_or  == 7
                            && newrel_and == Eq_Mask)
                            {
                                // (x<=y) + (x>=y) = 1 + (x=y)
                                relationships[c].relationship = Eq_Mask;
                                const_offset += 1;
                                return Suboptimal;
                            }
                            continue;
                        }
                    }
                    return Suboptimal;
                }
            }
            Comparison comp;
            comp.a = a;
            comp.b = b;
            comp.relationship = reltype;
            relationships.push_back(comp);
            return Ok;
        }
    };

    template<typename Value_t, typename CondType> /* ComparisonSet::ConditionType */
    bool ConstantFolding_LogicCommon(
        CodeTree<Value_t>& tree, CondType cond_type, bool is_logical)
    {
        bool should_regenerate = false;
        ComparisonSet<Value_t> comp;
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            typename ComparisonSetBase::RelationshipResult
                change = ComparisonSetBase::Ok;
            const CodeTree<Value_t>& atree = tree.GetParam(a);
            switch(atree.GetOpcode())
            {
                case cEqual:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSetBase::Eq_Mask, cond_type);
                    break;
                case cNEqual:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSetBase::Ne_Mask, cond_type);
                    break;
                case cLess:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSetBase::Lt_Mask, cond_type);
                    break;
                case cLessOrEq:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSetBase::Le_Mask, cond_type);
                    break;
                case cGreater:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSetBase::Gt_Mask, cond_type);
                    break;
                case cGreaterOrEq:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSetBase::Ge_Mask, cond_type);
                    break;
                case cNot:
                    change = comp.AddItem(atree.GetParam(0), true, cond_type);
                    break;
                case cNotNot:
                    change = comp.AddItem(atree.GetParam(0), false, cond_type);
                    break;
                default:
                    if(is_logical || IsLogicalValue(atree))
                        change = comp.AddItem(atree, false, cond_type);
            }
            switch(change)
            {
            ReplaceTreeWithZero:
                    tree.ReplaceWithImmed(0);
                    return true;
            ReplaceTreeWithOne:
                    tree.ReplaceWithImmed(1);
                    return true;
                case ComparisonSetBase::Ok: // ok
                    break;
                case ComparisonSetBase::BecomeZero: // whole set was invalidated
                    goto ReplaceTreeWithZero;
                case ComparisonSetBase::BecomeOne: // whole set was validated
                    goto ReplaceTreeWithOne;
                case ComparisonSetBase::Suboptimal: // something was changed
                    should_regenerate = true;
                    break;
            }
        }
        if(should_regenerate)
        {
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Before ConstantFolding_LogicCommon: "; DumpTree(tree);
            std::cout << "\n";
          #endif

            if(is_logical)
            {
                tree.DelParams(); // delete all params
            }
            else
            {
                // Delete only logical params
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                {
                    const CodeTree<Value_t>& atree = tree.GetParam(a);
                    if(IsLogicalValue(atree))
                        tree.DelParam(a);
                }
            }

            for(size_t a=0; a<comp.plain_set.size(); ++a)
            {
                if(comp.plain_set[a].negated)
                {
                    CodeTree<Value_t> r;
                    r.SetOpcode(cNot);
                    r.AddParamMove(comp.plain_set[a].value);
                    r.Rehash();
                    tree.AddParamMove(r);
                }
                else if(!is_logical)
                {
                    CodeTree<Value_t> r;
                    r.SetOpcode(cNotNot);
                    r.AddParamMove(comp.plain_set[a].value);
                    r.Rehash();
                    tree.AddParamMove(r);
                }
                else
                    tree.AddParamMove(comp.plain_set[a].value);
            }
            for(size_t a=0; a<comp.relationships.size(); ++a)
            {
                CodeTree<Value_t> r;
                r.SetOpcode(cNop); // dummy
                switch(comp.relationships[a].relationship)
                {
                    case ComparisonSetBase::Lt_Mask: r.SetOpcode( cLess ); break;
                    case ComparisonSetBase::Eq_Mask: r.SetOpcode( cEqual ); break;
                    case ComparisonSetBase::Gt_Mask: r.SetOpcode( cGreater ); break;
                    case ComparisonSetBase::Le_Mask: r.SetOpcode( cLessOrEq ); break;
                    case ComparisonSetBase::Ne_Mask: r.SetOpcode( cNEqual ); break;
                    case ComparisonSetBase::Ge_Mask: r.SetOpcode( cGreaterOrEq ); break;
                }
                r.AddParamMove(comp.relationships[a].a);
                r.AddParamMove(comp.relationships[a].b);
                r.Rehash();
                tree.AddParamMove(r);
            }
            if(comp.const_offset != 0)
                tree.AddParam( CodeTreeImmed( Value_t(comp.const_offset) ) );
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "After ConstantFolding_LogicCommon: "; DumpTree(tree);
            std::cout << "\n";
          #endif
            return true;
        }
        /*
        Note: One thing this does not yet do, is to detect chains
              such as x=y & y=z & x=z, which could be optimized
              to x=y & x=z.
        */
        return false;
    }

    /* ConstantFolding_AndLogic:
     *       (x > y) & (x >= y)   --> (x > y)
     *       (x <= y) & (x >= y)  --> (x = y)
     *       (x <= y) & (x > y)   --> 0
     *       !x & !!x             --> 0
     * etc.
     */
    template<typename Value_t>
    bool ConstantFolding_AndLogic(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cAnd || tree.GetOpcode() == cAbsAnd);
        return ConstantFolding_LogicCommon(tree, ComparisonSetBase::cond_and, true );
    }

    /* ConstantFolding_OrLogic:
     *       (x > y) | (x >= y)   --> (x >= y)
     *       (x <= y) | (x >= y)  --> 1
     *       (x <= y) | (x > y)   --> 1
     *       !x | !!x             --> 1
     * etc.
     */
    template<typename Value_t>
    bool ConstantFolding_OrLogic(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cOr || tree.GetOpcode() == cAbsOr);
        return ConstantFolding_LogicCommon(tree, ComparisonSetBase::cond_or, true );
    }

    /* ConstantFolding_AddLogic:
     *       (x <= y) + (x >= y)  --> (x = y) + 1
     *       !x + !!x             --> 1
     * etc.
     */
    template<typename Value_t>
    bool ConstantFolding_AddLogicItems(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cAdd);
        return ConstantFolding_LogicCommon(tree, ComparisonSetBase::cond_add, false );
    }

    /* ConstantFolding_MulLogic:
     *       (x <= y) * (x >= y)  --> (x = y)
     *       (x > y)  * (x >= y)  --> (x > y)
     *       !x * !!x             --> 0
     * etc.
     */
    template<typename Value_t>
    bool ConstantFolding_MulLogicItems(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cMul);
        return ConstantFolding_LogicCommon(tree, ComparisonSetBase::cond_mul, false );
    }
}
