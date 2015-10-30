#include <vector>
#include <map>
#include <algorithm>

#include "codetree.hh"
#include "../lib/functional.hh"

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /**************************************/
    /* GROUPING OF COMMON FACTORS / TERMS */
    /**************************************/
    struct CollectionSetBase
    {
        enum CollectionResult
        {
            Ok,
            Suboptimal
        };
    };

    template<typename Value_t>
    struct CollectionSet: public CollectionSetBase /* For optimizing Add,  Mul */
    {
        struct Collection
        {
            CodeTree<Value_t> value;
            CodeTree<Value_t> factor;
            bool factor_needs_rehashing;

            Collection() : value(),factor(), factor_needs_rehashing(false) { }
            Collection(const CodeTree<Value_t>& v,
                       const CodeTree<Value_t>& f)
                : value(v), factor(f), factor_needs_rehashing(false) { }
        };
        std::multimap<fphash_t, Collection> collections;

        typedef typename
            std::multimap<fphash_t, Collection>::iterator
            PositionType;

        CollectionSet() : collections() {}

        PositionType FindIdenticalValueTo(const CodeTree<Value_t>& value)
        {
            fphash_t hash = value.GetHash();
            for(PositionType
                i = collections.lower_bound(hash);
                i != collections.end() && i->first == hash;
                ++i)
            {
                if(value.IsIdenticalTo(i->second.value))
                    return i;
            }
            return collections.end();
        }
        bool Found(const PositionType& b) { return b != collections.end(); }

        CollectionResult AddCollectionTo
            (const CodeTree<Value_t>& factor,
             const PositionType& into_which)
        {
            Collection& c = into_which->second;
            if(c.factor_needs_rehashing)
                c.factor.AddParam(factor);
            else
            {
                CodeTree<Value_t> add;
                add.SetOpcode( cAdd );
                add.AddParamMove(c.factor);
                add.AddParam(factor);
                c.factor.swap(add);
                c.factor_needs_rehashing = true;
            }
            return Suboptimal;
        }

        CollectionResult AddCollection
            (const CodeTree<Value_t>& value,
             const CodeTree<Value_t>& factor)
        {
            const fphash_t hash = value.GetHash();
            PositionType i = collections.lower_bound(hash);
            for(; i != collections.end() && i->first == hash; ++i)
            {
                if(i->second.value.IsIdenticalTo(value))
                    return AddCollectionTo(factor, i);
            }
            collections.insert(
                i,
                std::make_pair( hash, Collection(value, factor) ) );
            return Ok;
        }

        CollectionResult AddCollection(const CodeTree<Value_t>& a)
        {
            return AddCollection(a, CodeTreeImmed(Value_t(1)) );
        }
    };

    template<typename Value_t>
    struct ConstantExponentCollection
    {
        typedef std::vector<CodeTree<Value_t> > TreeSet;
        typedef std::pair<Value_t, TreeSet> ExponentInfo;
        std::vector<ExponentInfo> data;

        ConstantExponentCollection(): data(){}

        void MoveToSet_Unique(const Value_t& exponent, TreeSet& source_set)
        {
            data.push_back( std::pair<Value_t, TreeSet >
                            (exponent, TreeSet() ) );
            data.back().second.swap(source_set);
        }
        void MoveToSet_NonUnique(const Value_t& exponent, TreeSet& source_set)
        {
            typename std::vector<ExponentInfo>::iterator i
                = std::lower_bound(data.begin(), data.end(), exponent, Compare1st());
            if(i != data.end() && i->first == exponent)
            {
                i->second.insert(i->second.end(), source_set.begin(), source_set.end());
            }
            else
            {
                //MoveToSet_Unique(exponent, source_set);
                data.insert(i,  std::pair<Value_t, TreeSet >
                                (exponent, source_set) );
            }
        }

        bool Optimize()
        {
            /* TODO: Group them such that:
             *
             *      x^3 *         z^2 becomes (x*z)^2 * x^1
             *      x^3 * y^2.5 * z^2 becomes (x*z*y)^2 * y^0.5 * x^1
             *                    rather than (x*y*z)^2 * (x*y)^0.5 * x^0.5
             *
             *      x^4.5 * z^2.5     becomes (z * x)^2.5 * x^2
             *                        becomes (x*z*x)^2 * (z*x)^0.5
             *                        becomes (z*x*x*z*x)^0.5 * (z*x*x)^1.5 -- buzz, bad.
             *
             */
            bool changed = false;
            std::sort( data.begin(), data.end(), Compare1st() );
        redo:
            /* Supposed algorithm:
             * For the smallest pair of data[] where the difference
             * between the two is a "neat value" (x*16 is positive integer),
             * do the combining as indicated above.
             */
            /*
             * NOTE: Hanged in Testbed test P44, looping the following
             *       (Var0 ^ 0.75) * ((1.5 * Var0) ^ 1.0)
             *     = (Var0 ^ 1.75) *  (1.5         ^ 1.0)
             *       Fixed by limiting to cases where (exp_a != 1.0).
             *
             * NOTE: Converting (x*z)^0.5 * x^16.5
             *              into x^17 * z^0.5
             *       is handled by code within CollectMulGroup().
             *       However, bacause it is prone for infinite looping,
             *       the use of "IsIdenticalTo(before)" is added at the
             *       end of ConstantFolding_MulGrouping().
             *
             *       This algorithm could make it into (x*z*x)^0.5 * x^16,
             *       but this is wrong, for it falsely includes x^evenint.. twice.
             */
            for(size_t a=0; a<data.size(); ++a)
            {
                Value_t exp_a = data[a].first;
                if(fp_equal(exp_a, Value_t(1))) continue;
                for(size_t b=a+1; b<data.size(); ++b)
                {
                    Value_t exp_b = data[b].first;
                    Value_t exp_diff = exp_b - exp_a;
                    if(exp_diff >= fp_abs(exp_a)) break;
                    Value_t exp_diff_still_probable_integer = exp_diff * Value_t(16);
                    if(isInteger(exp_diff_still_probable_integer)
                    && !(isInteger(exp_b) && !isInteger(exp_diff))
                      )
                    {
                        /* When input is x^3 * z^2,
                         * exp_a = 2
                         * a_set = z
                         * exp_b = 3
                         * b_set = x
                         * exp_diff = 3-2 = 1
                         */
                        TreeSet& a_set = data[a].second;
                        TreeSet& b_set = data[b].second;
          #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "Before ConstantExponentCollection iteration:\n";
                        Dump(std::cout);
          #endif
                        if(isEvenInteger(exp_b)
                        //&& !isEvenInteger(exp_diff)
                        && !isEvenInteger(exp_diff+exp_a))
                        {
                            CodeTree<Value_t> tmp2;
                            tmp2.SetOpcode( cMul );
                            tmp2.SetParamsMove(b_set);
                            tmp2.Rehash();
                            CodeTree<Value_t> tmp;
                            tmp.SetOpcode( cAbs );
                            tmp.AddParamMove(tmp2);
                            tmp.Rehash();
                            b_set.resize(1);
                            b_set[0].swap(tmp);
                        }

                        a_set.insert(a_set.end(), b_set.begin(), b_set.end());

                        TreeSet b_copy = b_set;
                        data.erase(data.begin() + b);
                        MoveToSet_NonUnique(exp_diff, b_copy);
                        changed = true;

          #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "After ConstantExponentCollection iteration:\n";
                        Dump(std::cout);
          #endif
                        goto redo;
                    }
                }
            }
            return changed;
        }

    #ifdef DEBUG_SUBSTITUTIONS
        void Dump(std::ostream& out)
        {
            for(size_t a=0; a<data.size(); ++a)
            {
                out.precision(12);
                out << data[a].first << ": ";
                for(size_t b=0; b<data[a].second.size(); ++b)
                {
                    if(b > 0) out << '*';
                    DumpTree(data[a].second[b], out);
                }
                out << std::endl;
            }
        }
    #endif
    };

    template<typename Value_t>
    static CodeTree<Value_t> CollectMulGroup_Item(
        CodeTree<Value_t>& value,
        bool& has_highlevel_opcodes)
    {
        switch(value.GetOpcode())
        {
            case cPow:
            {
                CodeTree<Value_t> exponent = value.GetParam(1);
                value.Become( value.GetParam(0) );
                return exponent;
            }
            /* - disabled to avoid clashes with powi
            case cCbrt:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTreeImmed( Value_t(1) / Value_t(3) );
            case cSqrt:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTreeImmed( Value_t(0.5) );
            */
            case cRSqrt:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTreeImmed( Value_t(-0.5) );
            case cInv:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTreeImmed( Value_t(-1) );
            default: break;
        }
        return CodeTreeImmed( Value_t(1) );
    }

    template<typename Value_t>
    static void CollectMulGroup(
        CollectionSet<Value_t>& mul,
        const CodeTree<Value_t>& tree,
        const CodeTree<Value_t>& factor,
        bool& should_regenerate,
        bool& has_highlevel_opcodes
    )
    {
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            CodeTree<Value_t> value(tree.GetParam(a));

            CodeTree<Value_t> exponent ( CollectMulGroup_Item(value, has_highlevel_opcodes) );

            if(!factor.IsImmed() || factor.GetImmed() != Value_t(1.0))
            {
                CodeTree<Value_t> new_exp;
                new_exp.SetOpcode(cMul);
                new_exp.AddParam( exponent );
                new_exp.AddParam( factor );
                new_exp.Rehash();
                exponent.swap( new_exp );
            }
        #if 0 /* FIXME: This does not work */
            if(value.GetOpcode() == cMul)
            {
                if(1)
                {
                    // Avoid erroneously converting
                    //          (x*z)^0.5 * z^2
                    // into     x^0.5 * z^2.5
                    // It should be x^0.5 * abs(z)^2.5, but this is not a good conversion.
                    bool exponent_is_even = exponent.IsImmed() && isEvenInteger(exponent.GetImmed());

                    for(size_t b=0; b<value.GetParamCount(); ++b)
                    {
                        bool tmp=false;
                        CodeTree<Value_t> val(value.GetParam(b));
                        CodeTree<Value_t> exp(CollectMulGroup_Item(val, tmp));
                        if(exponent_is_even
                        || (exp.IsImmed() && isEvenInteger(exp.GetImmed())))
                        {
                            CodeTree<Value_t> new_exp;
                            new_exp.SetOpcode(cMul);
                            new_exp.AddParam(exponent);
                            new_exp.AddParamMove(exp);
                            new_exp.ConstantFolding();
                            if(!new_exp.IsImmed() || !isEvenInteger(new_exp.GetImmed()))
                            {
                                goto cannot_adopt_mul;
                            }
                        }
                    }
                }
                CollectMulGroup(mul, value, exponent,
                                should_regenerate,
                                has_highlevel_opcodes);
            }
            else cannot_adopt_mul:
        #endif
            {
                if(mul.AddCollection(value, exponent) == CollectionSetBase::Suboptimal)
                    should_regenerate = true;
            }
        }
    }

    template<typename Value_t>
    bool ConstantFolding_MulGrouping(CodeTree<Value_t>& tree)
    {
        bool has_highlevel_opcodes = false;
        bool should_regenerate = false;
        CollectionSet<Value_t> mul;

        CollectMulGroup(mul, tree, CodeTreeImmed(Value_t(1)),
                        should_regenerate,
                        has_highlevel_opcodes);

        typedef std::pair<CodeTree<Value_t>/*exponent*/,
                          std::vector<CodeTree<Value_t> >/*base value (mul group)*/
                         > exponent_list;
        typedef std::multimap<fphash_t,/*exponent hash*/
                              exponent_list> exponent_map;
        exponent_map by_exponent;

        for(typename CollectionSet<Value_t>::PositionType
            j = mul.collections.begin();
            j != mul.collections.end();
            ++j)
        {
            CodeTree<Value_t>& value = j->second.value;
            CodeTree<Value_t>& exponent = j->second.factor;
            if(j->second.factor_needs_rehashing) exponent.Rehash();
            const fphash_t exponent_hash = exponent.GetHash();

            typename exponent_map::iterator i = by_exponent.lower_bound(exponent_hash);
            for(; i != by_exponent.end() && i->first == exponent_hash; ++i)
                if(i->second.first.IsIdenticalTo(exponent))
                {
                    if(!exponent.IsImmed() || !fp_equal(exponent.GetImmed(), Value_t(1)))
                        should_regenerate = true;
                    i->second.second.push_back(value);
                    goto skip_b;
                }
            by_exponent.insert(i, std::make_pair(exponent_hash,
                std::make_pair(exponent,
                               std::vector<CodeTree<Value_t> > (size_t(1), value)
                              )));
        skip_b:;
        }

    #ifdef FP_MUL_COMBINE_EXPONENTS
        ConstantExponentCollection<Value_t> by_float_exponent;
        for(typename exponent_map::iterator
            j,i = by_exponent.begin();
            i != by_exponent.end();
            i=j)
        {
            j=i; ++j;
            exponent_list& list = i->second;
            if(list.first.IsImmed())
            {
                Value_t exponent = list.first.GetImmed();
                if(!(exponent == Value_t(0)))
                    by_float_exponent.MoveToSet_Unique(exponent, list.second);
                by_exponent.erase(i);
            }
        }
        if(by_float_exponent.Optimize())
            should_regenerate = true;
    #endif

        if(should_regenerate)
        {
            CodeTree<Value_t> before = tree;
            before.CopyOnWrite();

          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Before ConstantFolding_MulGrouping: "; DumpTree(before);
            std::cout << "\n";
          #endif
            tree.DelParams();

            /* Group by exponents */
            /* First handle non-constant exponents */
            for(typename exponent_map::iterator
                i = by_exponent.begin();
                i != by_exponent.end();
                ++i)
            {
                exponent_list& list = i->second;
        #ifndef FP_MUL_COMBINE_EXPONENTS
                if(list.first.IsImmed())
                {
                    Value_t exponent = list.first.GetImmed();
                    if(exponent == Value_t(0)) continue;
                    if(fp_equal(exponent, Value_t(1) ))
                    {
                        tree.AddParamsMove(list.second);
                        continue;
                    }
                }
        #endif
                CodeTree<Value_t> mul;
                mul.SetOpcode(cMul);
                mul.SetParamsMove( list.second);
                mul.Rehash();

                if(has_highlevel_opcodes && list.first.IsImmed())
                {
                    if(list.first.GetImmed() == Value_t(1) / Value_t(3))
                    {
                        CodeTree<Value_t> cbrt;
                        cbrt.SetOpcode(cCbrt);
                        cbrt.AddParamMove(mul);
                        cbrt.Rehash();
                        tree.AddParamMove(cbrt);
                        continue;
                    }
                    if(list.first.GetImmed() == Value_t(0.5) )
                    {
                        CodeTree<Value_t> sqrt;
                        sqrt.SetOpcode(cSqrt);
                        sqrt.AddParamMove(mul);
                        sqrt.Rehash();
                        tree.AddParamMove(sqrt);
                        continue;
                    }
                    if(list.first.GetImmed() == Value_t(-0.5) )
                    {
                        CodeTree<Value_t> rsqrt;
                        rsqrt.SetOpcode(cRSqrt);
                        rsqrt.AddParamMove(mul);
                        rsqrt.Rehash();
                        tree.AddParamMove(rsqrt);
                        continue;
                    }
                    if(list.first.GetImmed() == Value_t(-1))
                    {
                        CodeTree<Value_t> inv;
                        inv.SetOpcode(cInv);
                        inv.AddParamMove(mul);
                        inv.Rehash();
                        tree.AddParamMove(inv);
                        continue;
                    }
                }
                CodeTree<Value_t> pow;
                pow.SetOpcode(cPow);
                pow.AddParamMove(mul);
                pow.AddParamMove( list.first );
                pow.Rehash();
                tree.AddParamMove(pow);
            }
        #ifdef FP_MUL_COMBINE_EXPONENTS
            by_exponent.clear();
            /* Then handle constant exponents */
            for(size_t a=0; a<by_float_exponent.data.size(); ++a)
            {
                Value_t exponent = by_float_exponent.data[a].first;
                if(fp_equal(exponent, Value_t(1)))
                {
                    tree.AddParamsMove(by_float_exponent.data[a].second);
                    continue;
                }
                CodeTree<Value_t> mul;
                mul.SetOpcode(cMul);
                mul.SetParamsMove( by_float_exponent.data[a].second );
                mul.Rehash();
                CodeTree<Value_t> pow;
                pow.SetOpcode(cPow);
                pow.AddParamMove(mul);
                pow.AddParam( CodeTreeImmed( exponent ) );
                pow.Rehash();
                tree.AddParamMove(pow);
            }
        #endif
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "After ConstantFolding_MulGrouping: "; DumpTree(tree);
            std::cout << "\n";
          #endif
            // return true;
            return !tree.IsIdenticalTo(before); // avoids infinite looping
        }
        return false;
    }

    template<typename Value_t>
    bool ConstantFolding_AddGrouping(CodeTree<Value_t>& tree)
    {
        bool should_regenerate = false;
        CollectionSet<Value_t> add;
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            if(tree.GetParam(a).GetOpcode() == cMul) continue;
            if(add.AddCollection(tree.GetParam(a)) == CollectionSetBase::Suboptimal)
                should_regenerate = true;
            // This catches x + x and x - x
        }
        std::vector<bool> remaining ( tree.GetParamCount() );
        size_t has_mulgroups_remaining = 0;
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            const CodeTree<Value_t>& mulgroup = tree.GetParam(a);
            if(mulgroup.GetOpcode() == cMul)
            {
                // This catches x + y*x*z, producing x*(1 + y*z)
                //
                // However we avoid changing 7 + 7*x into 7*(x+1),
                // because it may lead us into producing code such
                // as 20*x + 50*(x+1) + 10, which would be much
                // better expressed as 70*x + 60, and converting
                // back to that format would be needlessly hairy.
                for(size_t b=0; b<mulgroup.GetParamCount(); ++b)
                {
                    if(mulgroup.GetParam(b).IsImmed()) continue;
                    typename CollectionSet<Value_t>::PositionType c
                        = add.FindIdenticalValueTo(mulgroup.GetParam(b));
                    if(add.Found(c))
                    {
                        CodeTree<Value_t> tmp(mulgroup, typename CodeTree<Value_t>::CloneTag());
                        tmp.DelParam(b);
                        tmp.Rehash();
                        add.AddCollectionTo(tmp, c);
                        should_regenerate = true;
                        goto done_a;
                    }
                }
                remaining[a]  = true;
                has_mulgroups_remaining += 1;
            done_a:;
            }
        }

        if(has_mulgroups_remaining > 0)
        {
            if(has_mulgroups_remaining > 1) // is it possible to find a duplicate?
            {
                std::vector< std::pair<CodeTree<Value_t>, size_t> > occurance_counts;
                std::multimap<fphash_t, size_t> occurance_pos;
                bool found_dup = false;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                    if(remaining[a])
                    {
                        // This catches x*a + x*b, producing x*(a+b)
                        for(size_t b=0; b<tree.GetParam(a).GetParamCount(); ++b)
                        {
                            const CodeTree<Value_t>& p = tree.GetParam(a).GetParam(b);
                            const fphash_t   p_hash = p.GetHash();
                            for(std::multimap<fphash_t, size_t>::const_iterator
                                i = occurance_pos.lower_bound(p_hash);
                                i != occurance_pos.end() && i->first == p_hash;
                                ++i)
                            {
                                if(occurance_counts[i->second].first.IsIdenticalTo(p))
                                {
                                    occurance_counts[i->second].second += 1;
                                    found_dup = true;
                                    goto found_mulgroup_item_dup;
                                }
                            }
                            occurance_counts.push_back(std::make_pair(p, size_t(1)));
                            occurance_pos.insert(std::make_pair(p_hash, occurance_counts.size()-1));
                        found_mulgroup_item_dup:;
                        }
                    }
                if(found_dup)
                {
                    // Find the "x" to group by
                    CodeTree<Value_t> group_by; { size_t max = 0;
                    for(size_t p=0; p<occurance_counts.size(); ++p)
                        if(occurance_counts[p].second <= 1)
                            occurance_counts[p].second = 0;
                        else
                        {
                            occurance_counts[p].second *= occurance_counts[p].first.GetDepth();
                            if(occurance_counts[p].second > max)
                                { group_by = occurance_counts[p].first; max = occurance_counts[p].second; }
                        } }
                    // Collect the items for adding in the group (a+b)
                    CodeTree<Value_t> group_add;
                    group_add.SetOpcode(cAdd);

        #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "Duplicate across some trees: ";
                    DumpTree(group_by);
                    std::cout << " in ";
                    DumpTree(tree);
                    std::cout << "\n";
        #endif
                    for(size_t a=0; a<tree.GetParamCount(); ++a)
                        if(remaining[a])
                            for(size_t b=0; b<tree.GetParam(a).GetParamCount(); ++b)
                                if(group_by.IsIdenticalTo(tree.GetParam(a).GetParam(b)))
                                {
                                    CodeTree<Value_t> tmp(tree.GetParam(a), typename CodeTree<Value_t>::CloneTag());
                                    tmp.DelParam(b);
                                    tmp.Rehash();
                                    group_add.AddParamMove(tmp);
                                    remaining[a] = false;
                                    break;
                                }
                    group_add.Rehash();
                    CodeTree<Value_t> group;
                    group.SetOpcode(cMul);
                    group.AddParamMove(group_by);
                    group.AddParamMove(group_add);
                    group.Rehash();
                    add.AddCollection(group);
                    should_regenerate = true;
                }
            }

            // all remaining mul-groups.
            for(size_t a=0; a<tree.GetParamCount(); ++a)
                if(remaining[a])
                {
                    if(add.AddCollection(tree.GetParam(a)) == CollectionSetBase::Suboptimal)
                        should_regenerate = true;
                }
        }

        if(should_regenerate)
        {
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Before ConstantFolding_AddGrouping: "; DumpTree(tree);
            std::cout << "\n";
          #endif
            tree.DelParams();

            for(typename CollectionSet<Value_t>::PositionType
                j = add.collections.begin();
                j != add.collections.end();
                ++j)
            {
                CodeTree<Value_t>& value = j->second.value;
                CodeTree<Value_t>& coeff = j->second.factor;
                if(j->second.factor_needs_rehashing) coeff.Rehash();

                if(coeff.IsImmed())
                {
                    if(fp_equal(coeff.GetImmed(), Value_t(0)))
                        continue;
                    if(fp_equal(coeff.GetImmed(), Value_t(1)))
                    {
                        tree.AddParamMove(value);
                        continue;
                    }
                }
                CodeTree<Value_t> mul;
                mul.SetOpcode(cMul);
                mul.AddParamMove(value);
                mul.AddParamMove(coeff);
                mul.Rehash();
                tree.AddParamMove(mul);
            }
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "After ConstantFolding_AddGrouping: "; DumpTree(tree);
            std::cout << "\n";
          #endif
            return true;
        }
        return false;
    }
}
