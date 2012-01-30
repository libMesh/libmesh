#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include "grammar.hh"
#include "consts.hh"
#include "opcodename.hh"
#include "optimize.hh"

#include <stdio.h>

#include <algorithm>
#include <map>
#include <sstream>

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_Grammar;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace
{
    /* I have heard that std::equal_range() is practically worthless
     * due to the insane limitation that the two parameters for Comp() must
     * be of the same type. Hence we must reinvent the wheel and implement
     * our own here. This is practically identical to the one from
     * GNU libstdc++, except rewritten. -Bisqwit
     */
    template<typename It, typename T, typename Comp>
    std::pair<It, It>
    MyEqualRange(It first, It last, const T& val, Comp comp)
    {
      /*while(first != last && comp(*first, val)) ++first;
        while(first != last)
        {
            It temp = last; --temp;
            if(!comp(val, *temp)) break;
            last = temp;
        }
        return std::pair<It,It>(first,last);*/

        size_t len = last-first;
        while(len > 0)
        {
            size_t half = len/2;
            It middle(first); middle += half;
            if(comp(*middle, val))
            {
                first = middle;
                ++first;
                len = len - half - 1;
            }
            else if(comp(val, *middle))
            {
                len = half;
            }
            else
            {
                // The following implements this:
                // // left = lower_bound(first, middle, val, comp);
                It left(first);
              {///
                It& first2 = left;
                It last2(middle);
                size_t len2 = last2-first2;
                while(len2 > 0)
                {
                    size_t half2 = len2 / 2;
                    It middle2(first2); middle2 += half2;
                    if(comp(*middle2, val))
                    {
                        first2 = middle2;
                        ++first2;
                        len2 = len2 - half2 - 1;
                    }
                    else
                        len2 = half2;
                }
                // left = first2;  - not needed, already happens due to reference
              }///
                first += len;
                // The following implements this:
                // // right = upper_bound(++middle, first, val, comp);
                It right(++middle);
              {///
                It& first2 = right;
                It& last2 = first;
                size_t len2 = last2-first2;
                while(len2 > 0)
                {
                    size_t half2 = len2 / 2;
                    It middle2(first2); middle2 += half2;
                    if(comp(val, *middle2))
                        len2 = half2;
                    else
                    {
                        first2 = middle2;
                        ++first2;
                        len2 = len2 - half2 - 1;
                    }
                }
                // right = first2;  - not needed, already happens due to reference
              }///
                return std::pair<It,It> (left,right);
            }
        }
        return std::pair<It,It> (first,first);
    }

    /* A helper for std::equal_range */
    template<typename Value_t>
    struct OpcodeRuleCompare
    {
        bool operator() (const CodeTree<Value_t>& tree,
                         unsigned rulenumber) const
        {
            /* If this function returns true, len=half.
             */
            const Rule& rule = grammar_rules[rulenumber];
            return tree.GetOpcode() < rule.match_tree.subfunc_opcode;
        }
        bool operator() (unsigned rulenumber,
                         const CodeTree<Value_t>& tree) const
        {
            /* If this function returns true, rule will be excluded from the equal_range
             */
            const Rule& rule = grammar_rules[rulenumber];
            return rule.match_tree.subfunc_opcode < tree.GetOpcode();
        }
    };

    /* Test and apply a rule to a given CodeTree */
    template<typename Value_t>
    bool TestRuleAndApplyIfMatch(
        const Rule& rule,
        CodeTree<Value_t>& tree,
        bool from_logical_context)
    {
        MatchInfo<Value_t> info;

        MatchResultType found(false, MatchPositionSpecBaseP());

        if((rule.situation_flags & LogicalContextOnly)
        && !from_logical_context)
        {
            /* If the rule only applies in logical contexts,
             * but we do not have a logical context, fail the rule
             */
            goto fail;
        }
        if(FUNCTIONPARSERTYPES::IsIntType<Value_t>::result)
        {
            if(rule.situation_flags & NotForIntegers)
                goto fail;
        }
        else
        {
            if(rule.situation_flags & OnlyForIntegers)
                goto fail;
        }
        if(FUNCTIONPARSERTYPES::IsComplexType<Value_t>::result)
        {
            if(rule.situation_flags & NotForComplex)
                goto fail;
        }
        else
        {
            if(rule.situation_flags & OnlyForComplex)
                goto fail;
        }

        /*std::cout << "TESTING: ";
        DumpMatch(rule, *tree, info, false);*/

        for(;;)
        {
        #ifdef DEBUG_SUBSTITUTIONS
            //DumpMatch(rule, tree, info, "Testing");
        #endif
            found = TestParams(rule.match_tree, tree, found.specs, info, true);
            if(found.found) break;
            if(!&*found.specs)
            {
            fail:;
                // Did not match
        #ifdef DEBUG_SUBSTITUTIONS
                DumpMatch(rule, tree, info, false);
        #endif
                return false;
            }
        }
        // Matched
    #ifdef DEBUG_SUBSTITUTIONS
        DumpMatch(rule, tree, info, true);
    #endif
        SynthesizeRule(rule, tree, info);
        return true;
    }
}

namespace FPoptimizer_Optimize
{
    /* Apply the grammar to a given CodeTree */
    template<typename Value_t>
    bool ApplyGrammar(
        const Grammar& grammar,
        CodeTree<Value_t>& tree,
        bool from_logical_context)
    {
        if(tree.GetOptimizedUsing() == &grammar)
        {
#ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Already optimized:  ";
            DumpTree(tree);
            std::cout << "\n" << std::flush;
#endif
            return false;
        }

        /* First optimize all children */
        if(true)
        {
            bool changed = false;

            switch(tree.GetOpcode())
            {
                case cNot:
                case cNotNot:
                case cAnd:
                case cOr:
                    for(size_t a=0; a<tree.GetParamCount(); ++a)
                        if(ApplyGrammar( grammar, tree.GetParam(a), true))
                            changed = true;
                    break;
                case cIf:
                case cAbsIf:
                    if(ApplyGrammar( grammar, tree.GetParam(0), tree.GetOpcode() == cIf))
                        changed = true;
                    for(size_t a=1; a<tree.GetParamCount(); ++a)
                        if(ApplyGrammar( grammar, tree.GetParam(a), from_logical_context))
                            changed = true;
                    break;
                default:
                    for(size_t a=0; a<tree.GetParamCount(); ++a)
                        if(ApplyGrammar( grammar, tree.GetParam(a), false))
                            changed = true;
            }

            if(changed)
            {
                // Give the parent node a rerun at optimization
                tree.Mark_Incompletely_Hashed();
                return true;
            }
        }

        /* Figure out which rules _may_ match this tree */
        typedef const unsigned short* rulenumit;

        /*std::cout << "---TEST:";
        for(unsigned n=0; n<grammar.rule_count; ++n)
            std::cout << ' ' << (unsigned)grammar.rule_list[n];
        std::cout << "\n";*/

        std::pair<rulenumit, rulenumit> range =
            MyEqualRange(grammar.rule_list,
                         grammar.rule_list + grammar.rule_count,
                         tree,
                         OpcodeRuleCompare<Value_t> ());

        std::vector<unsigned short> rules;
        rules.reserve(range.second - range.first);
        for(rulenumit r = range.first; r != range.second; ++r)
        {
            //if(grammar_rules[*r].match_tree.subfunc_opcode != tree.GetOpcode()) continue;
            if(IsLogisticallyPlausibleParamsMatch(grammar_rules[*r].match_tree, tree))
                rules.push_back(*r);
        }
        range.first  = !rules.empty() ? &rules[0]                : 0;
        range.second = !rules.empty() ? &rules[rules.size()-1]+1 : 0;

        if(range.first != range.second)
        {
#ifdef DEBUG_SUBSTITUTIONS
            if(range.first != range.second)
            {
                std::cout << "Input (" << FP_GetOpcodeName(tree.GetOpcode())
                          << ")[" << tree.GetParamCount()
                          << "]";
                if(from_logical_context)
                    std::cout << "(Logical)";

                unsigned first=~unsigned(0), prev=~unsigned(0);
                const char* sep = ", rules ";
                for(rulenumit r = range.first; r != range.second; ++r)
                {
                    if(first==~unsigned(0)) first=prev=*r;
                    else if(*r == prev+1) prev=*r;
                    else
                    {
                        std::cout << sep << first; sep=",";
                        if(prev != first) std::cout << '-' << prev;
                        first = prev = *r;
                    }
                }
                if(first != ~unsigned(0))
                {
                    std::cout << sep << first;
                    if(prev != first) std::cout << '-' << prev;
                }
                std::cout << ": ";
                DumpTree(tree);
                std::cout << "\n" << std::flush;
            }
#endif

            bool changed = false;

            for(rulenumit r = range.first; r != range.second; ++r)
            {
            #ifndef DEBUG_SUBSTITUTIONS
                if(!IsLogisticallyPlausibleParamsMatch(grammar_rules[*r].match_tree, tree))
                    continue;
            #endif
                if(TestRuleAndApplyIfMatch(grammar_rules[*r], tree, from_logical_context))
                {
                    changed = true;
                    break;
                }
            }

            if(changed)
            {
    #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Changed." << std::endl;
                std::cout << "Output: ";
                DumpTree(tree);
                std::cout << "\n" << std::flush;
    #endif
                // Give the parent node a rerun at optimization
                tree.Mark_Incompletely_Hashed();
                return true;
            }
        }

        // No changes, consider the tree properly optimized.
        tree.SetOptimizedUsing(&grammar);
        return false;
    }

    // This function (void cast) helps avoid a type punning warning from GCC.
    template<typename Value_t>
    bool ApplyGrammar(const void* p, FPoptimizer_CodeTree::CodeTree<Value_t>& tree)
    {
        return ApplyGrammar( *(const Grammar*) p, tree);
    }

    template<typename Value_t>
    void ApplyGrammars(FPoptimizer_CodeTree::CodeTree<Value_t>& tree)
    {
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round1\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_round1, tree))
            { //std::cout << "Rerunning 1\n";
                tree.FixIncompleteHashes();
            }

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round2\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_round2, tree))
            { //std::cout << "Rerunning 2\n";
                tree.FixIncompleteHashes();
            }

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round3\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_round3, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }

        #ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_nonshortcut_logical_evaluation\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_nonshortcut_logical_evaluation, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }
        #endif

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round4\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_round4, tree))
            { //std::cout << "Rerunning 4\n";
                tree.FixIncompleteHashes();
            }

        #ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_shortcut_logical_evaluation\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_shortcut_logical_evaluation, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }
        #endif

        #ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_ignore_if_sideeffects\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_ignore_if_sideeffects, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }
        #endif

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_abslogical\n";
        #endif
        while(ApplyGrammar((const void*)&grammar_optimize_abslogical, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }

        #undef C
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_Optimize
{
#define FP_INSTANTIATE(type) \
    template void ApplyGrammars(FPoptimizer_CodeTree::CodeTree<type>& tree);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
