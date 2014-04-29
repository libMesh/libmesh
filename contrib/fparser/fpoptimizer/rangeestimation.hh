#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH

#include "codetree.hh"
#include "valuerange.hh"

namespace FPoptimizer_CodeTree
{
    enum TriTruthValue { IsAlways, IsNever, Unknown };

    /* This function calculates the minimum and maximum values
     * of the tree's result. If an estimate cannot be made,
     * -inf..+inf is assumed (min.known=max.known=false).
     */
    template<typename Value_t>
    range<Value_t> CalculateResultBoundaries(const CodeTree<Value_t>& tree);

    template<typename Value_t>
    bool IsLogicalValue(const CodeTree<Value_t>& tree);

    template<typename Value_t>
    TriTruthValue GetIntegerInfo(const CodeTree<Value_t>& tree);

    template<typename Value_t>
    inline TriTruthValue GetEvennessInfo(const CodeTree<Value_t>& tree)
    {
        if(!tree.IsImmed()) return Unknown;
        const Value_t& value = tree.GetImmed();
        if(FUNCTIONPARSERTYPES::isEvenInteger(value)) return IsAlways;
        if(FUNCTIONPARSERTYPES::isOddInteger(value)) return IsNever;
        return Unknown;
    }

    template<typename Value_t>
    inline TriTruthValue GetPositivityInfo(const CodeTree<Value_t>& tree)
    {
        range<Value_t> p = CalculateResultBoundaries(tree);
        if(p.min.known && p.min.val >= Value_t()) return IsAlways;
        if(p.max.known && p.max.val <  Value_t()) return IsNever;
        return Unknown;
    }

    template<typename Value_t>
    inline TriTruthValue GetLogicalValue(const CodeTree<Value_t>& tree, bool abs)
    {
        range<Value_t> p = CalculateResultBoundaries(tree);
        if(IsLogicalTrueValue(p, abs)) return IsAlways;
        if(IsLogicalFalseValue(p, abs)) return IsNever;
        return Unknown;
    }
}

#endif
