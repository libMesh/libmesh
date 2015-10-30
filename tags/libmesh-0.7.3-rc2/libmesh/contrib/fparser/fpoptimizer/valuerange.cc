#include "valuerange.hh"

#ifdef FP_SUPPORT_OPTIMIZER

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void range<Value_t>::set_abs()
    {
        using namespace FUNCTIONPARSERTYPES;
        bool has_negative = !min.known || min.val < Value_t();
        bool has_positive = !max.known || max.val > Value_t();
        bool crosses_axis = has_negative && has_positive;

        rangehalf<Value_t> newmax;              //  ..+inf
        if(min.known && max.known)              //  ..N
            newmax.set( fp_max(fp_abs(min.val), fp_abs(max.val)) );

        if(crosses_axis)
            min.set( Value_t() );               // 0..
        else
        {
            // Does not cross axis, so choose the smallest of known values
            // (Either value is known; otherwise it would cross axis)
            if(min.known && max.known)          // N..
                min.set( fp_min(fp_abs(min.val), fp_abs(max.val)) );
            else if(min.known)
                min.set( fp_abs(min.val) );
            else //if(max.known)
                min.set( fp_abs(max.val) );
        }
        max = newmax;
    }

    template<typename Value_t>
    void range<Value_t>::set_neg()
    {
        std::swap(min, max);
        min.val = -min.val;
        max.val = -max.val;
    }

    template<typename Value_t>
    bool IsLogicalTrueValue(const range<Value_t>& p, bool abs)
    {
        if(FUNCTIONPARSERTYPES::IsIntType<Value_t>::result)
        {
            if(p.min.known && p.min.val >= Value_t(1)) return true;
            if(!abs && p.max.known && p.max.val <= Value_t(-1)) return true;
        }
        else
        {
            if(p.min.known && p.min.val >= Value_t(0.5)) return true;
            if(!abs && p.max.known && p.max.val <= Value_t(-0.5)) return true;
        }
        return false;
    }

    template<typename Value_t>
    bool IsLogicalFalseValue(const range<Value_t>& p, bool abs)
    {
        if(FUNCTIONPARSERTYPES::IsIntType<Value_t>::result)
        {
            if(abs)
                return p.max.known && p.max.val < Value_t(1);
            else
                return p.min.known && p.max.known
                  && p.min.val > Value_t(-1) && p.max.val < Value_t(1);
        }
        else
        {
            if(abs)
                return p.max.known && p.max.val < Value_t(0.5);
            else
                return p.min.known && p.max.known
                   && p.min.val > Value_t(-0.5) && p.max.val < Value_t(0.5);
        }
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template struct range<type>; \
    template bool IsLogicalTrueValue(const range<type> &, bool); \
    template bool IsLogicalFalseValue(const range<type> &, bool);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
