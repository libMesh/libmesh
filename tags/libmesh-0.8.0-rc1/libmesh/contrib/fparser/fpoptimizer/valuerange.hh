#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH

#include "fparser.hh"
#include "extrasrc/fpaux.hh"

namespace FPoptimizer_CodeTree
{
    namespace rangeutil
    {
        template<unsigned Compare> struct Comp { };
        template<>struct Comp<FUNCTIONPARSERTYPES::cLess> {
            template<typename Value_t>
            inline bool operator() (const Value_t& a, const Value_t& b) { return a<b; }
        };
        template<>struct Comp<FUNCTIONPARSERTYPES::cLessOrEq> {
            template<typename Value_t>
            inline bool operator() (const Value_t& a, const Value_t& b) { return a<=b; }
        };
        template<>struct Comp<FUNCTIONPARSERTYPES::cGreater> {
            template<typename Value_t>
            inline bool operator() (const Value_t& a, const Value_t& b) { return a>b; }
        };
        template<>struct Comp<FUNCTIONPARSERTYPES::cGreaterOrEq> {
            template<typename Value_t>
            inline bool operator() (const Value_t& a, const Value_t& b) { return a>=b; }
        };
        template<>struct Comp<FUNCTIONPARSERTYPES::cEqual> {
            template<typename Value_t>
            inline bool operator() (const Value_t& a, const Value_t& b) { return a==b; }
        };
        template<>struct Comp<FUNCTIONPARSERTYPES::cNEqual> {
            template<typename Value_t>
            inline bool operator() (const Value_t& a, const Value_t& b) { return a!=b; }
        };
    }

    template<typename Value_t>
    struct rangehalf
    {
        Value_t val;
        bool    known;

        rangehalf(): val(), known(false) { }
        rangehalf(const Value_t& v) : val(v), known(true) { }

        inline void set(const Value_t& v) { known=true; val=v; }

        /////////

        /* If value is known, refine it using func.
         * Otherwise, use the model.
         *
         * Call like this:
         *      range := param
         *      range.min.set(fp_floor)
         *      If param is known, sets minimum to floor(param.min)
         *      Otherwise, sets min to unknown
         * Or:
         *      range := param
         *      range.min.set(fp_atan, -pihalf)
         *      If param is known, sets minimum to floor(param.min)
         *      Otherwise, sets min to -pihalf
         */
        void set
            (Value_t (*const func)(Value_t),
             rangehalf<Value_t> model = rangehalf<Value_t>())
        {
            if(known) val = func(val); else *this = model;
        }

        void set
            (Value_t (*const func)(const Value_t&),
             rangehalf<Value_t> model = rangehalf<Value_t>())
        {
            if(known) val = func(val); else *this = model;
        }

        /* Call like this:
         *      range := param
         *      range.min.set_if<cGreater>(-1, fp_asin, -pihalf)
         *      If param is known AND param.min > -1, sets minimum to asin(param.min)
         *      Otherwise, sets min to -pihalf
         *      The purpose of the condition is to ensure that the function
         *      is not being called with illegal values.
         */
        template<unsigned Compare>
        void set_if
            (Value_t v,
             Value_t (*const func)(Value_t),
             rangehalf<Value_t> model = rangehalf<Value_t>())
        {
            if(known && rangeutil::Comp<Compare>() (val,v))
                val = func(val);
            else
                *this = model;
        }
        template<unsigned Compare>
        void set_if
            (const Value_t& v,
             Value_t (*const func)(const Value_t&),
             rangehalf<Value_t> model = rangehalf<Value_t>())
        {
            if(known && rangeutil::Comp<Compare>() (val,v))
                val = func(val);
            else
                *this = model;
        }
    };

    /* range expresses the range of values that an expression can take. */
    template<typename Value_t>
    struct range
    {
        rangehalf<Value_t> min, max;

        /* Initializations */
        range() : min(),max() { }
        range(Value_t mi,Value_t ma): min(mi),max(ma) { }
        range(bool,Value_t ma): min(),max(ma) { }
        range(Value_t mi,bool): min(mi),max() { }

        /* Apply the abs() function to the range,
         * i.e. +3..+5 becomes +3..+5;
         *      -3..+5 becomes  0..+5;
         *      -3..-1 becomes  0..+1
         */
        void set_abs();

        /* Negate the range, i.e. -3..+5 becomes -5..+3 */
        void set_neg();
    };

    /* Analysis functions for a range */
    template<typename Value_t>
    bool IsLogicalTrueValue(const range<Value_t>& p, bool abs);

    template<typename Value_t>
    bool IsLogicalFalseValue(const range<Value_t>& p, bool abs);
}

#endif
