#include "rangeestimation.hh"
#include "consts.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_CodeTree;

//#define DEBUG_SUBSTITUTIONS_extra_verbose

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    range<Value_t> CalculateResultBoundaries(const CodeTree<Value_t>& tree)
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
    {
        using namespace FUNCTIONPARSERTYPES;
        range<Value_t> tmp = CalculateResultBoundaries_do(tree);
        std::cout << "Estimated boundaries: ";
        if(tmp.min.known) std::cout << tmp.min.val; else std::cout << "-inf";
        std::cout << " .. ";
        if(tmp.max.known) std::cout << tmp.max.val; else std::cout << "+inf";
        std::cout << ": ";
        DumpTree(tree);
        std::cout << std::endl;
        return tmp;
    }
    template<typename Value_t>
    range<Value_t> CodeTree<Value_t>::CalculateResultBoundaries_do(const CodeTree<Value_t>& tree)
#endif
    {
        static const range<Value_t> pihalf_limits
            (-fp_const_pihalf<Value_t>(),
              fp_const_pihalf<Value_t>());

        static const range<Value_t> pi_limits
            (-fp_const_pi<Value_t>(),
              fp_const_pi<Value_t>());

        static const range<Value_t> abs_pi_limits
            ( Value_t(0),
              fp_const_pi<Value_t>());

        static const range<Value_t> plusminus1_limits
            ( Value_t(-1),
              Value_t(1) );

        using namespace std;
        switch( tree.GetOpcode() )
        {
            case cImmed:
                return range<Value_t>(tree.GetImmed(), tree.GetImmed()); // a definite value.
            case cAnd:
            case cAbsAnd:
            case cOr:
            case cAbsOr:
            case cNot:
            case cAbsNot:
            case cNotNot:
            case cAbsNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
            {
                /* These operations always produce truth values (0 or 1) */
                /* Narrowing them down is a matter of performing Constant optimization */
                return range<Value_t>( Value_t(0), Value_t(1) );
            }
            case cAbs:
            {
                /* cAbs always produces a positive value */
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.set_abs();
                return m;
            }

            case cLog: /* Defined for 0.0 < x <= inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.template set_if<cGreater>(Value_t(0), fp_log); // No boundaries
                m.max.template set_if<cGreater>(Value_t(0), fp_log); // No boundaries
                return m;
            }

            case cLog2: /* Defined for 0.0 < x <= inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.template set_if<cGreater>(Value_t(0), fp_log2); // No boundaries
                m.max.template set_if<cGreater>(Value_t(0), fp_log2); // No boundaries
                return m;
            }

            case cLog10: /* Defined for 0.0 < x <= inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.template set_if<cGreater>(Value_t(0), fp_log10); // No boundaries
                m.max.template set_if<cGreater>(Value_t(0), fp_log10); // No boundaries
                return m;
            }

            case cAcosh: /* defined for             1.0 <= x <= inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.template set_if<cGreaterOrEq>(Value_t(1), fp_acosh); // No boundaries
                m.max.template set_if<cGreaterOrEq>(Value_t(1), fp_acosh); // No boundaries
                return m;
            }
            case cAsinh: /* defined for all values -inf <= x <= inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_asinh); // No boundaries
                m.max.set(fp_asinh); // No boundaries
                return m;
            }
            case cAtanh: /* defined for -1.0 <= x < 1, results within -inf..+inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.template set_if<cGreater> (Value_t(-1), fp_atanh);
                m.max.template set_if<cLess>    (Value_t( 1), fp_atanh);
                return m;
            }
            case cAcos: /* defined for -1.0 <= x <= 1, results within CONSTANT_PI..0 */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                return range<Value_t>( // Note that the range is flipped!
                    (m.max.known && (m.max.val) < Value_t(1))
                        ? fp_acos(m.max.val) : Value_t(0),
                    (m.min.known && (m.min.val) >= Value_t(-1))
                        ? fp_acos(m.min.val) : fp_const_pi<Value_t>()
                                          );
            }
            case cAsin: /* defined for -1.0 <= x < 1, results within -CONSTANT_PIHALF..CONSTANT_PIHALF */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                /* Assuming that x is never outside valid limits */
                m.min.template set_if<cGreater>(Value_t(-1), fp_asin, pihalf_limits.min.val);
                m.max.template set_if<cLess   >(Value_t( 1), fp_asin, pihalf_limits.max.val);
                return m;
            }
            case cAtan: /* defined for all values -inf <= x <= inf, results within -CONSTANT_PIHALF..CONSTANT_PIHALF */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_atan, pihalf_limits.min.val);
                m.max.set(fp_atan, pihalf_limits.max.val);
                return m;
            }
            case cAtan2: /* too complicated to estimate */
            {
                //range<Value_t> p0 = CalculateResultBoundaries( tree.GetParam(0) );
                //range<Value_t> p1 = CalculateResultBoundaries( tree.GetParam(1) );
                if(tree.GetParam(0).IsImmed()
                && fp_equal(tree.GetParam(0).GetImmed(), Value_t(0)))   // y == 0
                {
                    // Either 0.0 or CONSTANT_PI
                    return abs_pi_limits;
                }
                if(tree.GetParam(1).IsImmed()
                && fp_equal(tree.GetParam(1).GetImmed(), Value_t(0)))   // x == 0
                {
                    // Either -CONSTANT_PIHALF or +CONSTANT_PIHALF
                    return pihalf_limits;
                }
                // Anything else
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return pi_limits;
            }

            case cSin:
            {
                /* Quite difficult to estimate due to the cyclic nature of the function. */
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                bool covers_full_cycle
                    = !m.min.known || !m.max.known
                    || (m.max.val - m.min.val) >= (fp_const_twopi<Value_t>());
                if(covers_full_cycle)
                    return range<Value_t>(Value_t(-1), Value_t(1));
                Value_t min = fp_mod(m.min.val, fp_const_twopi<Value_t>()); if(min<Value_t(0)) min+=fp_const_twopi<Value_t>();
                Value_t max = fp_mod(m.max.val, fp_const_twopi<Value_t>()); if(max<Value_t(0)) max+=fp_const_twopi<Value_t>();
                if(max < min) max += fp_const_twopi<Value_t>();
                bool covers_plus1  = (min <= fp_const_pihalf<Value_t>() && max >= fp_const_pihalf<Value_t>());
                bool covers_minus1 = (min <= Value_t(1.5)*fp_const_pi<Value_t>() && max >= Value_t(1.5)*fp_const_pi<Value_t>());
                if(covers_plus1 && covers_minus1)
                    return range<Value_t>(Value_t(-1), Value_t(1));
                if(covers_minus1)
                    return range<Value_t>(Value_t(-1), fp_max(fp_sin(min), fp_sin(max)));
                if(covers_plus1)
                    return range<Value_t>(fp_min(fp_sin(min), fp_sin(max)), Value_t(1));
                return range<Value_t>(fp_min(fp_sin(min), fp_sin(max)),
                                           fp_max(fp_sin(min), fp_sin(max)));
            }
            case cCos:
            {
                /* Quite difficult to estimate due to the cyclic nature of the function. */
                /* cos(x) = sin(pi/2 - x) = sin(x + pi/2) */
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                if(m.min.known) m.min.val += fp_const_pihalf<Value_t>();/*for cCos*/
                if(m.max.known) m.max.val += fp_const_pihalf<Value_t>();/*for cCos*/
                bool covers_full_cycle
                    = !m.min.known || !m.max.known
                    || (m.max.val - m.min.val) >= (fp_const_twopi<Value_t>());
                if(covers_full_cycle)
                    return range<Value_t>(Value_t(-1), Value_t(1));
                Value_t min = fp_mod(m.min.val, fp_const_twopi<Value_t>()); if(min<Value_t(0)) min+=fp_const_twopi<Value_t>();
                Value_t max = fp_mod(m.max.val, fp_const_twopi<Value_t>()); if(max<Value_t(0)) max+=fp_const_twopi<Value_t>();
                if(max < min) max += fp_const_twopi<Value_t>();
                bool covers_plus1  = (min <= fp_const_pihalf<Value_t>() && max >= fp_const_pihalf<Value_t>());
                bool covers_minus1 = (min <= Value_t(1.5)*fp_const_pi<Value_t>() && max >= Value_t(1.5)*fp_const_pi<Value_t>());
                if(covers_plus1 && covers_minus1)
                    return range<Value_t>(Value_t(-1), Value_t(1));
                if(covers_minus1)
                    return range<Value_t>(Value_t(-1), fp_max(fp_sin(min), fp_sin(max)));
                if(covers_plus1)
                    return range<Value_t>(fp_min(fp_sin(min), fp_sin(max)), Value_t(1));
                return range<Value_t>(fp_min(fp_sin(min), fp_sin(max)),
                                           fp_max(fp_sin(min), fp_sin(max)));
            }
            case cTan:
            {
                /* Could be narrowed down from here,
                 * but it's too complicated due to
                 * the cyclic nature of the function */
                /* TODO: A resourceful programmer may add it later. */
                return range<Value_t>(); // (CONSTANT_NEG_INF, CONSTANT_POS_INF);
            }

            case cCeil:
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.max.set(fp_ceil); // ceil() may increase the value, may not decrease
                return m;
            }
            case cFloor:
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_floor); // floor() may decrease the value, may not increase
                return m;
            }
            case cTrunc:
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_floor); // trunc() may either increase or decrease the value
                m.max.set(fp_ceil); // for safety, we assume both
                return m;
            }
            case cInt:
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_floor); // int() may either increase or decrease the value
                m.max.set(fp_ceil); // for safety, we assume both
                return m;
            }
            case cSinh: /* defined for all values -inf <= x <= inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_sinh); // No boundaries
                m.max.set(fp_sinh); // No boundaries
                return m;
            }
            case cTanh: /* defined for all values -inf <= x <= inf, results within -1..1 */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_tanh, plusminus1_limits.min);
                m.max.set(fp_tanh, plusminus1_limits.max);
                return m;
            }
            case cCosh: /* defined for all values -inf <= x <= inf, results within 1..inf */
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                if(m.min.known)
                {
                    if(m.max.known) // max, min
                    {
                        if(m.min.val >= Value_t(0) && m.max.val >= Value_t(0)) // +x .. +y
                            { m.min.val = fp_cosh(m.min.val); m.max.val = fp_cosh(m.max.val); }
                        else if((m.min.val) < Value_t(0) && m.max.val >= Value_t(0)) // -x .. +y
                            { Value_t tmp = fp_cosh(m.min.val); m.max.val = fp_cosh(m.max.val);
                              if(tmp > m.max.val) m.max.val = tmp;
                              m.min.val = Value_t(1); }
                        else // -x .. -y
                            { m.min.val = fp_cosh(m.min.val); m.max.val = fp_cosh(m.max.val);
                              std::swap(m.min.val, m.max.val); }
                    }
                    else // min, no max
                    {
                        if(m.min.val >= Value_t(0)) // 0..inf -> 1..inf
                            { m.max.known = false; m.min.val = fp_cosh(m.min.val); }
                        else
                            { m.max.known = false; m.min.val = Value_t(1); } // Anything between 1..inf
                    }
                }
                else // no min
                {
                    m.min.known = true; m.min.val = Value_t(1); // always a lower boundary
                    if(m.max.known) // max, no min
                    {
                        m.min.val = fp_cosh(m.max.val); // n..inf
                        m.max.known = false; // No upper boundary
                    }
                    else // no max, no min
                        m.max.known = false; // No upper boundary
                }
                return m;
            }

            case cIf:
            case cAbsIf:
            {
                // No guess which branch is chosen. Produce a spanning min & max.
                range<Value_t> res1 = CalculateResultBoundaries( tree.GetParam(1) );
                range<Value_t> res2 = CalculateResultBoundaries( tree.GetParam(2) );
                if(!res2.min.known) res1.min.known = false; else if(res1.min.known && (res2.min.val) < res1.min.val) res1.min.val = res2.min.val;
                if(!res2.max.known) res1.max.known = false; else if(res1.max.known && (res2.max.val) > res1.max.val) res1.max.val = res2.max.val;
                return res1;
            }

            case cMin:
            {
                bool has_unknown_min = false;
                bool has_unknown_max = false;

                range<Value_t> result;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    range<Value_t> m = CalculateResultBoundaries( tree.GetParam(a) );
                    if(!m.min.known)
                        has_unknown_min = true;
                    else if(!result.min.known || (m.min.val) < result.min.val)
                        result.min.val = m.min.val;

                    if(!m.max.known)
                        has_unknown_max = true;
                    else if(!result.max.known || (m.max.val) < result.max.val)
                        result.max.val = m.max.val;
                }
                if(has_unknown_min) result.min.known = false;
                if(has_unknown_max) result.max.known = false;
                return result;
            }
            case cMax:
            {
                bool has_unknown_min = false;
                bool has_unknown_max = false;

                range<Value_t> result;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    range<Value_t> m = CalculateResultBoundaries( tree.GetParam(a) );
                    if(!m.min.known)
                        has_unknown_min = true;
                    else if(!result.min.known || m.min.val > result.min.val)
                        result.min.val = m.min.val;

                    if(!m.max.known)
                        has_unknown_max = true;
                    else if(!result.max.known || m.max.val > result.max.val)
                        result.max.val = m.max.val;
                }
                if(has_unknown_min) result.min.known = false;
                if(has_unknown_max) result.max.known = false;
                return result;
            }
            case cAdd:
            {
                /* It's complicated. Follow the logic below. */
                /* Note: This also deals with the following opcodes:
                 *       cNeg, cSub, cRSub
                 */
                range<Value_t> result(Value_t(0), Value_t(0));
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    range<Value_t> item = CalculateResultBoundaries( tree.GetParam(a) );

                    if(item.min.known) result.min.val += item.min.val;
                    else             result.min.known = false;
                    if(item.max.known) result.max.val += item.max.val;
                    else             result.max.known = false;

                    if(!result.min.known && !result.max.known) break; // hopeless
                }
                if(result.min.known && result.max.known
                && result.min.val > result.max.val) std::swap(result.min.val, result.max.val);
                return result;
            }
            case cMul:
            {
                /* It's very complicated. Follow the logic below. */
                struct Value
                {
                    enum ValueType { Finite, MinusInf, PlusInf };
                    ValueType valueType;
                    Value_t value;

                    Value(ValueType t): valueType(t), value(0) {}
                    Value(Value_t v): valueType(Finite), value(v) {}

                    bool isNegative() const
                    {
                        return valueType == MinusInf ||
                            (valueType == Finite && value < Value_t(0));
                    }

                    void operator*=(const Value& rhs)
                    {
                        if(valueType == Finite && rhs.valueType == Finite)
                            value *= rhs.value;
                        else
                            valueType = (isNegative() != rhs.isNegative() ?
                                         MinusInf : PlusInf);
                    }

                    bool operator<(const Value& rhs) const
                    {
                        return
                            (valueType == MinusInf && rhs.valueType != MinusInf) ||
                            (valueType == Finite &&
                             (rhs.valueType == PlusInf ||
                              (rhs.valueType == Finite && value < rhs.value)));
                    }
                };

                struct MultiplicationRange
                {
                    Value minValue, maxValue;

                    MultiplicationRange():
                        minValue(Value::PlusInf),
                        maxValue(Value::MinusInf) {}

                    void multiply(Value value1, const Value& value2)
                    {
                        value1 *= value2;
                        if(value1 < minValue) minValue = value1;
                        if(maxValue < value1) maxValue = value1;
                    }
                };

                range<Value_t> result(Value_t(1), Value_t(1));
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    range<Value_t> item = CalculateResultBoundaries( tree.GetParam(a) );
                    if(!item.min.known && !item.max.known) return range<Value_t>(); // hopeless

                    Value minValue0 = result.min.known ? Value(result.min.val) : Value(Value::MinusInf);
                    Value maxValue0 = result.max.known ? Value(result.max.val) : Value(Value::PlusInf);
                    Value minValue1 = item.min.known ? Value(item.min.val) : Value(Value::MinusInf);
                    Value maxValue1 = item.max.known ? Value(item.max.val) : Value(Value::PlusInf);

                    MultiplicationRange range;
                    range.multiply(minValue0, minValue1);
                    range.multiply(minValue0, maxValue1);
                    range.multiply(maxValue0, minValue1);
                    range.multiply(maxValue0, maxValue1);

                    if(range.minValue.valueType == Value::Finite)
                        result.min.val = range.minValue.value;
                    else result.min.known = false;

                    if(range.maxValue.valueType == Value::Finite)
                        result.max.val = range.maxValue.value;
                    else result.max.known = false;

                    if(!result.min.known && !result.max.known) break; // hopeless
                }
                if(result.min.known && result.max.known
                && result.min.val > result.max.val) std::swap(result.min.val, result.max.val);
                return result;
            }
            case cMod:
            {
                /* TODO: The boundaries of modulo operator could be estimated better. */

                range<Value_t> x = CalculateResultBoundaries( tree.GetParam(0) );
                range<Value_t> y = CalculateResultBoundaries( tree.GetParam(1) );

                if(y.max.known)
                {
                    if(y.max.val >= Value_t(0))
                    {
                        if(!x.min.known || (x.min.val) < Value_t(0))
                            return range<Value_t>(-y.max.val, y.max.val);
                        else
                            return range<Value_t>(Value_t(0), y.max.val);
                    }
                    else
                    {
                        if(!x.max.known || (x.max.val) >= Value_t(0))
                            return range<Value_t>(y.max.val, -y.max.val);
                        else
                            return range<Value_t>(y.max.val, fp_const_negativezero<Value_t>());
                    }
                }
                else
                    return range<Value_t>();
            }
            case cPow:
            {
                if(tree.GetParam(1).IsImmed() && tree.GetParam(1).GetImmed() == Value_t(0))
                {
                    // Note: This makes 0^0 evaluate into 1.
                    return range<Value_t>(Value_t(1), Value_t(1)); // x^0 = 1
                }
                if(tree.GetParam(0).IsImmed() && tree.GetParam(0).GetImmed() == Value_t(0))
                {
                    // Note: This makes 0^0 evaluate into 0.
                    return range<Value_t>(Value_t(0), Value_t(0)); // 0^x = 0
                }
                if(tree.GetParam(0).IsImmed() && fp_equal(tree.GetParam(0).GetImmed(), Value_t(1)))
                {
                    return range<Value_t>(Value_t(1), Value_t(1)); // 1^x = 1
                }
                if(tree.GetParam(1).IsImmed()
                && tree.GetParam(1).GetImmed() > Value_t(0)
                && GetEvennessInfo(tree.GetParam(1)) == IsAlways)
                {
                    // x ^ even_int_const always produces a non-negative value.
                    Value_t exponent = tree.GetParam(1).GetImmed();
                    range<Value_t> tmp = CalculateResultBoundaries( tree.GetParam(0) );
                    range<Value_t> result;
                    result.min.known = true;
                    result.min.val = 0;
                    if(tmp.min.known && tmp.min.val >= Value_t(0))
                        result.min.val = fp_pow(tmp.min.val, exponent);
                    else if(tmp.max.known && tmp.max.val <= Value_t(0))
                        result.min.val = fp_pow(tmp.max.val, exponent);

                    result.max.known = false;
                    if(tmp.min.known && tmp.max.known)
                    {
                        result.max.known = true;
                        result.max.val     = fp_max(fp_abs(tmp.min.val), fp_abs(tmp.max.val));
                        result.max.val     = fp_pow(result.max.val, exponent);
                    }
                    return result;
                }

                range<Value_t> p0 = CalculateResultBoundaries( tree.GetParam(0) );
                range<Value_t> p1 = CalculateResultBoundaries( tree.GetParam(1) );
                TriTruthValue p0_positivity =
                    (p0.min.known && (p0.min.val) >= Value_t(0)) ? IsAlways
                  : (p0.max.known && (p0.max.val) < Value_t(0) ? IsNever
                    : Unknown);
                TriTruthValue p1_evenness = GetEvennessInfo(tree.GetParam(1));

                /* If param0 IsAlways, the return value is also IsAlways */
                /* If param1 is even, the return value is IsAlways */
                /* If param1 is odd, the return value is same as param0's */
                /* If param0 is negative and param1 is not integer,
                 * the return value is imaginary (assumed Unknown)
                 *
                 * Illustrated in this truth table:
                 *  P=positive, N=negative
                 *  E=even, O=odd, U=not integer
                 *  *=unknown, X=invalid (unknown), x=maybe invalid (unknown)
                 *
                 *   param1: PE PO P* NE NO N* PU NU *
                 * param0:
                 *   PE      P  P  P  P  P  P  P  P  P
                 *   PO      P  P  P  P  P  P  P  P  P
                 *   PU      P  P  P  P  P  P  P  P  P
                 *   P*      P  P  P  P  P  P  P  P  P
                 *   NE      P  N  *  P  N  *  X  X  x
                 *   NO      P  N  *  P  N  *  X  X  x
                 *   NU      P  N  *  P  N  *  X  X  x
                 *   N*      P  N  *  P  N  *  X  X  x
                 *   *       P  *  *  P  *  *  x  x  *
                 *
                 * Note: This also deals with the following opcodes:
                 *       cSqrt  (param0, PU) (x^0.5)
                 *       cRSqrt (param0, NU) (x^-0.5)
                 *       cExp   (PU, param1) (CONSTANT_E^x)
                 */
                TriTruthValue result_positivity = Unknown;
                switch(p0_positivity)
                {
                    case IsAlways:
                        // e.g.   5^x = positive.
                        result_positivity = IsAlways;
                        break;
                    case IsNever:
                    {
                        result_positivity = p1_evenness;
                        break;
                    }
                    default:
                        switch(p1_evenness)
                        {
                            case IsAlways:
                                // e.g. x^( 4) = positive
                                // e.g. x^(-4) = positive
                                result_positivity = IsAlways;
                                break;
                            case IsNever:
                                break;
                            case Unknown:
                            {
                                /* If p1 is const non-integer,
                                 * assume the result is positive
                                 * though it may be NaN instead.
                                 */
                                if(tree.GetParam(1).IsImmed()
                                && !isInteger(tree.GetParam(1).GetImmed())
                                && tree.GetParam(1).GetImmed() >= Value_t(0))
                                {
                                    result_positivity = IsAlways;
                                }
                                break;
                            }
                        }
                }
                switch(result_positivity)
                {
                    case IsAlways:
                    {
                        /* The result is always positive.
                         * Figure out whether we know the minimum value. */
                        Value_t min = Value_t(0);
                        if(p0.min.known && p1.min.known)
                        {
                            min = fp_pow(p0.min.val, p1.min.val);
                            if(p0.min.val < Value_t(0) && (!p1.max.known || p1.max.val >= Value_t(0)) && min >= Value_t(0))
                                min = Value_t(0);
                        }
                        if(p0.min.known && p0.min.val >= Value_t(0) && p0.max.known && p1.max.known)
                        {
                            Value_t max = fp_pow(p0.max.val, p1.max.val);
                            if(min > max) std::swap(min, max);
                            return range<Value_t>(min, max);
                        }
                        return range<Value_t>(min, false);
                    }
                    case IsNever:
                    {
                        /* The result is always negative.
                         * TODO: Figure out whether we know the maximum value.
                         */
                        return range<Value_t>(false, fp_const_negativezero<Value_t>());
                    }
                    default:
                    {
                        /* It can be negative or positive.
                         * We know nothing about the boundaries. */
                        break;
                    }
                }
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.cc and thus
             * they will never occur in the calling context for the
             * most of the parsing context. They may however occur
             * at the late phase, so we deal with them.
             */
            case cNeg:
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.set_neg();
                return m;
            }
            case cSub: // converted into cAdd(x, cNeg(y))
            {
                CodeTree<Value_t> tmp, tmp2;
                tmp2.SetOpcode(cNeg);
                tmp2.AddParam(tree.GetParam(1));
                tmp.SetOpcode(cAdd);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParamMove(tmp2);
                return CalculateResultBoundaries(tmp);
            }
            case cInv: // converted into cPow x -1
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParam(CodeTreeImmed(Value_t(-1)));
                return CalculateResultBoundaries(tmp);
            }
            case cDiv: // converted into cPow y -1
            {
                CodeTree<Value_t> tmp, tmp2;
                tmp2.SetOpcode(cInv);
                tmp2.AddParam(tree.GetParam(1));
                tmp.SetOpcode(cMul);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParamMove(tmp2);
                return CalculateResultBoundaries(tmp);
            }
            case cRad: // converted into cMul x CONSTANT_RD
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cMul);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParam(CodeTreeImmed(fp_const_rad_to_deg<Value_t>()));
                return CalculateResultBoundaries(tmp);
            }
            case cDeg: // converted into cMul x CONSTANT_DR
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cMul);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParam(CodeTreeImmed(fp_const_deg_to_rad<Value_t>()));
                return CalculateResultBoundaries(tmp);
            }
            case cSqr: // converted into cMul x x    or cPow x 2
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParam(CodeTreeImmed(Value_t(2)));
                return CalculateResultBoundaries(tmp);
            }
            case cExp: // converted into cPow CONSTANT_E x
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(CodeTreeImmed(fp_const_e<Value_t>()));
                tmp.AddParam(tree.GetParam(0));
                return CalculateResultBoundaries(tmp);
            }
            case cExp2: // converted into cPow 2 x
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(CodeTreeImmed(Value_t(2)));
                tmp.AddParam(tree.GetParam(0));
                return CalculateResultBoundaries(tmp);
            }
            case cCbrt: // converted into cPow x 0.33333333
            {
                // However, contrary to x^(1/3), this allows
                // negative values for x, and produces those
                // as well.
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                m.min.set(fp_cbrt);
                m.max.set(fp_cbrt);
                return m;
            }
            case cSqrt: // converted into cPow x 0.5
            {
                range<Value_t> m = CalculateResultBoundaries( tree.GetParam(0) );
                if(m.min.known) m.min.val = (m.min.val) < Value_t(0) ? 0 : fp_sqrt(m.min.val);
                if(m.max.known) m.max.val = (m.max.val) < Value_t(0) ? 0 : fp_sqrt(m.max.val);
                return m;
            }
            case cRSqrt: // converted into cPow x -0.5
            {
                CodeTree<Value_t> tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(tree.GetParam(0));
                tmp.AddParam(CodeTreeImmed( Value_t(-0.5) ));
                return CalculateResultBoundaries(tmp);
            }
            case cHypot: // converted into cSqrt(cAdd(cMul(x x), cMul(y y)))
            {
                CodeTree<Value_t> xsqr, ysqr, add, sqrt;
                xsqr.AddParam(tree.GetParam(0)); xsqr.AddParam(CodeTreeImmed( Value_t(2) ));
                ysqr.AddParam(tree.GetParam(1)); ysqr.AddParam(CodeTreeImmed( Value_t(2) ));
                xsqr.SetOpcode(cPow); ysqr.SetOpcode(cPow);
                add.AddParamMove(xsqr); add.AddParamMove(ysqr);
                add.SetOpcode(cAdd); sqrt.AddParamMove(add);
                sqrt.SetOpcode(cSqrt);
                return CalculateResultBoundaries(sqrt);
            }
            case cLog2by: // converted into cMul y CONSTANT_L2I (cLog x)
            {
                CodeTree<Value_t> tmp, tmp2;
                tmp2.SetOpcode(cLog2);
                tmp2.AddParam(tree.GetParam(0));
                tmp.SetOpcode(cMul);
                tmp.AddParamMove(tmp2);
                tmp.AddParam(tree.GetParam(1));
                return CalculateResultBoundaries(tmp);
            }
            case cCot: // converted into 1 / cTan
            {
                CodeTree<Value_t> tmp, tmp2;
                tmp2.SetOpcode(cTan);
                tmp2.AddParam(tree.GetParam(0));
                tmp.SetOpcode(cInv);
                tmp.AddParamMove(tmp2);
                return CalculateResultBoundaries(tmp);
            }
            case cSec: // converted into 1 / cCos
            {
                CodeTree<Value_t> tmp, tmp2;
                tmp2.SetOpcode(cCos);
                tmp2.AddParam(tree.GetParam(0));
                tmp.SetOpcode(cInv);
                tmp.AddParamMove(tmp2);
                return CalculateResultBoundaries(tmp);
            }
            case cCsc: // converted into 1 / cSin
            {
                CodeTree<Value_t> tmp, tmp2;
                tmp2.SetOpcode(cSin);
                tmp2.AddParam(tree.GetParam(0));
                tmp.SetOpcode(cInv);
                tmp.AddParamMove(tmp2);
                return CalculateResultBoundaries(tmp);
            }
            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.cc and thus
             * they will never occur in the calling context:
             */
                break; /* Should never occur */

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
            case VarBegin:
                break; /* Should never occur */

            /* Complex functions */
            case cArg:
            case cConj:
            case cImag:
            case cReal:
            case cPolar:
                break; /* Should never occur */

            /* Opcodes that are completely unpredictable */
            case cPCall:
                break;
            case cFCall:
                break; // Cannot deduce
            case cEval:
                break; // Cannot deduce
        }
        return range<Value_t>(); /* Cannot deduce */
    }

    template<typename Value_t>
    TriTruthValue GetIntegerInfo(const CodeTree<Value_t>& tree)
    {
        switch(tree.GetOpcode())
        {
            case cImmed:
                return isInteger(tree.GetImmed()) ? IsAlways : IsNever;
            case cFloor:
            case cCeil:
            case cTrunc:
            case cInt:
                return IsAlways;
            case cAnd:
            case cOr:
            case cNot:
            case cNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
                /* These operations always produce truth values (0 or 1) */
                return IsAlways; /* 0 and 1 are both integers */
            case cIf:
            {
                TriTruthValue a = GetIntegerInfo(tree.GetParam(1));
                TriTruthValue b = GetIntegerInfo(tree.GetParam(2));
                if(a == b) return a;
                return Unknown;
            }
            case cAdd:
            case cMul:
            {
                // It's integer if all the components are integer
                // Otherwise, unknown whether it's integer
                // A confirmed non-integer does not necessarily
                // mean the result isn't an integer, because:
                // 0.5 + 0.5 = 1.0; sqrt(2) * sqrt(2) = 2.0
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                    if(GetIntegerInfo(tree.GetParam(a)) != IsAlways)
                        return Unknown;
                return IsAlways;
            }
            default:
                break;
        }
        return Unknown; /* Don't know whether it's integer. */
    }

    template<typename Value_t>
    bool IsLogicalValue(const CodeTree<Value_t>& tree)
    {
        switch(tree.GetOpcode())
        {
            case cImmed:
                return fp_equal(tree.GetImmed(), Value_t(0))
                    || fp_equal(tree.GetImmed(), Value_t(1));
            case cAnd:
            case cOr:
            case cNot:
            case cNotNot:
            case cAbsAnd:
            case cAbsOr:
            case cAbsNot:
            case cAbsNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
                /* These operations always produce truth values (0 or 1) */
                return true;
            case cMul:
            {
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                    if(!IsLogicalValue(tree.GetParam(a)))
                        return false;
                return true;
            }
            case cIf:
            case cAbsIf:
            {
                return IsLogicalValue(tree.GetParam(1))
                    && IsLogicalValue(tree.GetParam(2));
            }
            default:
                break;
        }
        return false; // Not a logical value.
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template range<type> CalculateResultBoundaries(const CodeTree<type> &); \
    template bool IsLogicalValue(const CodeTree<type> &); \
    template TriTruthValue GetIntegerInfo(const CodeTree<type> &);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
