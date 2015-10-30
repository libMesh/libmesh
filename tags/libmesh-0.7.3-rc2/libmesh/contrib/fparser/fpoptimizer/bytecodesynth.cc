#include "bytecodesynth.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include "opcodename.hh"
#include "codetree.hh"

using namespace FUNCTIONPARSERTYPES;

namespace FPoptimizer_ByteCode
{
    template<typename Value_t>
    struct SequenceOpCode
    {
        Value_t basevalue;
        unsigned op_flip;
        unsigned op_normal, op_normal_flip;
        unsigned op_inverse, op_inverse_flip;
    };

    template<typename Value_t>
    const SequenceOpCode<Value_t>
          SequenceOpcodes<Value_t>::AddSequence = { Value_t(0), cNeg, cAdd, cAdd, cSub, cRSub };

    template<typename Value_t>
    const SequenceOpCode<Value_t>
          SequenceOpcodes<Value_t>::MulSequence = { Value_t(1), cInv, cMul, cMul, cDiv, cRDiv };

    /*******/
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode ByteCode
#define mImmed Immed

    template<typename Value_t>
    void ByteCodeSynth<Value_t>::AddFunctionOpcode(unsigned opcode,
                                                   Specializer<false,false>)
    {
        int mStackPtr=0;
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
    }

    template<typename Value_t>
    void ByteCodeSynth<Value_t>::AddFunctionOpcode(unsigned opcode,
                                                   Specializer<true,false>)
    {
        int mStackPtr=0;
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
    }

#ifdef FP_SUPPORT_COMPLEX_NUMBERS
    template<typename Value_t>
    void ByteCodeSynth<Value_t>::AddFunctionOpcode(unsigned opcode,
                                                   Specializer<false,true>)
    {
        int mStackPtr=0;
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
    }

    template<typename Value_t>
    void ByteCodeSynth<Value_t>::AddFunctionOpcode(unsigned opcode,
                                                   Specializer<true,true>)
    {
        int mStackPtr=0;
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
    }
#endif

#undef findName
#undef mImmed
#undef mByteCode
#undef mData
#undef TryCompilePowi
    /*******/
}

using namespace FPoptimizer_ByteCode;

#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
namespace FPoptimizer_ByteCode
{
    #ifndef FP_GENERATING_POWI_TABLE
    extern const
    unsigned char powi_table[POWI_TABLE_SIZE];
    const
    #endif
    unsigned char powi_table[POWI_TABLE_SIZE] =
    {
          0,   1,   1,   1,   2,   1,   2,   1, /*   0 -   7 */
          4,   1,   2,   1,   4,   1,   2, 131, /*   8 -  15 */
          8,   1,   2,   1,   4,   1,   2,   1, /*  16 -  23 */
          8, 133,   2, 131,   4,   1,  15,   1, /*  24 -  31 */
         16,   1,   2,   1,   4,   1,   2, 131, /*  32 -  39 */
          8,   1,   2,   1,   4, 133,   2,   1, /*  40 -  47 */
         16,   1,  25, 131,   4,   1,  27,   5, /*  48 -  55 */
          8,   3,   2,   1,  30,   1,  31,   3, /*  56 -  63 */
         32,   1,   2,   1,   4,   1,   2,   1, /*  64 -  71 */
          8,   1,   2, 131,   4,   1,  39,   1, /*  72 -  79 */
         16, 137,   2,   1,   4, 133,   2, 131, /*  80 -  87 */
          8,   1,  45, 135,   4,  31,   2,   5, /*  88 -  95 */
         32,   1,   2, 131,  50,   1,  51,   1, /*  96 - 103 */
          8,   3,   2,   1,  54,   1,  55,   3, /* 104 - 111 */
         16,   1,  57, 133,   4, 137,   2, 135, /* 112 - 119 */
         60,   1,  61,   3,  62, 133,  63,   1, /* 120 - 127 */
        130,   1,   2,   1, 130,   1,   2, 131, /* 128 - 135 */
        130,   1,   2,   1, 130,   1,   2, 139, /* 136 - 143 */
        130,   1,   2, 131, 130,   1,  30,   1, /* 144 - 151 */
        130, 137,   2,  31, 130,   1,   2, 131, /* 152 - 159 */
        130,   1, 130,   1, 130, 133,   2,   1, /* 160 - 167 */
        130,   1, 130,   1,   2,   1, 130, 133, /* 168 - 175 */
        130,   1,   2,   1, 130,   1,   2,  61, /* 176 - 183 */
        130, 133,  62, 139, 130, 137, 130,   1, /* 184 - 191 */
        130,   1,   2, 131, 130,   1, 130,   1, /* 192 - 199 */
        130,   1,   2,   1, 130,   1,   2, 131, /* 200 - 207 */
        130,   1, 130,   1, 130, 131,   2, 133, /* 208 - 215 */
        130,   1,   2, 131, 130, 141, 130,   1, /* 216 - 223 */
        130, 133,   2,   1, 130,   1,   5, 135, /* 224 - 231 */
        130,   1, 130,   1,   2, 131, 130,   1, /* 232 - 239 */
        130,   1,   2, 131, 130, 133, 130, 141, /* 240 - 247 */
        130, 131, 130,   1, 130,   1,   2, 131  /* 248 - 255 */
    }; /* as in gcc, but custom-optimized for stack calculation */
}
static const int POWI_CACHE_SIZE = 256;

#define FPO(x) /**/
//#define FPO(x) x
//#include <stdio.h>


namespace
{
    class PowiCache
    {
    private:
        int cache[POWI_CACHE_SIZE];
        int cache_needed[POWI_CACHE_SIZE];

    public:
        PowiCache()
            : cache(), cache_needed() /* Assume we have no factors in the cache */
        {
            /* Decide which factors we would need multiple times.
             * Output:
             *   cache[]        = these factors were generated
             *   cache_needed[] = number of times these factors were desired
             */
            cache[1] = 1; // We have this value already.
        }

        bool Plan_Add(long value, int count)
        {
            if(value >= POWI_CACHE_SIZE) return false;
            //FPO(fprintf(stderr, "%ld will be needed %d times more\n", count, need_count));
            cache_needed[value] += count;
            return cache[value] != 0;
        }

        void Plan_Has(long value)
        {
            if(value < POWI_CACHE_SIZE)
                cache[value] = 1; // This value has been generated
        }

        void Start(size_t value1_pos)
        {
            for(int n=2; n<POWI_CACHE_SIZE; ++n)
                cache[n] = -1; /* Stack location for each component */

            Remember(1, value1_pos);

            DumpContents();
        }

        int Find(long value) const
        {
            if(value < POWI_CACHE_SIZE)
            {
                if(cache[value] >= 0)
                {
                    // found from the cache
                    FPO(fprintf(stderr, "* I found %ld from cache (%u,%d)\n",
                        value, (unsigned)cache[value], cache_needed[value]));
                    return cache[value];
                }
            }
            return -1;
        }

        void Remember(long value, size_t stackpos)
        {
            if(value >= POWI_CACHE_SIZE) return;

            FPO(fprintf(stderr, "* Remembering that %ld can be found at %u (%d uses remain)\n",
                value, (unsigned)stackpos, cache_needed[value]));
            cache[value] = (int) stackpos;
        }

        void DumpContents() const
        {
            FPO(for(int a=1; a<POWI_CACHE_SIZE; ++a)
                if(cache[a] >= 0 || cache_needed[a] > 0)
                {
                    fprintf(stderr, "== cache: sp=%d, val=%d, needs=%d\n",
                        cache[a], a, cache_needed[a]);
                })
        }

        int UseGetNeeded(long value)
        {
            if(value >= 0 && value < POWI_CACHE_SIZE)
                return --cache_needed[value];
            return 0;
        }
    };

    template<typename Value_t>
    size_t AssembleSequence_Subdivide(
        long count,
        PowiCache& cache,
        const SequenceOpCode<Value_t>& sequencing,
        ByteCodeSynth<Value_t>& synth);

    template<typename Value_t>
    void Subdivide_Combine(
        size_t apos, long aval,
        size_t bpos, long bval,
        PowiCache& cache,

        unsigned cumulation_opcode,
        unsigned cimulation_opcode_flip,

        ByteCodeSynth<Value_t>& synth);

    void PlanNtimesCache
        (long value,
         PowiCache& cache,
         int need_count,
         int recursioncount=0)
    {
        if(value < 1) return;

    #ifdef FP_GENERATING_POWI_TABLE
        if(recursioncount > 32) throw false;
    #endif

        if(cache.Plan_Add(value, need_count)) return;

        long half = 1;
        if(value < POWI_TABLE_SIZE)
        {
            half = powi_table[value];
            if(half & 128)
            {
                half &= 127;
                if(half & 64)
                    half = -(half & 63) - 1;

                FPO(fprintf(stderr, "value=%ld, half=%ld, otherhalf=%ld\n", value,half,value/half));

                PlanNtimesCache(half,      cache, 1, recursioncount+1);
                cache.Plan_Has(half);
                return;
            }
            else if(half & 64)
            {
                half = -(half & 63) - 1;
            }
        }
        else if(value & 1)
            half = value & ((1 << POWI_WINDOW_SIZE) - 1); // that is, value & 7
        else
            half = value / 2;

        long otherhalf = value-half;
        if(half > otherhalf || half<0) std::swap(half,otherhalf);

        FPO(fprintf(stderr, "value=%ld, half=%ld, otherhalf=%ld\n", value,half,otherhalf));

        if(half == otherhalf)
        {
            PlanNtimesCache(half,      cache, 2, recursioncount+1);
        }
        else
        {
            PlanNtimesCache(half,      cache, 1, recursioncount+1);
            PlanNtimesCache(otherhalf>0?otherhalf:-otherhalf,
                                       cache, 1, recursioncount+1);
        }
        cache.Plan_Has(value);
    }

    template<typename Value_t>
    size_t AssembleSequence_Subdivide(
        long value,
        PowiCache& cache,
        const SequenceOpCode<Value_t>& sequencing,
        ByteCodeSynth<Value_t>& synth)
    {
        int cachepos = cache.Find(value);
        if(cachepos >= 0)
        {
            // found from the cache
            return cachepos;
        }

        long half = 1;
        if(value < POWI_TABLE_SIZE)
        {
            half = powi_table[value];
            if(half & 128)
            {
                half &= 127;
                if(half & 64)
                    half = -(half & 63) - 1;

                FPO(fprintf(stderr, "* I want %ld, my plan is %ld * %ld\n", value, half, value/half));
                size_t half_pos = AssembleSequence_Subdivide(half, cache, sequencing, synth);
                if(cache.UseGetNeeded(half) > 0
                || half_pos != synth.GetStackTop()-1)
                {
                    synth.DoDup(half_pos);
                    cache.Remember(half, synth.GetStackTop()-1);
                }
                AssembleSequence(value/half, sequencing, synth);
                size_t stackpos = synth.GetStackTop()-1;
                cache.Remember(value, stackpos);
                cache.DumpContents();
                return stackpos;
            }
            else if(half & 64)
            {
                half = -(half & 63) - 1;
            }
        }
        else if(value & 1)
            half = value & ((1 << POWI_WINDOW_SIZE) - 1); // that is, value & 7
        else
            half = value / 2;

        long otherhalf = value-half;
        if(half > otherhalf || half<0) std::swap(half,otherhalf);

        FPO(fprintf(stderr, "* I want %ld, my plan is %ld + %ld\n", value, half, value-half));

        if(half == otherhalf)
        {
            size_t half_pos = AssembleSequence_Subdivide(half, cache, sequencing, synth);

            // self-cumulate the subdivide result
            Subdivide_Combine(half_pos,half, half_pos,half, cache,
                sequencing.op_normal, sequencing.op_normal_flip,
                synth);
        }
        else
        {
            long part1 = half;
            long part2 = otherhalf>0?otherhalf:-otherhalf;

            size_t part1_pos = AssembleSequence_Subdivide(part1, cache, sequencing, synth);
            size_t part2_pos = AssembleSequence_Subdivide(part2, cache, sequencing, synth);

            FPO(fprintf(stderr, "Subdivide(%ld: %ld, %ld)\n", value, half, otherhalf));

            Subdivide_Combine(part1_pos,part1, part2_pos,part2, cache,
                otherhalf>0 ? sequencing.op_normal      : sequencing.op_inverse,
                otherhalf>0 ? sequencing.op_normal_flip : sequencing.op_inverse_flip,
                synth);
        }

        size_t stackpos = synth.GetStackTop()-1;
        cache.Remember(value, stackpos);
        cache.DumpContents();
        return stackpos;
    }

    template<typename Value_t>
    void Subdivide_Combine(
        size_t apos, long aval,
        size_t bpos, long bval,
        PowiCache& cache,
        unsigned cumulation_opcode,
        unsigned cumulation_opcode_flip,
        ByteCodeSynth<Value_t>& synth)
    {
        /*FPO(fprintf(stderr, "== making result for (sp=%u, val=%d, needs=%d) and (sp=%u, val=%d, needs=%d), stacktop=%u\n",
            (unsigned)apos, aval, aval>=0 ? cache_needed[aval] : -1,
            (unsigned)bpos, bval, bval>=0 ? cache_needed[bval] : -1,
            (unsigned)synth.GetStackTop()));*/

        // Figure out whether we can trample a and b
        int a_needed = cache.UseGetNeeded(aval);
        int b_needed = cache.UseGetNeeded(bval);

        bool flipped = false;

        #define DUP_BOTH() do { \
            if(apos < bpos) { size_t tmp=apos; apos=bpos; bpos=tmp; flipped=!flipped; } \
            FPO(fprintf(stderr, "-> dup(%u) dup(%u) op\n", (unsigned)apos, (unsigned)bpos)); \
            synth.DoDup(apos); \
            synth.DoDup(apos==bpos ? synth.GetStackTop()-1 : bpos); } while(0)
        #define DUP_ONE(p) do { \
            FPO(fprintf(stderr, "-> dup(%u) op\n", (unsigned)p)); \
            synth.DoDup(p); \
        } while(0)

        if(a_needed > 0)
        {
            if(b_needed > 0)
            {
                // If they must both be preserved, make duplicates
                // First push the one that is at the larger stack
                // address. This increases the odds of possibly using cDup.
                DUP_BOTH();

                //SCENARIO 1:
                // Input:  x B A x x
                // Temp:   x B A x x A B
                // Output: x B A x x R
                //SCENARIO 2:
                // Input:  x A B x x
                // Temp:   x A B x x B A
                // Output: x A B x x R
            }
            else
            {
                // A must be preserved, but B can be trampled over

                // SCENARIO 1:
                //  Input:  x B x x A
                //   Temp:  x B x x A A B   (dup both, later first)
                //  Output: x B x x A R
                // SCENARIO 2:
                //  Input:  x A x x B
                //   Temp:  x A x x B A
                //  Output: x A x x R       -- only commutative cases
                // SCENARIO 3:
                //  Input:  x x x B A
                //   Temp:  x x x B A A B   (dup both, later first)
                //  Output: x x x B A R
                // SCENARIO 4:
                //  Input:  x x x A B
                //   Temp:  x x x A B A     -- only commutative cases
                //  Output: x x x A R
                // SCENARIO 5:
                //  Input:  x A B x x
                //   Temp:  x A B x x A B   (dup both, later first)
                //  Output: x A B x x R

                // if B is not at the top, dup both.
                if(bpos != synth.GetStackTop()-1)
                    DUP_BOTH();    // dup both
                else
                {
                    DUP_ONE(apos); // just dup A
                    flipped=!flipped;
                }
            }
        }
        else if(b_needed > 0)
        {
            // B must be preserved, but A can be trampled over
            // This is a mirror image of the a_needed>0 case, so I'll cut the chase
            if(apos != synth.GetStackTop()-1)
                DUP_BOTH();
            else
                DUP_ONE(bpos);
        }
        else
        {
            // Both can be trampled over.
            // SCENARIO 1:
            //  Input:  x B x x A
            //   Temp:  x B x x A B
            //  Output: x B x x R
            // SCENARIO 2:
            //  Input:  x A x x B
            //   Temp:  x A x x B A
            //  Output: x A x x R       -- only commutative cases
            // SCENARIO 3:
            //  Input:  x x x B A
            //  Output: x x x R         -- only commutative cases
            // SCENARIO 4:
            //  Input:  x x x A B
            //  Output: x x x R
            // SCENARIO 5:
            //  Input:  x A B x x
            //   Temp:  x A B x x A B   (dup both, later first)
            //  Output: x A B x x R
            // SCENARIO 6:
            //  Input:  x x x C
            //   Temp:  x x x C C   (c is both A and B)
            //  Output: x x x R

            if(apos == bpos && apos == synth.GetStackTop()-1)
                DUP_ONE(apos); // scenario 6
            else if(apos == synth.GetStackTop()-1 && bpos == synth.GetStackTop()-2)
            {
                FPO(fprintf(stderr, "-> op\n")); // scenario 3
                flipped=!flipped;
            }
            else if(apos == synth.GetStackTop()-2 && bpos == synth.GetStackTop()-1)
                FPO(fprintf(stderr, "-> op\n")); // scenario 4
            else if(apos == synth.GetStackTop()-1)
                DUP_ONE(bpos); // scenario 1
            else if(bpos == synth.GetStackTop()-1)
            {
                DUP_ONE(apos); // scenario 2
                flipped=!flipped;
            }
            else
                DUP_BOTH(); // scenario 5
        }
        // Add them together.
        synth.AddOperation(flipped ? cumulation_opcode_flip : cumulation_opcode, 2);
    }

    template<typename Value_t>
    void LightWeight(
        long count,
        const SequenceOpCode<Value_t>& sequencing,
        ByteCodeSynth<Value_t>& synth)
    {
        while(count < 256)
        {
            int half = FPoptimizer_ByteCode::powi_table[count];
            if(half & 128)
            {
                half &= 127;
                LightWeight(half,       sequencing, synth);
                count /= half;
            }
            else break;
        }
        if(count == 1) return;
        if(!(count & 1))
        {
            synth.AddOperation(cSqr, 1);
            LightWeight(count/2, sequencing, synth);
        }
        else
        {
            synth.DoDup(synth.GetStackTop()-1);
            LightWeight(count-1, sequencing, synth);
            synth.AddOperation(cMul, 2);
        }
    }
}

namespace FPoptimizer_ByteCode
{
    template<typename Value_t>
    void AssembleSequence(
        long count,
        const SequenceOpCode<Value_t>& sequencing,
        ByteCodeSynth<Value_t>& synth)
    {
        if(count == 0)
            synth.PushImmed(sequencing.basevalue);
        else
        {
            bool needs_flip = false;
            if(count < 0)
            {
                needs_flip = true;
                count = -count;
            }

            if(false)
                LightWeight(count,sequencing,synth);
            else if(count > 1)
            {
                /* To prevent calculating the same factors over and over again,
                 * we use a cache. */
                PowiCache cache;
                PlanNtimesCache(count, cache, 1);

                size_t stacktop_desired = synth.GetStackTop();

                cache.Start( synth.GetStackTop()-1 );

                FPO(fprintf(stderr, "Calculating result for %ld...\n", count));
                size_t res_stackpos = AssembleSequence_Subdivide(
                    count, cache, sequencing,
                    synth);

                size_t n_excess = synth.GetStackTop() - stacktop_desired;
                if(n_excess > 0 || res_stackpos != stacktop_desired-1)
                {
                    // Remove the cache values
                    synth.DoPopNMov(stacktop_desired-1, res_stackpos);
                }
            }

            if(needs_flip)
                synth.AddOperation(sequencing.op_flip, 1);
        }
    }
}

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_ByteCode
{
#define FP_INSTANTIATE(type) \
    template struct SequenceOpcodes<type>; \
    template void ByteCodeSynth<type>::AddFunctionOpcode(unsigned); \
    template void ByteCodeSynth<type>::AddFunctionOpcode(unsigned, \
                 Specializer< bool(FUNCTIONPARSERTYPES::IsIntType<type>::result), \
                              bool(FUNCTIONPARSERTYPES::IsComplexType<type>::result) \
                           > ); \
    template void AssembleSequence( \
        long count, \
        const SequenceOpCode<type>& sequencing, \
        ByteCodeSynth<type>& synth);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
