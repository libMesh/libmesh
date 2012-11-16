#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>

#include "codetree.hh"

#ifndef FP_GENERATING_POWI_TABLE
enum { MAX_POWI_BYTECODE_LENGTH = 20 };
#else
enum { MAX_POWI_BYTECODE_LENGTH = 999 };
#endif
enum { MAX_MULI_BYTECODE_LENGTH = 3 };

namespace FPoptimizer_ByteCode
{
    template<typename Value_t>
    class ByteCodeSynth
    {
    public:
        ByteCodeSynth()
            : ByteCode(), Immed(), StackState(), StackTop(0), StackMax(0)
        {
            /* estimate the initial requirements as such */
            ByteCode.reserve(64);
            Immed.reserve(8);
            StackState.reserve(16);
        }

        void Pull(std::vector<unsigned>& bc,
                  std::vector<Value_t>&   imm,
                  size_t& StackTop_max)
        {
            /* The bitmask 0x80000000u was added to each non-opcode
             * value within ByteCode[] (opcode parameters) to prevent
             * them being interpreted as opcodes by fp_opcode_add.inc.
             * fparser uses cNop for the same purpose.
             */
            for(unsigned a=0; a<ByteCode.size(); ++a)
            {
                ByteCode[a] &= ~0x80000000u;
            }
            ByteCode.swap(bc);
            Immed.swap(imm);
            StackTop_max = StackMax;
        }

        size_t GetByteCodeSize() const { return ByteCode.size(); }
        size_t GetStackTop()     const { return StackTop; }

        void PushVar(unsigned varno)
        {
            ByteCode.push_back(varno);
            SetStackTop(StackTop+1);
        }

        void PushImmed(Value_t immed)
        {
            using namespace FUNCTIONPARSERTYPES;
            ByteCode.push_back(cImmed);
            Immed.push_back(immed);
            SetStackTop(StackTop+1);
        }

        void StackTopIs(const FPoptimizer_CodeTree::CodeTree<Value_t>& tree, int offset = 0)
        {
            if((int)StackTop > offset)
            {
                StackState[StackTop-1-offset].first = true;
                StackState[StackTop-1-offset].second = tree;
            }
        }

        bool IsStackTop(const FPoptimizer_CodeTree::CodeTree<Value_t>& tree, int offset = 0) const
        {
            return (int)StackTop > offset
               && StackState[StackTop-1-offset].first
               && StackState[StackTop-1-offset].second.IsIdenticalTo(tree);
        }

        inline void EatNParams(unsigned eat_count)
        {
            StackTop -= eat_count;
        }

        void ProducedNParams(unsigned produce_count)
        {
            SetStackTop(StackTop + produce_count);
        }

        void DoPopNMov(size_t targetpos, size_t srcpos)
        {
            using namespace FUNCTIONPARSERTYPES;
            ByteCode.push_back(cPopNMov);
            ByteCode.push_back( 0x80000000u | (unsigned) targetpos);
            ByteCode.push_back( 0x80000000u | (unsigned) srcpos);

            SetStackTop(srcpos+1);
            StackState[targetpos] = StackState[srcpos];
            SetStackTop(targetpos+1);
        }

        void DoDup(size_t src_pos)
        {
            using namespace FUNCTIONPARSERTYPES;
            if(src_pos == StackTop-1)
            {
                ByteCode.push_back(cDup);
            }
            else
            {
                ByteCode.push_back(cFetch);
                ByteCode.push_back( 0x80000000u | (unsigned) src_pos);
            }
            SetStackTop(StackTop + 1);
            StackState[StackTop-1] = StackState[src_pos];
        }

#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
        template<int/*defer*/>
        void Dump()
        {
            std::ostream& o = std::cout;
            o << "Stack state now(" << StackTop << "):\n";
            for(size_t a=0; a<StackTop; ++a)
            {
                o << a << ": ";
                if(StackState[a].first)
                {
                    const FPoptimizer_CodeTree::CodeTree<Value_t>
                        & tree = StackState[a].second;
                    o << '[' << std::hex << (void*)(&tree.GetParams())
                             << std::dec
                             << ',' << tree.GetRefCount()
                             << ']';
                    DumpTree(tree, o);
                }
                else
                    o << "?";
                o << "\n";
            }
            o << std::flush;
        }
#endif

        size_t FindPos(const FPoptimizer_CodeTree::CodeTree<Value_t>& tree) const
        {
            for(size_t a=StackTop; a-->0; )
                if(StackState[a].first && StackState[a].second.IsIdenticalTo(tree))
                    return a;
            return ~size_t(0);
        }

        bool Find(const FPoptimizer_CodeTree::CodeTree<Value_t>& tree) const
        {
            return FindPos(tree) != ~size_t(0);
        }

        bool FindAndDup(const FPoptimizer_CodeTree::CodeTree<Value_t>& tree)
        {
            size_t pos = FindPos(tree);
            if(pos != ~size_t(0))
            {
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Found duplicate at [" << pos <<"]: ";
                DumpTree(tree);
                std::cout << " -- issuing cDup or cFetch\n";
            #endif
                DoDup(pos);
                return true;
            }
            return false;
        }

        struct IfData
        {
            size_t ofs;
        };

        void SynthIfStep1(IfData& ifdata, FUNCTIONPARSERTYPES::OPCODE op)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // the If condition was popped.

            ifdata.ofs = ByteCode.size();
            ByteCode.push_back(op);
            ByteCode.push_back(0x80000000u); // code index
            ByteCode.push_back(0x80000000u); // Immed index
        }
        void SynthIfStep2(IfData& ifdata)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // ignore the pushed then-branch result.

            ByteCode[ifdata.ofs+1] = 0x80000000u | unsigned( ByteCode.size()+2 );
            ByteCode[ifdata.ofs+2] = 0x80000000u | unsigned( Immed.size()      );

            ifdata.ofs = ByteCode.size();
            ByteCode.push_back(cJump);
            ByteCode.push_back(0x80000000u); // code index
            ByteCode.push_back(0x80000000u); // Immed index
        }
        void SynthIfStep3(IfData& ifdata)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // ignore the pushed else-branch result.

            ByteCode.back() |= 0x80000000u;
            // ^Necessary for guarding against if(x,1,2)+1 being changed
            //  into if(x,1,3) by fp_opcode_add.inc

            ByteCode[ifdata.ofs+1] = 0x80000000u | unsigned( ByteCode.size()-1 );
            ByteCode[ifdata.ofs+2] = 0x80000000u | unsigned( Immed.size()      );

            SetStackTop(StackTop+1); // one or the other was pushed.

            /* Threading jumps:
             * If there are any cJumps that point
             * to the cJump instruction we just changed,
             * change them to point to this target as well.
             * This screws up PrintByteCode() majorly.
             */
            for(size_t a=0; a<ifdata.ofs; ++a)
            {
                if(ByteCode[a]   == cJump
                && ByteCode[a+1] == (0x80000000u | (ifdata.ofs-1)))
                {
                    ByteCode[a+1] = 0x80000000u | unsigned( ByteCode.size()-1 );
                    ByteCode[a+2] = 0x80000000u | unsigned( Immed.size()      );
                }
                switch(ByteCode[a])
                {
                    case cAbsIf:
                    case cIf:
                    case cJump:
                    case cPopNMov: a += 2; break;
                    case cFCall:
                    case cPCall:
                    case cFetch: a += 1; break;
                    default: break;
                }
            }
        }

    protected:
        void SetStackTop(size_t value)
        {
            StackTop = value;
            if(StackTop > StackMax)
            {
                StackMax = StackTop;
                StackState.resize(StackMax);
            }
        }

    protected:
        std::vector<unsigned> ByteCode;
        std::vector<Value_t>   Immed;

        std::vector<
            std::pair<bool/*known*/,
                      FPoptimizer_CodeTree::CodeTree<Value_t>/*tree*/>
                   > StackState;
        size_t StackTop;
        size_t StackMax;
    private:
        void incStackPtr()
        {
            if(StackTop+2 > StackMax) StackState.resize(StackMax=StackTop+2);
        }

        template<bool IsIntType, bool IsComplexType>
        struct Specializer { };
    public:
        void AddOperation(unsigned opcode, unsigned eat_count, unsigned produce_count = 1)
        {
            EatNParams(eat_count);
            AddFunctionOpcode(opcode);
            ProducedNParams(produce_count);
        }

        void AddFunctionOpcode(unsigned opcode, Specializer<false,false>);
        void AddFunctionOpcode(unsigned opcode, Specializer<false,true>);
        void AddFunctionOpcode(unsigned opcode, Specializer<true,false>);
        void AddFunctionOpcode(unsigned opcode, Specializer<true,true>);
        inline void AddFunctionOpcode(unsigned opcode)
        {
            AddFunctionOpcode
                (opcode,
                 Specializer< bool(FUNCTIONPARSERTYPES::IsIntType<Value_t>::result),
                              bool(FUNCTIONPARSERTYPES::IsComplexType<Value_t>::result)
                           > ()
                         );
        }
    };

    template<typename Value_t>
    struct SequenceOpCode;
    template<typename Value_t>
    struct SequenceOpcodes
    {
        /* Multiplication implemented with adds */
        static const SequenceOpCode<Value_t> AddSequence;
        /* Exponentiation implemented with muls */
        static const SequenceOpCode<Value_t> MulSequence;
    };

    /* Generate a sequence that multiplies or exponentifies the
     * last operand in the stack by the given constant integer
     * amount (positive or negative).
     */
    template<typename Value_t>
    void AssembleSequence(
        long count,
        const SequenceOpCode<Value_t>& sequencing,
        ByteCodeSynth<Value_t>& synth);
}

#endif
