#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH

#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>

#include "hash.hh"
#include "../lib/autoptr.hh"

namespace FPoptimizer_Grammar
{
    struct Grammar;
}

namespace FPoptimizer_ByteCode
{
    template<typename Value_t>
    class ByteCodeSynth;
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    class CodeTree;

    template<typename Value_t>
    struct CodeTreeData;

    template<typename Value_t>
    class CodeTree
    {
        typedef FPOPT_autoptr<CodeTreeData<Value_t> > DataP;
        DataP data;

    public:
        CodeTree();
        ~CodeTree();

        struct OpcodeTag { };
        explicit CodeTree(FUNCTIONPARSERTYPES::OPCODE o, OpcodeTag); // produce an opcode
        struct FuncOpcodeTag { };
        explicit CodeTree(FUNCTIONPARSERTYPES::OPCODE o, unsigned f, FuncOpcodeTag);
        struct ImmedTag { };
        explicit CodeTree(const Value_t& v, ImmedTag); // produce an immed
#ifdef __GXX_EXPERIMENTAL_CXX0X__
        explicit CodeTree(Value_t&& v, ImmedTag); // produce an immed
#endif
        struct VarTag { };
        explicit CodeTree(unsigned varno, VarTag); // produce a var reference
        struct CloneTag { };
        explicit CodeTree(const CodeTree& b, CloneTag);

        /* Generates a CodeTree from the given bytecode */
        void GenerateFrom(
            const typename FunctionParserBase<Value_t>::Data& data,
            bool keep_powi = false);

        void GenerateFrom(
            const typename FunctionParserBase<Value_t>::Data& data,
            const std::vector<CodeTree>& var_trees,
            bool keep_powi = false);

        void SynthesizeByteCode(
            std::vector<unsigned>& byteCode,
            std::vector<Value_t>&   immed,
            size_t& stacktop_max);

        void SynthesizeByteCode(
            FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
            bool MustPopTemps=true) const;

        size_t SynthCommonSubExpressions(
            FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth) const;

        void SetParams(const std::vector<CodeTree>& RefParams);
        void SetParamsMove(std::vector<CodeTree>& RefParams);

        CodeTree GetUniqueRef();
        // ^use this when CodeTree tmp=x; tmp.CopyOnWrite(); does not do exactly what you want

#ifdef __GXX_EXPERIMENTAL_CXX0X__
        void SetParams(std::vector<CodeTree>&& RefParams);
#endif
        void SetParam(size_t which, const CodeTree& b);
        void SetParamMove(size_t which, CodeTree& b);
        void AddParam(const CodeTree& param);
        void AddParamMove(CodeTree& param);
        void AddParams(const std::vector<CodeTree>& RefParams);
        void AddParamsMove(std::vector<CodeTree>& RefParams);
        void AddParamsMove(std::vector<CodeTree>& RefParams, size_t replacing_slot);
        void DelParam(size_t index);
        void DelParams();

        void Become(const CodeTree& b);

        inline size_t GetParamCount() const { return GetParams().size(); }
        inline CodeTree& GetParam(size_t n) { return GetParams()[n]; }
        inline const CodeTree& GetParam(size_t n) const { return GetParams()[n]; }
        inline void SetOpcode(FUNCTIONPARSERTYPES::OPCODE o) { data->Opcode = o; }
        inline FUNCTIONPARSERTYPES::OPCODE GetOpcode() const { return data->Opcode; }
        inline FUNCTIONPARSERTYPES::fphash_t GetHash() const { return data->Hash; }
        inline const std::vector<CodeTree>& GetParams() const { return data->Params; }
        inline std::vector<CodeTree>& GetParams() { return data->Params; }
        inline size_t GetDepth() const         { return data->Depth; }
        inline const Value_t& GetImmed() const { return data->Value; }
        inline unsigned GetVar() const    { return data->Var_or_Funcno; }
        inline unsigned GetFuncNo() const { return data->Var_or_Funcno; }
        inline bool IsDefined() const { return GetOpcode() != FUNCTIONPARSERTYPES::cNop; }

        inline bool    IsImmed() const { return GetOpcode() == FUNCTIONPARSERTYPES::cImmed; }
        inline bool      IsVar() const { return GetOpcode() == FUNCTIONPARSERTYPES::VarBegin; }
        inline unsigned GetRefCount() const { return data->RefCount; }

        void ReplaceWithImmed(const Value_t& i);
        void Rehash(bool constantfolding = true);
        void Sort();
        inline void Mark_Incompletely_Hashed() { data->Depth=0; }
        inline bool Is_Incompletely_Hashed() const { return data->Depth == 0; }

        inline const FPoptimizer_Grammar::Grammar* GetOptimizedUsing() const
            { return data->OptimizedUsing; }
        inline void SetOptimizedUsing(const FPoptimizer_Grammar::Grammar* g)
            { data->OptimizedUsing = g; }

        bool RecreateInversionsAndNegations(bool prefer_base2 = false);
        void FixIncompleteHashes();

        void swap(CodeTree& b) { data.swap(b.data); }
        bool IsIdenticalTo(const CodeTree& b) const;
        void CopyOnWrite();
    };

    template<typename Value_t>
    struct CodeTreeData
    {
        int RefCount;

        /* Describing the codetree node */
        FUNCTIONPARSERTYPES::OPCODE Opcode;
        Value_t Value;          // In case of cImmed:   value of the immed
        unsigned Var_or_Funcno; // In case of VarBegin: variable number
                                // In case of cFCall or cPCall: function number

        // Parameters for the function
        //  These use the sign:
        //   For cAdd: operands to add together (0 to n)
        //             sign indicates that the value is negated before adding (0-x)
        //   For cMul: operands to multiply together (0 to n)
        //             sign indicates that the value is inverted before multiplying (1/x)
        //   For cAnd: operands to bitwise-and together (0 to n)
        //             sign indicates that the value is inverted before anding (!x)
        //   For cOr:  operands to bitwise-or together (0 to n)
        //             sign indicates that the value is inverted before orring (!x)
        //  These don't use the sign (sign is always false):
        //   For cMin: operands to select the minimum of
        //   For cMax: operands to select the maximum of
        //   For cImmed, not used
        //   For VarBegin, not used
        //   For cIf:  operand 1 = condition, operand 2 = yes-branch, operand 3 = no-branch
        //   For anything else: the parameters required by the operation/function
        std::vector<CodeTree<Value_t> > Params;

        /* Internal operation */
        FUNCTIONPARSERTYPES::fphash_t      Hash;
        size_t        Depth;
        const FPoptimizer_Grammar::Grammar* OptimizedUsing;

        CodeTreeData();
        CodeTreeData(const CodeTreeData& b);
        explicit CodeTreeData(FUNCTIONPARSERTYPES::OPCODE o);
        explicit CodeTreeData(FUNCTIONPARSERTYPES::OPCODE o, unsigned f);
        explicit CodeTreeData(const Value_t& i);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
        explicit CodeTreeData(Value_t&& i);
        CodeTreeData(CodeTreeData&& b);
#endif

        bool IsIdenticalTo(const CodeTreeData& b) const;
        void Sort();
        void Recalculate_Hash_NoRecursion();

    private:
        void operator=(const CodeTreeData& b);
    };

    /* Utility functions for creating different kind of CodeTrees */
    template<typename Value_t>
    static inline CodeTree<Value_t> CodeTreeImmed(const Value_t& i)
    {
        return CodeTree<Value_t> (i, typename CodeTree<Value_t>::ImmedTag());
    }

#ifdef __GXX_EXPERIMENTAL_CXX0X__
    template<typename Value_t>
    static inline CodeTree<Value_t> CodeTreeImmed(Value_t&& i)
    {
        return CodeTree<Value_t> (std::move(i),
                                  typename CodeTree<Value_t>::ImmedTag());
    }
#endif

    template<typename Value_t>
    static inline CodeTree<Value_t> CodeTreeOp(FUNCTIONPARSERTYPES::OPCODE opcode)
    {
        return CodeTree<Value_t> (opcode, typename CodeTree<Value_t>::OpcodeTag());
    }

    template<typename Value_t>
    static inline CodeTree<Value_t> CodeTreeFuncOp(FUNCTIONPARSERTYPES::OPCODE opcode, unsigned f)
    {
        return CodeTree<Value_t> (opcode, f, typename CodeTree<Value_t>::FuncOpcodeTag());
    }

    template<typename Value_t>
    static inline CodeTree<Value_t> CodeTreeVar(unsigned varno)
    {
        return CodeTree<Value_t> (varno, typename CodeTree<Value_t>::VarTag());
    }

    /* Debugging functions */
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
    template<typename Value_t>
    void DumpHashes(const CodeTree<Value_t>& tree, std::ostream& o = std::cout);

    template<typename Value_t>
    void DumpTree(const CodeTree<Value_t>& tree, std::ostream& o = std::cout);

    template<typename Value_t>
    void DumpTreeWithIndent(const CodeTree<Value_t>& tree, std::ostream& o = std::cout, const std::string& indent = "\\");
#endif
}

#endif

#endif
