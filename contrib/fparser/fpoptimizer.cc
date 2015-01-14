/***************************************************************************\
|* Function Parser for C++ v4.5.1                                          *|
|*-------------------------------------------------------------------------*|
|* Function optimizer                                                      *|
|*-------------------------------------------------------------------------*|
|* Copyright: Joel Yliluoma                                                *|
|*                                                                         *|
|* This library is distributed under the terms of the                      *|
|* GNU Lesser General Public License version 3.                            *|
|* (See lgpl.txt and gpl.txt for the license text.)                        *|
\***************************************************************************/

/* NOTE:
 This file contains generated code (from the optimizer sources) and is
 not intended to be modified by hand. If you want to modify the optimizer,
 download the development version of the library.
*/

#include "fpconfig.hh"
#ifdef FP_SUPPORT_OPTIMIZER
#include "fparser.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"
#line 1 "lib/crc32.hh"
/* crc32 */

#ifdef _MSC_VER

 typedef unsigned int crc32_t;

#else

 #include <stdint.h>
 typedef uint_least32_t crc32_t;

#endif

namespace crc32
{
    enum { startvalue = 0xFFFFFFFFUL, poly = 0xEDB88320UL };

    /* This code constructs the CRC32 table at compile-time,
     * avoiding the need for a huge explicitly written table of magical numbers. */
    template<crc32_t crc> // One byte of a CRC32 (eight bits):
    struct b8
    {
        enum { b1 = (crc & 1) ? (poly ^ (crc >> 1)) : (crc >> 1),
               b2 = (b1  & 1) ? (poly ^ (b1  >> 1)) : (b1  >> 1),
               b3 = (b2  & 1) ? (poly ^ (b2  >> 1)) : (b2  >> 1),
               b4 = (b3  & 1) ? (poly ^ (b3  >> 1)) : (b3  >> 1),
               b5 = (b4  & 1) ? (poly ^ (b4  >> 1)) : (b4  >> 1),
               b6 = (b5  & 1) ? (poly ^ (b5  >> 1)) : (b5  >> 1),
               b7 = (b6  & 1) ? (poly ^ (b6  >> 1)) : (b6  >> 1),
               res= (b7  & 1) ? (poly ^ (b7  >> 1)) : (b7  >> 1) };
    };
    inline crc32_t update(crc32_t crc, unsigned/* char */b) // __attribute__((pure))
    {
        // Four values of the table
        #define B4(n) b8<n>::res,b8<n+1>::res,b8<n+2>::res,b8<n+3>::res
        // Sixteen values of the table
        #define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
        // The whole table, index by steps of 16
        static const crc32_t table[256] =
        { R(0x00),R(0x10),R(0x20),R(0x30), R(0x40),R(0x50),R(0x60),R(0x70),
          R(0x80),R(0x90),R(0xA0),R(0xB0), R(0xC0),R(0xD0),R(0xE0),R(0xF0) };
        #undef R
        #undef B4
        return ((crc >> 8) /* & 0x00FFFFFF*/) ^ table[/*(unsigned char)*/(crc^b)&0xFF];
    }
    inline crc32_t calc_upd(crc32_t c, const unsigned char* buf, size_t size)
    {
        crc32_t value = c;
        for(size_t p=0; p<size; ++p) value = update(value, buf[p]);
        return value;
    }
    inline crc32_t calc(const unsigned char* buf, size_t size)
    {
        return calc_upd(startvalue, buf, size);
    }
}

#line 1 "lib/autoptr.hh"
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH

template<typename Ref>
class FPOPT_autoptr
{
public:
    FPOPT_autoptr()                   : p(0)   { }
    FPOPT_autoptr(Ref*        b) : p(b)   { Birth(); }
    FPOPT_autoptr(const FPOPT_autoptr& b) : p(b.p) { Birth(); }

    inline Ref& operator* () const { return *p; }
    inline Ref* operator->() const { return p; }

    FPOPT_autoptr& operator= (Ref*        b) { Set(b); return *this; }
    FPOPT_autoptr& operator= (const FPOPT_autoptr& b) { Set(b.p); return *this; }
#ifdef FP_SUPPORT_CXX11_MOVE
    FPOPT_autoptr(FPOPT_autoptr&& b)      : p(b.p) { b.p = 0; }
    FPOPT_autoptr& operator= (FPOPT_autoptr&& b) { if(p != b.p) { Forget(); p=b.p; b.p=0; }
                                                   return *this; }
#endif

    ~FPOPT_autoptr() { Forget(); }

    bool isNull() const { return p == 0; }

    void UnsafeSetP(Ref* newp) { p = newp; }
    void swap(FPOPT_autoptr<Ref>& b) { Ref* tmp=p; p=b.p; b.p=tmp; }

private:
    inline static void Have(Ref* p2);
    inline void Forget();
    inline void Birth();
    inline void Set(Ref* p2);
private:
    Ref* p;
};

//
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Forget()
{
    if(!p) return;
    p->RefCount -= 1;
    if(!p->RefCount) delete p;
    //assert(p->RefCount >= 0);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Have(Ref* p2)
{
    if(p2) ++(p2->RefCount);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Birth()
{
    Have(p);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Set(Ref* p2)
{
    Have(p2);
    Forget();
    p = p2;
}

#endif

#line 1 "lib/functional.hh"
#include <utility>

struct Compare2ndRev
{
    template<typename T>
    inline bool operator() (const T& a, const T& b) const
    {
        return a.second > b.second;
    }
};

struct Compare1st
{
    template<typename T1, typename T2>
    inline bool operator() (const std::pair<T1,T2>& a,
                            const std::pair<T1,T2>& b) const
    {
        return a.first < b.first;
    }

    template<typename T1, typename T2>
    inline bool operator() (const std::pair<T1,T2>& a, T1 b) const
    {
        return a.first < b;
    }

    template<typename T1, typename T2>
    inline bool operator() (T1 a, const std::pair<T1,T2>& b) const
    {
        return a < b.first;
    }
};

#line 1 "fpoptimizer/hash.hh"
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH

#ifdef _MSC_VER

typedef unsigned long long fphash_value_t;
#define FPHASH_CONST(x) x##ULL

#else

#include <stdint.h>
typedef uint_fast64_t fphash_value_t;
#define FPHASH_CONST(x) x##ULL

#endif

namespace FUNCTIONPARSERTYPES
{
    struct fphash_t
    {
        fphash_value_t hash1, hash2;

        fphash_t() : hash1(0), hash2(0) { }
        fphash_t(const fphash_value_t& a,
                 const fphash_value_t& b) : hash1(a), hash2(b) { }

        bool operator==(const fphash_t& rhs) const
        { return hash1 == rhs.hash1 && hash2 == rhs.hash2; }

        bool operator!=(const fphash_t& rhs) const
        { return hash1 != rhs.hash1 || hash2 != rhs.hash2; }

        bool operator<(const fphash_t& rhs) const
        { return hash1 != rhs.hash1 ? hash1 < rhs.hash1 : hash2 < rhs.hash2; }
    };
}

#endif

#line 1 "fpoptimizer/codetree.hh"
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH

// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fpaux.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>

// line removed for fpoptimizer.cc: #include "hash.hh"
// line removed for fpoptimizer.cc: #include "../lib/autoptr.hh"

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
#ifdef FP_SUPPORT_CXX11_MOVE
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

#ifdef FP_SUPPORT_CXX11_MOVE
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
#ifdef FP_SUPPORT_CXX11_MOVE
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

#ifdef FP_SUPPORT_CXX11_MOVE
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

#line 1 "fpoptimizer/grammar.hh"
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH

#include <iostream>

// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    class CodeTree;
}

namespace FPoptimizer_Grammar
{
    enum ImmedConstraint_Value
    {
        ValueMask = 0x07,
        Value_AnyNum     = 0x0, // any value
        Value_EvenInt    = 0x1, // any even integer (0,2,4, etc)
        Value_OddInt     = 0x2, // any odd integer (1,3, etc)
        Value_IsInteger  = 0x3, // any integer-value (excludes e.g. 0.2)
        Value_NonInteger = 0x4, // any non-integer (excludes e.g. 1 or 5)
        Value_Logical    = 0x5  // a result of cNot,cNotNot,cAnd,cOr or comparators
    };
    enum ImmedConstraint_Sign
    {
        SignMask  = 0x18,
        Sign_AnySign     = 0x00, // - or +
        Sign_Positive    = 0x08, // positive only
        Sign_Negative    = 0x10, // negative only
        Sign_NoIdea      = 0x18  // where sign cannot be guessed
    };
    enum ImmedConstraint_Oneness
    {
        OnenessMask   = 0x60,
        Oneness_Any      = 0x00,
        Oneness_One      = 0x20, // +1 or -1
        Oneness_NotOne   = 0x40  // anything but +1 or -1
    };
    enum ImmedConstraint_Constness
    {
        ConstnessMask = 0x180,
        Constness_Any    = 0x00,
        Constness_Const  = 0x80,
        Constness_NotConst=0x100
    };
    enum Modulo_Mode
    {
        Modulo_None    = 0,
        Modulo_Radians = 1
    };
    enum Situation_Flags
    {
        LogicalContextOnly = 0x01,
        NotForIntegers     = 0x02,
        OnlyForIntegers    = 0x04,
        OnlyForComplex     = 0x08,
        NotForComplex      = 0x10
    };

    /* The param_opcode field of the ParamSpec has the following
     * possible values (from enum SpecialOpcode):
     *   NumConstant:
     *      this describes a specific constant value (constvalue)
     *      that must be matched / synthesized.
     *   ParamHolder:
     *      this describes any node
     *      that must be matched / synthesized.
     *      "index" is the ID of the NamedHolder:
     *      In matching, all NamedHolders having the same ID
     *      must match the identical node.
     *      In synthesizing, the node matched by
     *      a NamedHolder with this ID must be synthesized.
     *    SubFunction:
     *      this describes a subtree
     *      that must be matched / synthesized.
     *      The subtree is described in subfunc_opcode,param_begin..+param_count.
     *      If the type is GroupFunction, the tree is expected
     *      to yield a constant value which is tested.
     */
    enum SpecialOpcode
    {
        NumConstant,        // Holds a particular value (syntax-time constant)
        ParamHolder,        // Holds a particular named param
        SubFunction         // Holds an opcode and the params
    };

    enum ParamMatchingType
    {
        PositionalParams, // this set of params in this order
        SelectedParams,   // this set of params in any order
        AnyParams,        // these params are included
        GroupFunction     // this function represents a constant value
    };

    enum RuleType
    {
        ProduceNewTree, // replace self with the first (and only) from replaced_param
        ReplaceParams   // replace indicate params with replaced_params
    };

#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif

    /* A ParamSpec object describes
     * either a parameter (leaf, node) that must be matched,
     * or a parameter (leaf, node) that must be synthesized.
     */
    typedef std::pair<SpecialOpcode, const void*> ParamSpec;

    template<typename Value_t>
    ParamSpec ParamSpec_Extract(unsigned paramlist, unsigned index);

    template<typename Value_t>
    bool ParamSpec_Compare(const void* a, const void* b, SpecialOpcode type);

    unsigned ParamSpec_GetDepCode(const ParamSpec& b);

    struct ParamSpec_ParamHolder
    {
        unsigned index       : 8; // holder ID
        unsigned constraints : 9; // constraints
        unsigned depcode     :15;
    } PACKED_GRAMMAR_ATTRIBUTE;

    template<typename Value_t>
    struct ParamSpec_NumConstant
    {
        Value_t     constvalue; // the value
        unsigned    modulo;     // modulo mode
    };// PACKED_GRAMMAR_ATTRIBUTE;

    struct ParamSpec_SubFunctionData
    {
        /* Expected parameters (leaves) of the tree: */
        unsigned param_count         : 2;
        unsigned param_list          : 30;
        /* The opcode that the tree must have when SubFunction */
        FUNCTIONPARSERTYPES::OPCODE subfunc_opcode : 8;

        /* When matching, type describes the method of matching.
         *
         *               Sample input tree:      (cOr 2 3)  (cOr 2 4) (cOr 3 2) (cOr 4 2 3) (cOr 2)
         * Possible methods:
         *    PositionalParams, e.g. (cOr [2 3]):  match     no match  no match  no match   no match
         *      The nodes described here are
         *      to be matched, in this order.
         *    SelectedParams,   e.g. (cOr {2 3}):  match     no match   match    no match   no match
         *      The nodes described here are
         *      to be matched, in any order.
         *    AnyParams,        e.g. (cOr 2 3  ):  match     no match   match     match     no match
         *      At least the nodes described here
         *      are to be matched, in any order.
         * When synthesizing, the type is ignored.
         */
        ParamMatchingType match_type : 3; /* When SubFunction */
        /* Note: match_type needs 2, but we specify 3 because
         * otherwise Microsoft VC++ borks things up
         * as it interprets the value as signed.
         */
        /* Optional restholder index for capturing the rest of parameters (0=not used)
         * Only valid when match_type = AnyParams
         */
        unsigned restholder_index : 5;
    } PACKED_GRAMMAR_ATTRIBUTE; // size: 2+30+6+2+8=48 bits=6 bytes

    struct ParamSpec_SubFunction
    {
        ParamSpec_SubFunctionData data;
        unsigned constraints : 9; // constraints
        unsigned depcode     : 7;
    } PACKED_GRAMMAR_ATTRIBUTE; // 8 bytes

    /* Theoretical minimal sizes in each param_opcode cases:
     * Assume param_opcode needs 3 bits.
     *    NumConstant:   3 + 64              (or 3+4 if just use index to clist[])
     *    ParamHolder:   3 + 7 + 2           (7 for constraints, 2 for immed index)
     *    SubFunction:   3 + 7 + 2 + 2 + 3*9 = 41
     */

    /* A rule describes a pattern for matching
     * and the method how to reconstruct the
     * matched node(s) in the tree.
     */
    struct Rule
    {
        /* If the rule matched, this field describes how to perform
         * the replacement.
         *   When type==ProduceNewTree,
         *       the source tree is replaced entirely with
         *       the new tree described at repl_param_begin[0].
         *   When type==ReplaceParams,
         *       the matching leaves in the source tree are removed
         *       and new leaves are constructedfrom the trees
         *       described at repl_param_begin[0..repl_param_count].
         *       Other leaves remain intact.
         */
        RuleType  ruletype         : 2;
        unsigned  situation_flags  : 5;

        /* The replacement parameters (if NewTree, begin[0] represents the new tree) */
        unsigned  repl_param_count : 2+9; /* Assumed to be 1 when type == ProduceNewTree */
        unsigned  repl_param_list  : 30;

        /* The function that we must match. Always a SubFunction. */
        ParamSpec_SubFunctionData match_tree;
    } PACKED_GRAMMAR_ATTRIBUTE; // size: 2+1+13+2+30 + 48 = 96 bits = 12 bytes

    /* Grammar is a set of rules for tree substitutions. */
    struct Grammar
    {
        /* The rules of this grammar */
        unsigned rule_count;
        unsigned short rule_list[999]; // maximum limit...
        /* Note: Changing the limit has no effect to performance of
         * fparser. The limit is only actually used within grammar_parser.
         * A too low limit causes a memory corruption during the parse.
         * A too high limit just may cause inconvenience.
         * The actual grammar items linked to fparser are optimized for size,
         * and the size of the Grammar object may be considerably smaller
         * than what is indicated by this prototype.
         */
    };

    extern "C" {
        extern const Rule      grammar_rules[];
        /* 
        extern const Grammar   grammar_optimize_round1;
        extern const Grammar   grammar_optimize_round2;
        extern const Grammar   grammar_optimize_round3;
        extern const Grammar   grammar_optimize_round4;
        extern const Grammar   grammar_optimize_recreate;
        extern const Grammar   grammar_optimize_shortcut_logical_evaluation;
        extern const Grammar   grammar_optimize_nonshortcut_logical_evaluation;
        extern const Grammar   grammar_optimize_ignore_if_sideeffects;
        extern const Grammar   grammar_optimize_abslogical;
        extern const Grammar   grammar_optimize_base2_expand;
S */
    }

    template<typename Value_t>
    void DumpParam(const ParamSpec& p, std::ostream& o = std::cout);

    template<typename Value_t>
    void DumpParams(unsigned paramlist, unsigned count, std::ostream& o = std::cout);
}

#endif

#line 1 "fpoptimizer/consts.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fpaux.hh"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/*
#define CONSTANT_L10B  0.3010299956639811952137 // log10(2)
#define CONSTANT_L10BI 3.3219280948873623478703 // 1/log10(2)
#define CONSTANT_LB10  CONSTANT_L10BI          // log2(10)
#define CONSTANT_LB10I CONSTANT_L10B           // 1/log2(10)
*/

#define CONSTANT_POS_INF     HUGE_VAL  // positive infinity, from math.h
#define CONSTANT_NEG_INF   (-HUGE_VAL) // negative infinity

namespace FUNCTIONPARSERTYPES
{
    template<typename Value_t>
    inline Value_t fp_const_pihalf() // CONSTANT_PIHALF
    {
        return fp_const_pi<Value_t>() * Value_t(0.5);
    }
    template<typename Value_t>
    inline Value_t fp_const_twopi() // CONSTANT_TWOPI
    {
        Value_t result( fp_const_pi<Value_t>() );
        result += result;
        return result;
    }
    template<typename Value_t>
    inline Value_t fp_const_twoe() // CONSTANT_2E
    {
        Value_t result( fp_const_e<Value_t>() );
        result += result;
        return result;
    }
    template<typename Value_t>
    inline Value_t fp_const_twoeinv() // CONSTANT_2EI
    {
        Value_t result( fp_const_einv<Value_t>() );
        result += result;
        return result;
    }

    template<typename Value_t>
    inline Value_t fp_const_negativezero()
    {
        return -Epsilon<Value_t>::value;
    }
}

#line 1 "fpoptimizer/optimize.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "grammar.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>
#include <iostream>

//#define DEBUG_SUBSTITUTIONS

namespace FPoptimizer_Optimize
{
    using namespace FPoptimizer_Grammar;
    using namespace FPoptimizer_CodeTree;
    using namespace FUNCTIONPARSERTYPES;

    /* This struct collects information regarding the matching process so far */
    template<typename Value_t>
    class MatchInfo
    {
    public:
        std::vector<std::pair<bool,std::vector<CodeTree<Value_t> >
                             > > restholder_matches;
        std::vector<CodeTree<Value_t> > paramholder_matches;
        std::vector<unsigned> matched_params;
    public:
        MatchInfo(): restholder_matches(), paramholder_matches(), matched_params() {}
    public:
        /* These functions save data from matching */
        bool SaveOrTestRestHolder(
            unsigned restholder_index,
            const std::vector<CodeTree<Value_t> >& treelist)
        {
            if(restholder_matches.size() <= restholder_index)
            {
                restholder_matches.resize(restholder_index+1);
                restholder_matches[restholder_index].first  = true;
                restholder_matches[restholder_index].second = treelist;
                return true;
            }
            if(restholder_matches[restholder_index].first == false)
            {
                restholder_matches[restholder_index].first  = true;
                restholder_matches[restholder_index].second = treelist;
                return true;
            }
            const std::vector<CodeTree<Value_t> >& found =
                restholder_matches[restholder_index].second;
            if(treelist.size() != found.size())
                return false;
            for(size_t a=0; a<treelist.size(); ++a)
                if(!treelist[a].IsIdenticalTo(found[a]))
                    return false;
            return true;
        }

        void SaveRestHolder(
            unsigned restholder_index,
            std::vector<CodeTree<Value_t> >& treelist)
        {
            if(restholder_matches.size() <= restholder_index)
                restholder_matches.resize(restholder_index+1);
            restholder_matches[restholder_index].first = true;
            restholder_matches[restholder_index].second.swap(treelist);
        }

        bool SaveOrTestParamHolder(
            unsigned paramholder_index,
            const CodeTree<Value_t>& treeptr)
        {
            if(paramholder_matches.size() <= paramholder_index)
            {
                paramholder_matches.reserve(paramholder_index+1);
                paramholder_matches.resize(paramholder_index);
                paramholder_matches.push_back(treeptr);
                return true;
            }
            if(!paramholder_matches[paramholder_index].IsDefined())
            {
                paramholder_matches[paramholder_index] = treeptr;
                return true;
            }
            return treeptr.IsIdenticalTo(paramholder_matches[paramholder_index]);
        }

        void SaveMatchedParamIndex(unsigned index)
        {
            matched_params.push_back(index);
        }

        /* These functions retrieve the data from matching
         * for use when synthesizing the resulting tree.
         */
        const CodeTree<Value_t>& GetParamHolderValueIfFound( unsigned paramholder_index ) const
        {
            static const CodeTree<Value_t> dummytree;
            if(paramholder_matches.size() <= paramholder_index)
                return dummytree;
            return paramholder_matches[paramholder_index];
        }

        const CodeTree<Value_t>& GetParamHolderValue( unsigned paramholder_index ) const
            { return paramholder_matches[paramholder_index]; }

        bool HasRestHolder(unsigned restholder_index) const
            { return restholder_matches.size() > restholder_index
                  && restholder_matches[restholder_index].first == true; }

        const std::vector<CodeTree<Value_t> >& GetRestHolderValues( unsigned restholder_index ) const
        {
            static const std::vector<CodeTree<Value_t> > empty_result;
            if(restholder_matches.size() <= restholder_index)
                return empty_result;
            return restholder_matches[restholder_index].second;
        }

        const std::vector<unsigned>& GetMatchedParamIndexes() const
            { return matched_params; }

        /* */
        void swap(MatchInfo<Value_t>& b)
        {
            restholder_matches.swap(b.restholder_matches);
            paramholder_matches.swap(b.paramholder_matches);
            matched_params.swap(b.matched_params);
        }
        MatchInfo<Value_t>& operator=(const MatchInfo<Value_t>& b)
        {
            restholder_matches = b.restholder_matches;
            paramholder_matches = b.paramholder_matches;
            matched_params = b.matched_params;
            return *this;
        }
    };

    class MatchPositionSpecBase;

    // For iterating through match candidates
    typedef FPOPT_autoptr<MatchPositionSpecBase> MatchPositionSpecBaseP;

    class MatchPositionSpecBase
    {
    public:
        int RefCount;
    public:
        MatchPositionSpecBase() : RefCount(0) { }
        virtual ~MatchPositionSpecBase() { }
    };
    struct MatchResultType
    {
        bool found;
        MatchPositionSpecBaseP specs;

        MatchResultType(bool f) : found(f), specs() { }
        MatchResultType(bool f,
                        const MatchPositionSpecBaseP& s) : found(f), specs(s) { }
    };

    /* Synthesize the given grammatic rule's replacement into the codetree */
    template<typename Value_t>
    void SynthesizeRule(
        const Rule& rule,
        CodeTree<Value_t>& tree,
        MatchInfo<Value_t>& info);

    /* Test the given parameter to a given CodeTree */
    template<typename Value_t>
    MatchResultType TestParam(
        const ParamSpec& parampair,
        const CodeTree<Value_t> & tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo<Value_t>& info);

    /* Test the list of parameters to a given CodeTree */
    template<typename Value_t>
    MatchResultType TestParams(
        const ParamSpec_SubFunctionData& model_tree,
        const CodeTree<Value_t> & tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo<Value_t>& info,
        bool TopLevel);

    template<typename Value_t>
    bool ApplyGrammar(const Grammar& grammar,
                      FPoptimizer_CodeTree::CodeTree<Value_t> & tree,
                      bool from_logical_context = false);

    template<typename Value_t>
    void ApplyGrammars(FPoptimizer_CodeTree::CodeTree<Value_t>& tree);

    template<typename Value_t>
    bool IsLogisticallyPlausibleParamsMatch(
        const ParamSpec_SubFunctionData& params,
        const CodeTree<Value_t>& tree);
}

namespace FPoptimizer_Grammar
{
    template<typename Value_t>
    void DumpMatch(const Rule& rule,
                   const FPoptimizer_CodeTree::CodeTree<Value_t> & tree,
                   const FPoptimizer_Optimize::MatchInfo<Value_t>& info,
                   bool DidMatch,
                   std::ostream& o = std::cout);

    template<typename Value_t>
    void DumpMatch(const Rule& rule,
                   const FPoptimizer_CodeTree::CodeTree<Value_t> & tree,
                   const FPoptimizer_Optimize::MatchInfo<Value_t>& info,
                   const char* whydump,
                   std::ostream& o = std::cout);
}

#endif

#line 1 "fpoptimizer/opcodename.hh"
// line removed for fpoptimizer.cc: #include "grammar.hh"
#include <string>

const std::string FP_GetOpcodeName(FPoptimizer_Grammar::SpecialOpcode opcode, bool pad=false);
const std::string FP_GetOpcodeName(FUNCTIONPARSERTYPES::OPCODE opcode,        bool pad=false);

#line 1 "fpoptimizer/opcodename.cc"
#include <string>
#include <sstream>
#include <assert.h>

#include <iostream>

// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

// line removed for fpoptimizer.cc: #include "grammar.hh"
// line removed for fpoptimizer.cc: #include "opcodename.hh"

using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

const std::string FP_GetOpcodeName(FPoptimizer_Grammar::SpecialOpcode opcode, bool pad)
{
#if 1
    /* Symbolic meanings for the opcodes? */
    const char* p = 0;
    switch( opcode )
    {
        case NumConstant:   p = "NumConstant"; break;
        case ParamHolder:   p = "ParamHolder"; break;
        case SubFunction:   p = "SubFunction"; break;
    }
    std::ostringstream tmp;
    //if(!p) std::cerr << "o=" << opcode << "\n";
    assert(p);
    tmp << p;
    if(pad) while(tmp.str().size() < 12) tmp << ' ';
    return tmp.str();
#else
    /* Just numeric meanings */
    std::ostringstream tmp;
    tmp << opcode;
    if(pad) while(tmp.str().size() < 5) tmp << ' ';
    return tmp.str();
#endif
}

const std::string FP_GetOpcodeName(FUNCTIONPARSERTYPES::OPCODE opcode,        bool pad)
{
#if 1
    /* Symbolic meanings for the opcodes? */
    const char* p = 0;
    switch(opcode)
    {
        case cAbs: p = "cAbs"; break;
        case cAcos: p = "cAcos"; break;
        case cAcosh: p = "cAcosh"; break;
        case cArg: p = "cArg"; break;
        case cAsin: p = "cAsin"; break;
        case cAsinh: p = "cAsinh"; break;
        case cAtan: p = "cAtan"; break;
        case cAtan2: p = "cAtan2"; break;
        case cAtanh: p = "cAtanh"; break;
        case cCbrt: p = "cCbrt"; break;
        case cCeil: p = "cCeil"; break;
        case cConj: p = "cConj"; break;
        case cCos: p = "cCos"; break;
        case cCosh: p = "cCosh"; break;
        case cCot: p = "cCot"; break;
        case cCsc: p = "cCsc"; break;
        case cExp: p = "cExp"; break;
        case cExp2: p = "cExp2"; break;
        case cFloor: p = "cFloor"; break;
        case cHypot: p = "cHypot"; break;
        case cIf: p = "cIf"; break;
        case cImag: p = "cImag"; break;
        case cInt: p = "cInt"; break;
        case cLog: p = "cLog"; break;
        case cLog2: p = "cLog2"; break;
        case cLog10: p = "cLog10"; break;
        case cMax: p = "cMax"; break;
        case cMin: p = "cMin"; break;
        case cPolar: p = "cPolar"; break;
        case cPow: p = "cPow"; break;
        case cReal: p = "cReal"; break;
        case cSec: p = "cSec"; break;
        case cSin: p = "cSin"; break;
        case cSinh: p = "cSinh"; break;
        case cSqrt: p = "cSqrt"; break;
        case cTan: p = "cTan"; break;
        case cTanh: p = "cTanh"; break;
        case cTrunc: p = "cTrunc"; break;
        case cImmed: p = "cImmed"; break;
        case cJump: p = "cJump"; break;
        case cNeg: p = "cNeg"; break;
        case cAdd: p = "cAdd"; break;
        case cSub: p = "cSub"; break;
        case cMul: p = "cMul"; break;
        case cDiv: p = "cDiv"; break;
        case cMod: p = "cMod"; break;
        case cEqual: p = "cEqual"; break;
        case cNEqual: p = "cNEqual"; break;
        case cLess: p = "cLess"; break;
        case cLessOrEq: p = "cLessOrEq"; break;
        case cGreater: p = "cGreater"; break;
        case cGreaterOrEq: p = "cGreaterOrEq"; break;
        case cNot: p = "cNot"; break;
        case cAnd: p = "cAnd"; break;
        case cOr: p = "cOr"; break;
        case cDeg: p = "cDeg"; break;
        case cRad: p = "cRad"; break;
        case cFCall: p = "cFCall"; break;
        case cPCall: p = "cPCall"; break;
#ifdef FP_SUPPORT_OPTIMIZER
        case cFetch: p = "cFetch"; break;
        case cPopNMov: p = "cPopNMov"; break;
        case cLog2by: p = "cLog2by"; break;
        case cNop: p = "cNop"; break;
#endif
        case cSinCos: p = "cSinCos"; break;
        case cSinhCosh: p = "cSinhCosh"; break;
        case cAbsNot: p = "cAbsNot"; break;
        case cAbsNotNot: p = "cAbsNotNot"; break;
        case cAbsAnd: p = "cAbsAnd"; break;
        case cAbsOr: p = "cAbsOr"; break;
        case cAbsIf: p = "cAbsIf"; break;
        case cDup: p = "cDup"; break;
        case cInv: p = "cInv"; break;
        case cSqr: p = "cSqr"; break;
        case cRDiv: p = "cRDiv"; break;
        case cRSub: p = "cRSub"; break;
        case cNotNot: p = "cNotNot"; break;
        case cRSqrt: p = "cRSqrt"; break;
        case VarBegin: p = "VarBegin"; break;
    }
    std::ostringstream tmp;
    //if(!p) std::cerr << "o=" << opcode << "\n";
    assert(p);
    tmp << p;
    if(pad) while(tmp.str().size() < 12) tmp << ' ';
    return tmp.str();
#else
    /* Just numeric meanings */
    std::ostringstream tmp;
    tmp << opcode;
    if(pad) while(tmp.str().size() < 5) tmp << ' ';
    return tmp.str();
#endif
}

#line 1 "fpoptimizer/bytecodesynth.hh"
// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>

// line removed for fpoptimizer.cc: #include "codetree.hh"

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

#line 1 "fpoptimizer/bytecodesynth.cc"
// line removed for fpoptimizer.cc: #include "bytecodesynth.hh"

#ifdef FP_SUPPORT_OPTIMIZER

// line removed for fpoptimizer.cc: #include "opcodename.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"

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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif

# include "extrasrc/fp_opcode_add.inc"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif

# include "extrasrc/fp_opcode_add.inc"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif

# include "extrasrc/fp_opcode_add.inc"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

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
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif

# include "extrasrc/fp_opcode_add.inc"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

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

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
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
 */

#endif

#line 1 "fpoptimizer/valuerange.hh"
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH

// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fpaux.hh"

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

#line 1 "fpoptimizer/rangeestimation.hh"
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH

// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "valuerange.hh"

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

#line 1 "fpoptimizer/constantfolding.hh"
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH

// line removed for fpoptimizer.cc: #include "codetree.hh"

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void ConstantFolding(CodeTree<Value_t>& tree);
}

#endif

#line 1 "fpoptimizer/logic_boolgroups.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"

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

#line 1 "fpoptimizer/logic_collections.hh"
#include <vector>
#include <map>
#include <algorithm>

// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "../lib/functional.hh"

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

#line 1 "fpoptimizer/logic_ifoperations.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /**************************************/
    /* IF OPERATIONS                      */
    /**************************************/

    template<typename Value_t>
    bool ConstantFolding_IfOperations(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cIf || tree.GetOpcode() == cAbsIf);

        // If the If() condition begins with a cNot,
        // remove the cNot and swap the branches.
        for(;;)
        {
            if(tree.GetParam(0).GetOpcode() == cNot) // TEST 20/ifnot
            {
                tree.SetOpcode(cIf);
                tree.GetParam(0).Become( tree.GetParam(0).GetParam(0) );
                tree.GetParam(1).swap(tree.GetParam(2));
            }
            else if(tree.GetParam(0).GetOpcode() == cAbsNot) // TEST 20/ifabsnot
            {
                tree.SetOpcode(cAbsIf);
                tree.GetParam(0).Become( tree.GetParam(0).GetParam(0) );
                tree.GetParam(1).swap(tree.GetParam(2));
            }
            else break;
        }

        // If the sub-expression evaluates to approx. zero, yield param3.
        // If the sub-expression evaluates to approx. nonzero, yield param2.
        switch(GetLogicalValue(tree.GetParam(0), tree.GetOpcode()==cAbsIf))
        { // TEST 20/ifconst
            case IsAlways: // true
                tree.Become(tree.GetParam(1));
                return true; // rerun optimization (opcode changed)
            case IsNever: // false
                tree.Become(tree.GetParam(2));
                return true; // rerun optimization (opcode changed)
            case Unknown: default: ;
        }

        if(tree.GetParam(0).GetOpcode() == cIf
        || tree.GetParam(0).GetOpcode() == cAbsIf)
        {
            // TEST 20/ififconst
            //     if(if(x, a,b), c,d)
            //  -> if(x, if(a, c,d), if(b, c,d))
            // when either a or b is constantly true/false
            CodeTree<Value_t> cond = tree.GetParam(0);
            CodeTree<Value_t> truth_a;
            truth_a.SetOpcode(cond.GetOpcode() == cIf ? cNotNot : cAbsNotNot);
            truth_a.AddParam(cond.GetParam(1));
            ConstantFolding(truth_a);
            CodeTree<Value_t> truth_b;
            truth_b.SetOpcode(cond.GetOpcode() == cIf ? cNotNot : cAbsNotNot);
            truth_b.AddParam(cond.GetParam(2));
            ConstantFolding(truth_b);
            if(truth_a.IsImmed() || truth_b.IsImmed())
            {
                CodeTree<Value_t> then_tree;
                then_tree.SetOpcode(cond.GetOpcode());
                then_tree.AddParam(cond.GetParam(1));
                then_tree.AddParam(tree.GetParam(1));
                then_tree.AddParam(tree.GetParam(2));
                then_tree.Rehash();
                CodeTree<Value_t> else_tree;
                else_tree.SetOpcode(cond.GetOpcode());
                else_tree.AddParam(cond.GetParam(2));
                else_tree.AddParam(tree.GetParam(1));
                else_tree.AddParam(tree.GetParam(2));
                else_tree.Rehash();
                tree.SetOpcode(cond.GetOpcode());
                tree.SetParam(0, cond.GetParam(0));
                tree.SetParamMove(1, then_tree);
                tree.SetParamMove(2, else_tree);
                return true; // rerun cIf optimization
            }
        }
        if(tree.GetParam(1).GetOpcode() == tree.GetParam(2).GetOpcode()
        && (tree.GetParam(1).GetOpcode() == cIf
         || tree.GetParam(1).GetOpcode() == cAbsIf))
        {
            CodeTree<Value_t>& leaf1 = tree.GetParam(1);
            CodeTree<Value_t>& leaf2 = tree.GetParam(2);
            if(leaf1.GetParam(0).IsIdenticalTo(leaf2.GetParam(0))
            && (leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(1))
             || leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(2))))
            {
            // TEST 20/ifmerge
            //     if(x, if(y,a,b), if(y,c,d))
            // ->  if(y, if(x,a,c), if(x,b,d))
            // when either a,c are identical or b,d are identical
                CodeTree<Value_t> then_tree;
                then_tree.SetOpcode(tree.GetOpcode());
                then_tree.AddParam(tree.GetParam(0));
                then_tree.AddParam(leaf1.GetParam(1));
                then_tree.AddParam(leaf2.GetParam(1));
                then_tree.Rehash();
                CodeTree<Value_t> else_tree;
                else_tree.SetOpcode(tree.GetOpcode());
                else_tree.AddParam(tree.GetParam(0));
                else_tree.AddParam(leaf1.GetParam(2));
                else_tree.AddParam(leaf2.GetParam(2));
                else_tree.Rehash();
                tree.SetOpcode(leaf1.GetOpcode());
                tree.SetParam(0, leaf1.GetParam(0));
                tree.SetParamMove(1, then_tree);
                tree.SetParamMove(2, else_tree);
                return true; // rerun cIf optimization
            // cIf [x (cIf [y a z]) (cIf [y z b])] : (cXor x y) z (cIf[x a b])
            // ^ if only we had cXor opcode.
            }
            if(leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(1))
            && leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(2)))
            {
                // TEST 20/ifmerge2
                //    if(x, if(y,a,b), if(z,a,b))
                // -> if( if(x, y,z), a,b)
                CodeTree<Value_t> cond_tree;
                cond_tree.SetOpcode(tree.GetOpcode());
                cond_tree.AddParamMove(tree.GetParam(0));
                cond_tree.AddParam(leaf1.GetParam(0));
                cond_tree.AddParam(leaf2.GetParam(0));
                cond_tree.Rehash();
                tree.SetOpcode(leaf1.GetOpcode());
                tree.SetParamMove(0, cond_tree);
                tree.SetParam(2, leaf1.GetParam(2));
                tree.SetParam(1, leaf1.GetParam(1));
                return true; // rerun cIf optimization
            }
            if(leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(2))
            && leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(1)))
            {
                // TEST 20/ifmerge2b
                //    if(x, if(y,a,b), if(z,b,a))
                // -> if( if(x, y,!z), a,b)
                CodeTree<Value_t> not_tree;
                not_tree.SetOpcode(leaf2.GetOpcode() == cIf ? cNot : cAbsNot);
                not_tree.AddParam(leaf2.GetParam(0));
                not_tree.Rehash();
                CodeTree<Value_t> cond_tree;
                cond_tree.SetOpcode(tree.GetOpcode());
                cond_tree.AddParamMove(tree.GetParam(0));
                cond_tree.AddParam(leaf1.GetParam(0));
                cond_tree.AddParamMove(not_tree);
                cond_tree.Rehash();
                tree.SetOpcode(leaf1.GetOpcode());
                tree.SetParamMove(0, cond_tree);
                tree.SetParam(2, leaf1.GetParam(2));
                tree.SetParam(1, leaf1.GetParam(1));
                return true; // rerun cIf optimization
            }
        }

        CodeTree<Value_t>& branch1 = tree.GetParam(1);
        CodeTree<Value_t>& branch2 = tree.GetParam(2);

        if(branch1.IsIdenticalTo(branch2))
        {
            // TEST 20/ifnop
            // If both branches of an If() are identical, the test becomes unnecessary
            // NOTE: Possible side-effect on condition removed
            tree.Become(tree.GetParam(1));
            return true; // rerun optimization (opcode changed)
        }

        const OPCODE op1 = branch1.GetOpcode();
        const OPCODE op2 = branch2.GetOpcode();
        if(op1 == op2)
        {
            // TEST 20/if_extract_sin
            // TEST 20/if_extract_abs
            // If both branches apply the same unary function to different values,
            // extract the function. E.g. if(x,sin(a),sin(b)) -> sin(if(x,a,b))
            if(branch1.GetParamCount() == 1)
            {
                CodeTree<Value_t> changed_if;
                changed_if.SetOpcode(tree.GetOpcode());
                changed_if.AddParamMove(tree.GetParam(0));
                changed_if.AddParam(branch1.GetParam(0));
                changed_if.AddParam(branch2.GetParam(0));
                changed_if.Rehash();
                tree.SetOpcode(op1);
                tree.DelParams();
                tree.AddParamMove(changed_if);
                return true; // rerun optimization (opcode changed)
            }
            // TEST 20/if_extract_div
            // If both branches apply the same binary function to a set of parameters
            // where only one parameter differs, extract the function.
            // E.g. if(x, y/2, z/2) --> if(x, y,z)/2 (integer mode)
            if(branch1.GetParamCount() == 2
            && branch2.GetParamCount() == 2)
            {
                if(branch1.GetParam(0).IsIdenticalTo(branch2.GetParam(0)))
                {
                    CodeTree<Value_t> param0 = branch1.GetParam(0);
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.AddParamMove(tree.GetParam(0));
                    changed_if.AddParam(branch1.GetParam(1));
                    changed_if.AddParam(branch2.GetParam(1));
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.DelParams();
                    tree.AddParamMove(param0);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
                if(branch1.GetParam(1).IsIdenticalTo(branch2.GetParam(1)))
                {
                    CodeTree<Value_t> param1 = branch1.GetParam(1);
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.AddParamMove(tree.GetParam(0));
                    changed_if.AddParam(branch1.GetParam(0));
                    changed_if.AddParam(branch2.GetParam(0));
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.DelParams();
                    tree.AddParamMove(changed_if);
                    tree.AddParamMove(param1);
                    return true; // rerun optimization (opcode changed)
                }
            }
            // TEST 20/if_extract_add
            // TEST 20/if_extract_mul
            // TEST 20/if_extract_min
            if(op1 == cAdd    || op1 == cMul
            || op1 == cAnd    || op1 == cOr
            || op1 == cAbsAnd || op1 == cAbsOr
            || op1 == cMin    || op1 == cMax)
            {
                // If the two commutative groups contain one
                // or more identical values, extract them.
                std::vector<CodeTree<Value_t> > overlap;
                for(size_t a=branch1.GetParamCount(); a-- > 0; )
                {
                    for(size_t b=branch2.GetParamCount(); b-- > 0; )
                    {
                        if(branch1.GetParam(a).IsIdenticalTo(branch2.GetParam(b)))
                        {
                            if(overlap.empty()) { branch1.CopyOnWrite(); branch2.CopyOnWrite(); }
                            overlap.push_back(branch1.GetParam(a));
                            branch2.DelParam(b);
                            branch1.DelParam(a);
                            break;
                        }
                    }
                }
                if(!overlap.empty())
                {
                    branch1.Rehash();
                    branch2.Rehash();
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.SetParamsMove(overlap);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
            }
        }
        // if(x, y+z, y) -> if(x, z,0)+y
        if(op1 == cAdd
        || op1 == cMul
        || (op1 == cAnd && IsLogicalValue(branch2))
        || (op1 == cOr  && IsLogicalValue(branch2))
          )
        {
            // TEST 20/if_extract_add1
            // TEST 20/if_extract_mul1
            // TEST 20/if_extract_and1_l
            // TEST 20/if_extract_and1_nl
            // TEST 20/if_extract_or1_l
            // TEST 20/if_extract_or1_nl
            for(size_t a=branch1.GetParamCount(); a-- > 0; )
                if(branch1.GetParam(a).IsIdenticalTo(branch2))
                {
                    branch1.CopyOnWrite();
                    branch1.DelParam(a);
                    branch1.Rehash();
                    CodeTree<Value_t> branch2_backup = branch2;
                    branch2 = CodeTreeImmed( Value_t( (op1==cAdd||op1==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.AddParamMove(branch2_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, y&z, !!y) -> if(x, z, 1) & y
        if((op1 == cAnd || op1 == cOr) && op2 == cNotNot)
        {
            CodeTree<Value_t>& branch2op = branch2.GetParam(0);
            for(size_t a=branch1.GetParamCount(); a-- > 0; )
                if(branch1.GetParam(a).IsIdenticalTo(branch2op))
                {
                    branch1.CopyOnWrite();
                    branch1.DelParam(a);
                    branch1.Rehash();
                    CodeTree<Value_t> branch2_backup = branch2op;
                    branch2 = CodeTreeImmed( Value_t( (op1==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.AddParamMove(branch2_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, y, y+z) -> if(x, 0,z)+y
        if(op2 == cAdd
        || op2 == cMul
        || (op2 == cAnd && IsLogicalValue(branch1))
        || (op2 == cOr  && IsLogicalValue(branch1))
          )
        {
            // TEST 20/if_extract_add2
            // TEST 20/if_extract_mul2
            // TEST 20/if_extract_and2_l
            // TEST 20/if_extract_and2_nl
            // TEST 20/if_extract_or2_l
            // TEST 20/if_extract_or2_nl
            for(size_t a=branch2.GetParamCount(); a-- > 0; )
                if(branch2.GetParam(a).IsIdenticalTo(branch1))
                {
                    branch2.CopyOnWrite();
                    branch2.DelParam(a);
                    branch2.Rehash();
                    CodeTree<Value_t> branch1_backup = branch1;
                    branch1 = CodeTreeImmed( Value_t( (op2==cAdd||op2==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op2);
                    tree.AddParamMove(branch1_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, !!y, y&z) -> if(x, 1, z) & y
        if((op2 == cAnd || op2 == cOr) && op1 == cNotNot)
        {
            CodeTree<Value_t>& branch1op = branch1.GetParam(0);
            for(size_t a=branch2.GetParamCount(); a-- > 0; )
                if(branch2.GetParam(a).IsIdenticalTo(branch1op))
                {
                    branch2.CopyOnWrite();
                    branch2.DelParam(a);
                    branch2.Rehash();
                    CodeTree<Value_t> branch1_backup = branch1op;
                    branch1 = CodeTreeImmed( Value_t( (op2==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op2);
                    tree.AddParamMove(branch1_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        return false; // No changes
    }
}

#line 1 "fpoptimizer/logic_powoperations.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"

#include <limits>

/****
#ifdef _MSC_VER
#include <float.h>
#define isinf(x) (!_finite(x))
#endif
*/

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /**************************************/
    /* OPERATIONS DONE TO POW()           */
    /**************************************/

    template<typename Value_t>
    int maxFPExponent()
    {
        return std::numeric_limits<Value_t>::max_exponent;
    }

    template<typename Value_t>
    bool fPExponentIsTooLarge(Value_t base, Value_t exponent)
    {
        if(base < Value_t(0)) return true;
        if(fp_equal(base, Value_t(0)) || fp_equal(base, Value_t(1)))
            return false;
        return exponent >= Value_t(maxFPExponent<Value_t>()) / fp_log2(base);
    }

    template<typename Value_t>
    bool ConstantFolding_PowOperations(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cPow);

        if(tree.GetParam(0).IsImmed()
        && tree.GetParam(1).IsImmed())
        {
            Value_t const_value = fp_pow(tree.GetParam(0).GetImmed(),
                                        tree.GetParam(1).GetImmed());
            tree.ReplaceWithImmed(const_value);
            return false;
        }
        if(tree.GetParam(1).IsImmed()
        && fp_equal(tree.GetParam(1).GetImmed(), Value_t(1)))
        {
            // Used to be: float(getimmed()) == 1.0
            // Conversion through a float type value gets rid of
            // awkward abs(x)^1 generated from exp(log(x^6)/6),
            // without sacrificing as much precision as fp_equal() does.
            // x^1 = x
            tree.Become(tree.GetParam(0));
            return true; // rerun optimization (opcode changed)
        }
        if(tree.GetParam(0).IsImmed()
        && fp_equal(tree.GetParam(0).GetImmed(), Value_t(1)))
        {
            // 1^x = 1
            tree.ReplaceWithImmed(1);
            return false;
        }

        // 5^(20*x) = (5^20)^x
        if(tree.GetParam(0).IsImmed()
        && tree.GetParam(1).GetOpcode() == cMul)
        {
            bool changes = false;
            Value_t base_immed = tree.GetParam(0).GetImmed();
            CodeTree<Value_t> mulgroup = tree.GetParam(1);
            for(size_t a=mulgroup.GetParamCount(); a-->0; )
                if(mulgroup.GetParam(a).IsImmed())
                {
                    Value_t imm = mulgroup.GetParam(a).GetImmed();
                    //if(imm >= 0.0)
                    {
                        /****
                        Value_t new_base_immed = fp_pow(base_immed, imm);
                        if(isinf(new_base_immed)
                        || fp_equal(new_base_immed, Value_t(0)))
                        {
                            // It produced an infinity. Do not change.
                            break;
                        }
                        */
                        if(fPExponentIsTooLarge(base_immed, imm))
                            break;

                        Value_t new_base_immed = fp_pow(base_immed, imm);
                        if(fp_equal(new_base_immed, Value_t(0)) || fp_equal(new_base_immed, Value_t(1)))
                            break;

                        if(!changes)
                        {
                            changes = true;
                            mulgroup.CopyOnWrite();
                        }
                        base_immed = new_base_immed;
                        mulgroup.DelParam(a);
                        break; //
                    }
                }
            if(changes)
            {
                mulgroup.Rehash();
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before pow-mul change: "; DumpTree(tree);
                std::cout << "\n";
            #endif
                tree.GetParam(0).Become(CodeTreeImmed<Value_t> (base_immed));
                tree.GetParam(1).Become(mulgroup);
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After pow-mul change: "; DumpTree(tree);
                std::cout << "\n";
            #endif
            }
        }
        // (x*20)^2 = x^2 * 20^2
        if(tree.GetParam(1).IsImmed()
        && tree.GetParam(0).GetOpcode() == cMul)
        {
            Value_t exponent_immed = tree.GetParam(1).GetImmed();
            Value_t factor_immed   = 1.0;
            bool changes = false;
            CodeTree<Value_t>& mulgroup = tree.GetParam(0);
            for(size_t a=mulgroup.GetParamCount(); a-->0; )
                if(mulgroup.GetParam(a).IsImmed())
                {
                    Value_t imm = mulgroup.GetParam(a).GetImmed();
                    //if(imm >= 0.0)
                    {
                        /****
                        Value_t new_factor_immed = fp_pow(imm, exponent_immed);
                        if(isinf(new_factor_immed)
                        || fp_equal(new_factor_immed, Value_t(0)))
                        {
                            // It produced an infinity. Do not change.
                            break;
                        }
                        */
                        if(fPExponentIsTooLarge(imm, exponent_immed))
                            break;

                        Value_t new_factor_immed = fp_pow(imm, exponent_immed);
                        if(fp_equal(new_factor_immed, Value_t(0)))
                            break;

                        if(!changes)
                        {
                            changes = true;
                            mulgroup.CopyOnWrite();
                        }
                        factor_immed *= new_factor_immed;
                        mulgroup.DelParam(a);
                        break; //
                    }
                }
            if(changes)
            {
                mulgroup.Rehash();
                CodeTree<Value_t> newpow;
                newpow.SetOpcode(cPow);
                newpow.SetParamsMove(tree.GetParams());
                newpow.Rehash(false);
                tree.SetOpcode(cMul);
                tree.AddParamMove(newpow);
                tree.AddParam( CodeTreeImmed<Value_t>(factor_immed) );
                return true; // rerun optimization (opcode changed)
            }
        }

        // (x^3)^2 = x^6
        // NOTE: If 3 is even and 3*2 is not, x must be changed to abs(x).
        if(tree.GetParam(0).GetOpcode() == cPow
        && tree.GetParam(1).IsImmed()
        && tree.GetParam(0).GetParam(1).IsImmed())
        {
            Value_t a = tree.GetParam(0).GetParam(1).GetImmed();
            Value_t b = tree.GetParam(1).GetImmed();
            Value_t c = a * b; // new exponent
            if(isEvenInteger(a) // a is an even int?
            && !isEvenInteger(c)) // c is not?
            {
                CodeTree<Value_t> newbase;
                newbase.SetOpcode(cAbs);
                newbase.AddParam(tree.GetParam(0).GetParam(0));
                newbase.Rehash();
                tree.SetParamMove(0, newbase);
            }
            else
                tree.SetParam(0, tree.GetParam(0).GetParam(0));
            tree.SetParam(1, CodeTreeImmed<Value_t>(c));
        }
        return false; // No changes that require a rerun
    }
}

#line 1 "fpoptimizer/logic_comparisons.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /*****************************************/
    /* RANGE-BASED TREATMENTS TO COMPARISONS */
    /*****************************************/

    struct RangeComparisonData
    {
        enum Decision
        {
            MakeFalse=0,
            MakeTrue=1,
            MakeNEqual=2,
            MakeEqual=3,
            MakeNotNotP0=4,
            MakeNotNotP1=5,
            MakeNotP0=6,
            MakeNotP1=7,
            Unchanged=8
        };
        enum WhatDoWhenCase
        {
            Never =0,
            Eq0   =1, // val==0
            Eq1   =2, // val==1
            Gt0Le1=3, // val>0 && val<=1
            Ge0Lt1=4  // val>=0 && val<1
        };

        Decision if_identical; // What to do when operands are identical
        Decision if_always[4]; // What to do if Always <, <=, >, >=
        struct { Decision what : 4; WhatDoWhenCase when : 4; }
            p0_logical_a, p1_logical_a,
            p0_logical_b, p1_logical_b;

        template<typename Value_t>
        Decision Analyze(const CodeTree<Value_t>& a, const CodeTree<Value_t>& b) const
        {
            if(a.IsIdenticalTo(b))
                return if_identical;

            range<Value_t> p0 = CalculateResultBoundaries(a);
            range<Value_t> p1 = CalculateResultBoundaries(b);
            if(p0.max.known && p1.min.known)
            {
                if(p0.max.val <  p1.min.val && if_always[0] != Unchanged)
                    return if_always[0]; // p0 < p1
                if(p0.max.val <= p1.min.val && if_always[1] != Unchanged)
                    return if_always[1]; // p0 <= p1
            }
            if(p0.min.known && p1.max.known)
            {
                if(p0.min.val >  p1.max.val && if_always[2] != Unchanged)
                    return if_always[2]; // p0 > p1
                if(p0.min.val >= p1.max.val && if_always[3] != Unchanged)
                    return if_always[3]; // p0 >= p1
            }

            if(IsLogicalValue(a))
            {
                if(p0_logical_a.what != Unchanged)
                    if(TestCase(p0_logical_a.when, p1)) return p0_logical_a.what;
                if(p0_logical_b.what != Unchanged)
                    if(TestCase(p0_logical_b.when, p1)) return p0_logical_b.what;
            }
            if(IsLogicalValue(b))
            {
                if(p1_logical_a.what != Unchanged)
                    if(TestCase(p1_logical_a.when, p0)) return p1_logical_a.what;
                if(p1_logical_b.what != Unchanged)
                    if(TestCase(p1_logical_b.when, p0)) return p1_logical_b.what;
            }
            return Unchanged;
        }

        template<typename Value_t>
        static bool TestCase(WhatDoWhenCase when, const range<Value_t>& p)
        {
            if(!p.min.known || !p.max.known) return false;
            switch(when)
            {
                case Eq0: return p.min.val==Value_t(0.0) && p.max.val==p.min.val;
                case Eq1: return p.min.val==Value_t(1.0) && p.max.val==p.max.val;
                case Gt0Le1: return p.min.val>Value_t(0) && p.max.val<=Value_t(1);
                case Ge0Lt1: return p.min.val>=Value_t(0) && p.max.val<Value_t(1);
                default:;
            }
            return false;
        }
    };

    namespace RangeComparisonsData
    {
        static const RangeComparisonData Data[6] =
        {
            // cEqual:
            // Case:      p0 == p1  Antonym: p0 != p1
            // Synonym:   p1 == p0  Antonym: p1 != p0
            { RangeComparisonData::MakeTrue,  // If identical: always true
              {RangeComparisonData::MakeFalse,  // If Always p0 < p1: always false
               RangeComparisonData::Unchanged,
               RangeComparisonData::MakeFalse,  // If Always p0 > p1: always false
               RangeComparisonData::Unchanged},
             // NotNot(p0) if p1==1    NotNot(p1) if p0==1
             //    Not(p0) if p1==0       Not(p1) if p0==0
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotP0, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotP1, RangeComparisonData::Eq0}
            },
            // cNEqual:
            // Case:      p0 != p1  Antonym: p0 == p1
            // Synonym:   p1 != p0  Antonym: p1 == p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeTrue,  // If Always p0 < p1: always true
               RangeComparisonData::Unchanged,
               RangeComparisonData::MakeTrue,  // If Always p0 > p1: always true
               RangeComparisonData::Unchanged},
             // NotNot(p0) if p1==0    NotNot(p1) if p0==0
             //    Not(p0) if p1==1       Not(p1) if p0==1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotP0, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotP1, RangeComparisonData::Eq1}
            },
            // cLess:
            // Case:      p0 < p1   Antonym: p0 >= p1
            // Synonym:   p1 > p0   Antonym: p1 <= p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeTrue,  // If Always p0  < p1: always true
               RangeComparisonData::MakeNEqual,
               RangeComparisonData::MakeFalse, // If Always p0 > p1: always false
               RangeComparisonData::MakeFalse},// If Always p0 >= p1: always false
             // Not(p0)   if p1>0 & p1<=1    --   NotNot(p1) if p0>=0 & p0<1
              {RangeComparisonData::MakeNotP0,    RangeComparisonData::Gt0Le1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cLessOrEq:
            // Case:      p0 <= p1  Antonym: p0 > p1
            // Synonym:   p1 >= p0  Antonym: p1 < p0
            { RangeComparisonData::MakeTrue,   // If identical: always true
              {RangeComparisonData::Unchanged, // If Always p0  < p1: ?
               RangeComparisonData::MakeTrue,  // If Always p0 <= p1: always true
               RangeComparisonData::MakeFalse, // If Always p0  > p1: always false
               RangeComparisonData::MakeEqual},// If Never  p0  < p1:  use cEqual
             // Not(p0)    if p1>=0 & p1<1   --   NotNot(p1) if p0>0 & p0<=1
              {RangeComparisonData::MakeNotP0,    RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Gt0Le1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cGreater:
            // Case:      p0 >  p1  Antonym: p0 <= p1
            // Synonym:   p1 <  p0  Antonym: p1 >= p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeFalse, // If Always p0  < p1: always false
               RangeComparisonData::MakeFalse, // If Always p0 <= p1: always false
               RangeComparisonData::MakeTrue,  // If Always p0  > p1: always true
               RangeComparisonData::MakeNEqual},
             // NotNot(p0) if p1>=0 & p1<1   --   Not(p1)   if p0>0 & p0<=1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::MakeNotP1,    RangeComparisonData::Gt0Le1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cGreaterOrEq:
            // Case:      p0 >= p1  Antonym: p0 < p1
            // Synonym:   p1 <= p0  Antonym: p1 > p0
            { RangeComparisonData::MakeTrue,   // If identical: always true
              {RangeComparisonData::MakeFalse, // If Always p0  < p1: always false
               RangeComparisonData::MakeEqual, // If Always p0 >= p1: always true
               RangeComparisonData::Unchanged, // If always p0  > p1: ?
               RangeComparisonData::MakeTrue}, // If Never  p0  > p1:  use cEqual
             // NotNot(p0) if p1>0 & p1<=1   --   Not(p1)    if p0>=0 & p0<1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Gt0Le1},
              {RangeComparisonData::MakeNotP1,    RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            }
        };
    }

    template<typename Value_t>
    bool ConstantFolding_Comparison(CodeTree<Value_t>& tree)
    {
        using namespace RangeComparisonsData;

        assert(tree.GetOpcode() >= cEqual && tree.GetOpcode() <= cGreaterOrEq);

        switch(Data[tree.GetOpcode()-cEqual].
            Analyze(tree.GetParam(0), tree.GetParam(1)))
        {
            case RangeComparisonData::MakeFalse:
                tree.ReplaceWithImmed(0); return true;
            case RangeComparisonData::MakeTrue:
                tree.ReplaceWithImmed(1); return true;
            case RangeComparisonData::MakeEqual:  tree.SetOpcode(cEqual); return true;
            case RangeComparisonData::MakeNEqual: tree.SetOpcode(cNEqual); return true;
            case RangeComparisonData::MakeNotNotP0: tree.SetOpcode(cNotNot); tree.DelParam(1); return true;
            case RangeComparisonData::MakeNotNotP1: tree.SetOpcode(cNotNot); tree.DelParam(0); return true;
            case RangeComparisonData::MakeNotP0: tree.SetOpcode(cNot); tree.DelParam(1); return true;
            case RangeComparisonData::MakeNotP1: tree.SetOpcode(cNot); tree.DelParam(0); return true;
            case RangeComparisonData::Unchanged:;
        }

        // Any reversible functions:
        //   sin(x)  -> ASIN: Not doable, x can be cyclic
        //   asin(x) -> SIN: doable.
        //                   Invalid combinations are caught by
        //                   range-estimation. Threshold is at |pi/2|.
        //   acos(x) -> COS: doable.
        //                   Invalid combinations are caught by
        //                   range-estimation. Note that though
        //                   the range is contiguous, it is direction-flipped.
        //    log(x) -> EXP: no problem
        //   exp2, exp10: Converted to cPow, done by grammar.
        //   atan(x) -> TAN: doable.
        //                   Invalid combinations are caught by
        //                   range-estimation. Threshold is at |pi/2|.
        //   sinh(x) -> ASINH: no problem
        //   tanh(x) -> ATANH: no problem, but atanh is limited to -1..1
        //                     Invalid combinations are caught by
        //                     range-estimation, but the exact value
        //                     of 1.0 still needs checking, because
        //                     it involves infinity.
        if(tree.GetParam(1).IsImmed())
            switch(tree.GetParam(0).GetOpcode())
            {
                case cAsin:
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_sin(tree.GetParam(1).GetImmed())));
                    return true;
                case cAcos:
                    // -1..+1 --> pi..0 (polarity-flipping)
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_cos(tree.GetParam(1).GetImmed())));
                    tree.SetOpcode( tree.GetOpcode()==cLess ? cGreater
                                  : tree.GetOpcode()==cLessOrEq ? cGreaterOrEq
                                  : tree.GetOpcode()==cGreater ? cLess
                                  : tree.GetOpcode()==cGreaterOrEq ? cLessOrEq
                                  : tree.GetOpcode() );
                    return true;
                case cAtan:
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_tan(tree.GetParam(1).GetImmed())));
                    return true;
                case cLog:
                    // Different logarithms have a constant-multiplication,
                    // which is no problem.
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_exp(tree.GetParam(1).GetImmed())));
                    return true;
                case cSinh:
                    tree.SetParam(0, tree.GetParam(0).GetParam(0));
                    tree.SetParam(1, CodeTreeImmed(fp_asinh(tree.GetParam(1).GetImmed())));
                    return true;
                case cTanh:
                    if(fp_less(fp_abs(tree.GetParam(1).GetImmed()), Value_t(1)))
                    {
                        tree.SetParam(0, tree.GetParam(0).GetParam(0));
                        tree.SetParam(1, CodeTreeImmed(fp_atanh(tree.GetParam(1).GetImmed())));
                        return true;
                    }
                    break;
                default: break;
            }

        return false;
    }
}

#line 1 "fpoptimizer/codetree.cc"
#include <list>
#include <algorithm>

// line removed for fpoptimizer.cc: #include "rangeestimation.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh" // for DEBUG_SUBSTITUTIONS
// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "consts.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
#ifdef DEBUG_SUBSTITUTIONS
    void OutFloatHex(std::ostream& o, double d)
    {
        union { double d; uint_least64_t h; } data;
        data.d = d;
        o << "(" << std::hex << data.h << std::dec << ")";
    }
  #ifdef FP_SUPPORT_FLOAT_TYPE
    void OutFloatHex(std::ostream& o, float f)
    {
        union { float f; uint_least32_t h; } data;
        data.f = f;
        o << "(" << std::hex << data.h << std::dec << ")";
    }
  #endif
  #ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
    void OutFloatHex(std::ostream& o, long double ld)
    {
        union { long double ld;
                struct { uint_least64_t a; unsigned short b; } s; } data;
        data.ld = ld;
        o << "(" << std::hex << data.s.b << data.s.a << std::dec << ")";
    }
  #endif
  #ifdef FP_SUPPORT_LONG_INT_TYPE
    void OutFloatHex(std::ostream& o, long ld)
    {
        o << "(" << std::hex << ld << std::dec << ")";
    }
  #endif

#endif
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    CodeTree<Value_t>::CodeTree()
        : data(new CodeTreeData<Value_t> ()) // sets opcode to cNop
    {
    }

    template<typename Value_t>
    CodeTree<Value_t>::CodeTree(const Value_t& i, typename CodeTree<Value_t>::ImmedTag)
        : data(new CodeTreeData<Value_t>(i))
    {
        data->Recalculate_Hash_NoRecursion();
    }

#ifdef FP_SUPPORT_CXX11_MOVE
    template<typename Value_t>
    CodeTree<Value_t>::CodeTree(Value_t&& i, typename CodeTree<Value_t>::ImmedTag)
        : data(new CodeTreeData<Value_t>(std::move(i)))
    {
        data->Recalculate_Hash_NoRecursion();
    }
#endif

    template<typename Value_t>
    CodeTree<Value_t>::CodeTree(unsigned v, typename CodeTree<Value_t>::VarTag)
        : data(new CodeTreeData<Value_t> (VarBegin, v))
    {
        data->Recalculate_Hash_NoRecursion();
    }

    template<typename Value_t>
    CodeTree<Value_t>::CodeTree(FUNCTIONPARSERTYPES::OPCODE o, typename CodeTree<Value_t>::OpcodeTag)
        : data(new CodeTreeData<Value_t> (o))
    {
        data->Recalculate_Hash_NoRecursion();
    }

    template<typename Value_t>
    CodeTree<Value_t>::CodeTree(FUNCTIONPARSERTYPES::OPCODE o, unsigned f, typename CodeTree<Value_t>::FuncOpcodeTag)
        : data(new CodeTreeData<Value_t> (o, f))
    {
        data->Recalculate_Hash_NoRecursion();
    }

    template<typename Value_t>
    CodeTree<Value_t>::CodeTree(const CodeTree<Value_t>& b, typename CodeTree<Value_t>::CloneTag)
        : data(new CodeTreeData<Value_t>(*b.data))
    {
    }

    template<typename Value_t>
    CodeTree<Value_t>::~CodeTree()
    {
    }

    template<typename Value_t>
    void CodeTree<Value_t>::ReplaceWithImmed(const Value_t& i)
    {
      #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Replacing "; DumpTree(*this);
        if(IsImmed())
            OutFloatHex(std::cout, GetImmed());
        std::cout << " with const value " << i;
        OutFloatHex(std::cout, i);
        std::cout << "\n";
      #endif
        data = new CodeTreeData<Value_t> (i);
    }

    template<typename Value_t>
    struct ParamComparer
    {
        bool operator() (const CodeTree<Value_t>& a, const CodeTree<Value_t>& b) const
        {
            if(a.GetDepth() != b.GetDepth())
                return a.GetDepth() < b.GetDepth();
            return a.GetHash() < b.GetHash();
        }
    };

    template<typename Value_t>
    void CodeTreeData<Value_t>::Sort()
    {
        /* If the tree is commutative, order the parameters
         * in a set order in order to make equality tests
         * efficient in the optimizer
         */
        switch(Opcode)
        {
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd: case cAbsAnd:
            case cOr: case cAbsOr:
            case cHypot:
            case cEqual:
            case cNEqual:
                std::sort(Params.begin(), Params.end(), ParamComparer<Value_t>());
                break;
            case cLess:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cGreater; }
                break;
            case cLessOrEq:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cGreaterOrEq; }
                break;
            case cGreater:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cLess; }
                break;
            case cGreaterOrEq:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cLessOrEq; }
                break;
            /*
            case cDiv:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cRDiv; }
                break;
            case cRDiv:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cDiv; }
                break;
            case cSub:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cRSub; }
                break;
            case cRSub:
                if(ParamComparer<Value_t>() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cSub; }
                break;
            */
            default:
                break;
        }
    }

    template<typename Value_t>
    void CodeTree<Value_t>::AddParam(const CodeTree<Value_t>& param)
    {
        //std::cout << "AddParam called\n";
        data->Params.push_back(param);
    }
    template<typename Value_t>
    void CodeTree<Value_t>::AddParamMove(CodeTree<Value_t>& param)
    {
        data->Params.push_back(CodeTree<Value_t>());
        data->Params.back().swap(param);
    }
    template<typename Value_t>
    void CodeTree<Value_t>::SetParam(size_t which, const CodeTree<Value_t>& b)
    {
        DataP slot_holder ( data->Params[which].data );
        data->Params[which] = b;
    }
    template<typename Value_t>
    void CodeTree<Value_t>::SetParamMove(size_t which, CodeTree<Value_t>& b)
    {
        DataP slot_holder ( data->Params[which].data );
        data->Params[which].swap(b);
    }

    template<typename Value_t>
    void CodeTree<Value_t>::AddParams(const std::vector<CodeTree<Value_t> >& RefParams)
    {
        data->Params.insert(data->Params.end(), RefParams.begin(), RefParams.end());
    }
    template<typename Value_t>
    void CodeTree<Value_t>::AddParamsMove(std::vector<CodeTree<Value_t> >& RefParams)
    {
        size_t endpos = data->Params.size(), added = RefParams.size();
        data->Params.resize(endpos + added, CodeTree<Value_t>());
        for(size_t p=0; p<added; ++p)
            data->Params[endpos+p].swap( RefParams[p] );
    }
    template<typename Value_t>
    void CodeTree<Value_t>::AddParamsMove(std::vector<CodeTree<Value_t> >& RefParams, size_t replacing_slot)
    {
        DataP slot_holder ( data->Params[replacing_slot].data );
        DelParam(replacing_slot);
        AddParamsMove(RefParams);
    /*
        const size_t n_added = RefParams.size();
        const size_t oldsize = data->Params.size();
        const size_t newsize = oldsize + n_added - 1;
        if(RefParams.empty())
            DelParam(replacing_slot);
        else
        {
            //    0 1 2 3 4 5 6 7 8 9 10 11
            //    a a a a X b b b b b
            //    a a a a Y Y Y b b b b  b
            //
            //   replacing_slot = 4
            //   n_added = 3
            //   oldsize = 10
            //   newsize = 12
            //   tail_length = 5

            data->Params.resize(newsize);
            data->Params[replacing_slot].data = 0;
            const size_t tail_length = oldsize - replacing_slot -1;
            for(size_t tail=0; tail<tail_length; ++tail)
                data->Params[newsize-1-tail].data.UnsafeSetP(
                &*data->Params[newsize-1-tail-(n_added-1)].data);
            for(size_t head=1; head<n_added; ++head)
                data->Params[replacing_slot+head].data.UnsafeSetP( 0 );
            for(size_t p=0; p<n_added; ++p)
                data->Params[replacing_slot+p].swap( RefParams[p] );
        }
    */
    }

    template<typename Value_t>
    void CodeTree<Value_t>::SetParams(const std::vector<CodeTree<Value_t> >& RefParams)
    {
        //std::cout << "SetParams called" << (do_clone ? ", clone" : ", no clone") << "\n";
        std::vector<CodeTree<Value_t> > tmp(RefParams);
        data->Params.swap(tmp);
    }

    template<typename Value_t>
    void CodeTree<Value_t>::SetParamsMove(std::vector<CodeTree<Value_t> >& RefParams)
    {
        data->Params.swap(RefParams);
        RefParams.clear();
    }

#ifdef FP_SUPPORT_CXX11_MOVE
    template<typename Value_t>
    void CodeTree<Value_t>::SetParams(std::vector<CodeTree<Value_t> >&& RefParams)
    {
        //std::cout << "SetParams&& called\n";
        SetParamsMove(RefParams);
    }
#endif

    template<typename Value_t>
    void CodeTree<Value_t>::DelParam(size_t index)
    {
        std::vector<CodeTree<Value_t> >& Params = data->Params;
        //std::cout << "DelParam(" << index << ") called\n";
    #ifdef FP_SUPPORT_CXX11_MOVE
        /* rvalue reference semantics makes this optimal */
        Params.erase( Params.begin() + index );
    #else
        /* This labor evades the need for refcount +1/-1 shuffling */
        Params[index].data = 0;
        for(size_t p=index; p+1<Params.size(); ++p)
            Params[p].data.UnsafeSetP( &*Params[p+1].data );
        Params[Params.size()-1].data.UnsafeSetP( 0 );
        Params.resize(Params.size()-1);
    #endif
    }

    template<typename Value_t>
    void CodeTree<Value_t>::DelParams()
    {
        data->Params.clear();
    }

    template<typename Value_t>
    bool CodeTree<Value_t>::IsIdenticalTo(const CodeTree<Value_t>& b) const
    {
        //if((!&*data) != (!&*b.data)) return false;
        if(&*data == &*b.data) return true;
        return data->IsIdenticalTo(*b.data);
    }

    template<typename Value_t>
    bool CodeTreeData<Value_t>::IsIdenticalTo(const CodeTreeData<Value_t>& b) const
    {
        if(Hash   != b.Hash) return false; // a quick catch-all
        if(Opcode != b.Opcode) return false;
        switch(Opcode)
        {
            case cImmed:   return fp_equal(Value, b.Value);
            case VarBegin: return Var_or_Funcno == b.Var_or_Funcno;
            case cFCall:
            case cPCall:   if(Var_or_Funcno != b.Var_or_Funcno) return false; break;
            default: break;
        }
        if(Params.size() != b.Params.size()) return false;
        for(size_t a=0; a<Params.size(); ++a)
        {
            if(!Params[a].IsIdenticalTo(b.Params[a])) return false;
        }
        return true;
    }

    template<typename Value_t>
    void CodeTree<Value_t>::Become(const CodeTree<Value_t>& b)
    {
        if(&b != this && &*data != &*b.data)
        {
            DataP tmp = b.data;
            CopyOnWrite();
            data.swap(tmp);
        }
    }

    template<typename Value_t>
    void CodeTree<Value_t>::CopyOnWrite()
    {
        if(GetRefCount() > 1)
            data = new CodeTreeData<Value_t>(*data);
    }

    template<typename Value_t>
    CodeTree<Value_t> CodeTree<Value_t>::GetUniqueRef()
    {
        if(GetRefCount() > 1)
            return CodeTree<Value_t>(*this, CloneTag());
        return *this;
    }

    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData()
        : RefCount(0),
          Opcode(cNop),
          Value(), Var_or_Funcno(),
          Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
    }

    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData(const CodeTreeData& b)
        : RefCount(0),
          Opcode(b.Opcode),
          Value(b.Value),
          Var_or_Funcno(b.Var_or_Funcno),
          Params(b.Params),
          Hash(b.Hash),
          Depth(b.Depth),
          OptimizedUsing(b.OptimizedUsing)
    {
    }

    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData(const Value_t& i)
        : RefCount(0),
          Opcode(cImmed),
          Value(i), Var_or_Funcno(),
          Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
    }

#ifdef FP_SUPPORT_CXX11_MOVE
    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData(CodeTreeData<Value_t>&& b)
        : RefCount(0),
          Opcode(b.Opcode),
          Value(std::move(b.Value)),
          Var_or_Funcno(b.Var_or_Funcno),
          Params(std::move(b.Params)),
          Hash(b.Hash),
          Depth(b.Depth),
          OptimizedUsing(b.OptimizedUsing)
    {
    }

    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData(Value_t&& i)
        : RefCount(0),
          Opcode(cImmed),
          Value(std::move(i)), Var_or_Funcno(),
          Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
    }
#endif

    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData(FUNCTIONPARSERTYPES::OPCODE o)
        : RefCount(0),
          Opcode(o),
          Value(), Var_or_Funcno(),
          Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
    }

    template<typename Value_t>
    CodeTreeData<Value_t>::CodeTreeData(FUNCTIONPARSERTYPES::OPCODE o, unsigned f)
        : RefCount(0),
          Opcode(o),
          Value(), Var_or_Funcno(f),
          Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
    }
}

/*
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template class CodeTree<type>; \
    template struct CodeTreeData<type>;
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/debug.cc"
// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "opcodename.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <sstream>
#include <string>
#include <map>
#include <set>
#include <iostream>

using namespace FUNCTIONPARSERTYPES;

#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
namespace
{
    template<typename Value_t>
    void DumpHashesFrom(
        const FPoptimizer_CodeTree::CodeTree<Value_t>& tree,
        std::map<fphash_t, std::set<std::string> >& done,
        std::ostream& o)
    {
        for(size_t a=0; a<tree.GetParamCount(); ++a)
            DumpHashesFrom(tree.GetParam(a), done, o);

        std::ostringstream buf;
        DumpTree(tree, buf);
        done[tree.GetHash()].insert(buf.str());
    }
}
#endif

namespace FPoptimizer_CodeTree
{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
    template<typename Value_t>
    void DumpHashes(const CodeTree<Value_t>& tree, std::ostream& o)
    {
        std::map<fphash_t, std::set<std::string> > done;
        DumpHashesFrom(tree, done, o);

        for(std::map<fphash_t, std::set<std::string> >::const_iterator
            i = done.begin();
            i != done.end();
            ++i)
        {
            const std::set<std::string>& flist = i->second;
            if(flist.size() != 1) o << "ERROR - HASH COLLISION?\n";
            for(std::set<std::string>::const_iterator
                j = flist.begin();
                j != flist.end();
                ++j)
            {
                o << '[' << std::hex << i->first.hash1
                              << ',' << i->first.hash2
                              << ']' << std::dec;
                o << ": " << *j << "\n";
            }
        }
    }

    template<typename Value_t>
    void DumpTree(const CodeTree<Value_t>& tree, std::ostream& o)
    {
        //o << "/*" << tree.Depth << "*/";
        const char* sep2 = " ";
        /*
        o << '[' << std::hex << tree.Hash.hash1
                      << ',' << tree.Hash.hash2
                      << ']' << std::dec;
        */
        switch(tree.GetOpcode())
        {
            case cImmed:
                o << tree.GetImmed();
                /*
                o << "(" << std::hex
                  << *(const uint_least64_t*)&tree.GetImmed()
                  << std::dec << ")";
                */
                return;
            case VarBegin: o << "Var" << (tree.GetVar() - VarBegin); return;
            case cAdd: sep2 = " +"; break;
            case cMul: sep2 = " *"; break;
            case cAnd: sep2 = " &"; break;
            case cOr: sep2 = " |"; break;
            case cPow: sep2 = " ^"; break;
            default:
                sep2 = " ";
                o << FP_GetOpcodeName(tree.GetOpcode());
                if(tree.GetOpcode() == cFCall || tree.GetOpcode() == cPCall)
                    o << ':' << tree.GetFuncNo();
        }
        o << '(';
        if(tree.GetParamCount() <= 1 && sep2[1]) o << (sep2+1) << ' ';
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            if(a > 0) o << ' ';

            DumpTree(tree.GetParam(a), o);

            if(a+1 < tree.GetParamCount()) o << sep2;
        }
        o << ')';
    }

    template<typename Value_t>
    void DumpTreeWithIndent(const CodeTree<Value_t>& tree, std::ostream& o, const std::string& indent)
    {
        o << '[' << std::hex << (void*)(&tree.GetParams())
                 << std::dec
                 << ',' << tree.GetRefCount()
                 << ']';
        o << indent << '_';

        switch(tree.GetOpcode())
        {
            case cImmed:   o << "cImmed " << tree.GetImmed(); o << '\n'; return;
            case VarBegin: o << "VarBegin " << (tree.GetVar() - VarBegin); o << '\n'; return;
            default:
                o << FP_GetOpcodeName(tree.GetOpcode());
                if(tree.GetOpcode() == cFCall || tree.GetOpcode() == cPCall)
                    o << ':' << tree.GetFuncNo();
                o << '\n';
        }
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            std::string ind = indent;
            for(size_t p=0; p < ind.size(); p+=2)
                if(ind[p] == '\\')
                    ind[p] = ' ';
            ind += (a+1 < tree.GetParamCount()) ? " |" : " \\";
            DumpTreeWithIndent(tree.GetParam(a), o, ind);
        }
        o << std::flush;
    }
#endif
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
#define FP_INSTANTIATE(type) \
    template void DumpHashes(const CodeTree<type>& tree, std::ostream& o); \
    template void DumpTree(const CodeTree<type>& tree, std::ostream& o); \
    template void DumpTreeWithIndent(const CodeTree<type>& tree, std::ostream& o, const std::string& indent);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
#endif
}
 */

#endif

#line 1 "fpoptimizer/grammar.cc"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"
// line removed for fpoptimizer.cc: #include "grammar.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"

// line removed for fpoptimizer.cc: #include "opcodename.hh"
using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

#include <cctype>

namespace FPoptimizer_Grammar
{
    unsigned ParamSpec_GetDepCode(const ParamSpec& b)
    {
        switch(b.first)
        {
            case ParamHolder:
            {
                const ParamSpec_ParamHolder* s = (const ParamSpec_ParamHolder*) b.second;
                return s->depcode;
            }
            case SubFunction:
            {
                const ParamSpec_SubFunction* s = (const ParamSpec_SubFunction*) b.second;
                return s->depcode;
            }
            default: break;
        }
        return 0;
    }

    template<typename Value_t>
    void DumpParam(const ParamSpec& parampair, std::ostream& o)
    {
        static const char ParamHolderNames[][2]  = {"%","&", "x","y","z","a","b","c"};

        //o << "/*p" << (&p-pack.plist) << "*/";
        unsigned constraints = 0;
        switch(parampair.first)
        {
            case NumConstant:
              { const ParamSpec_NumConstant<Value_t>& param = *(const ParamSpec_NumConstant<Value_t>*) parampair.second;
                using namespace FUNCTIONPARSERTYPES;
                o.precision(12);
                o << param.constvalue; break; }
            case ParamHolder:
              { const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                o << ParamHolderNames[param.index];
                constraints = param.constraints;
                break; }
            case SubFunction:
              { const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                constraints = param.constraints;
                if(param.data.match_type == GroupFunction)
                {
                    if(param.data.subfunc_opcode == cNeg)
                        { o << "-"; DumpParams<Value_t>(param.data.param_list, param.data.param_count, o); }
                    else if(param.data.subfunc_opcode == cInv)
                        { o << "/"; DumpParams<Value_t>(param.data.param_list, param.data.param_count, o); }
                    else
                    {
                        std::string opcode = FP_GetOpcodeName
                            ( (FUNCTIONPARSERTYPES::OPCODE) param.data.subfunc_opcode)
                            .substr(1);
                        for(size_t a=0; a<opcode.size(); ++a) opcode[a] = (char) std::toupper(opcode[a]);
                        o << opcode << "( ";
                        DumpParams<Value_t>(param.data.param_list, param.data.param_count, o);
                        o << " )";
                    }
                }
                else
                {
                    o << '(' << FP_GetOpcodeName
                        ( (FUNCTIONPARSERTYPES::OPCODE) param.data.subfunc_opcode)
                       << ' ';
                    if(param.data.match_type == PositionalParams) o << '[';
                    if(param.data.match_type == SelectedParams) o << '{';
                    DumpParams<Value_t>(param.data.param_list, param.data.param_count, o);
                    if(param.data.restholder_index != 0)
                        o << " <" << param.data.restholder_index << '>';
                    if(param.data.match_type == PositionalParams) o << "]";
                    if(param.data.match_type == SelectedParams) o << "}";
                    o << ')';
                }
                break; }
        }
        switch( ImmedConstraint_Value(constraints & ValueMask) )
        {
            case ValueMask: break;
            case Value_AnyNum: break;
            case Value_EvenInt:   o << "@E"; break;
            case Value_OddInt:    o << "@O"; break;
            case Value_IsInteger: o << "@I"; break;
            case Value_NonInteger:o << "@F"; break;
            case Value_Logical:   o << "@L"; break;
        }
        switch( ImmedConstraint_Sign(constraints & SignMask) )
        {
            case SignMask: break;
            case Sign_AnySign: break;
            case Sign_Positive:   o << "@P"; break;
            case Sign_Negative:   o << "@N"; break;
        }
        switch( ImmedConstraint_Oneness(constraints & OnenessMask) )
        {
            case OnenessMask: break;
            case Oneness_Any: break;
            case Oneness_One:     o << "@1"; break;
            case Oneness_NotOne:  o << "@M"; break;
        }
        switch( ImmedConstraint_Constness(constraints & ConstnessMask) )
        {
            case ConstnessMask: break;
            case Constness_Const:
                if(parampair.first == ParamHolder)
                {
                    const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                    if(param.index < 2) break;
                }
                o << "@C";
                break;
            case Constness_NotConst: o << "@V"; break;
            case Oneness_Any: break;
        }
    }

    template<typename Value_t>
    void DumpParams(unsigned paramlist, unsigned count, std::ostream& o)
    {
        for(unsigned a=0; a<count; ++a)
        {
            if(a > 0) o << ' ';
            const ParamSpec& param = ParamSpec_Extract<Value_t> (paramlist,a);
            DumpParam<Value_t> (param, o);
            unsigned depcode = ParamSpec_GetDepCode(param);
            if(depcode != 0)
                o << "@D" << depcode;
        }
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
#include <complex>
namespace FPoptimizer_Grammar
{
#define FP_INSTANTIATE(type) \
    template void DumpParams<type>(unsigned, unsigned, std::ostream& ); \
    template void DumpParam<type>(const ParamSpec&, std::ostream& );
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#line 1 "fpoptimizer/grammar_data.cc"
/* This file is automatically generated. Do not edit... */
// line removed for fpoptimizer.cc: #include "../fpoptimizer/consts.hh"
// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"
#include <algorithm>

/* 
#define grammar_optimize_abslogical grammar_optimize_abslogical_tweak
#define grammar_optimize_ignore_if_sideeffects grammar_optimize_ignore_if_sideeffects_tweak
#define grammar_optimize_nonshortcut_logical_evaluation grammar_optimize_nonshortcut_logical_evaluation_tweak
#define grammar_optimize_recreate grammar_optimize_recreate_tweak
#define grammar_optimize_round1 grammar_optimize_round1_tweak
#define grammar_optimize_round2 grammar_optimize_round2_tweak
#define grammar_optimize_round3 grammar_optimize_round3_tweak
#define grammar_optimize_round4 grammar_optimize_round4_tweak
#define grammar_optimize_shortcut_logical_evaluation grammar_optimize_shortcut_logical_evaluation_tweak
// line removed for fpoptimizer.cc: #include "../fpoptimizer/grammar.hh"
#undef grammar_optimize_abslogical
#undef grammar_optimize_ignore_if_sideeffects
#undef grammar_optimize_nonshortcut_logical_evaluation
#undef grammar_optimize_recreate
#undef grammar_optimize_round1
#undef grammar_optimize_round2
#undef grammar_optimize_round3
#undef grammar_optimize_round4
#undef grammar_optimize_shortcut_logical_evaluation
S */

using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

namespace
{
    const ParamSpec_ParamHolder plist_p[37] =
    {
    /* 0	*/ {2, 0, 0x0}, /* x */
    /* 1	*/ {2, 0, 0x4}, /* x */
    /* 2	*/ {2, Sign_Positive, 0x0}, /* x@P */
    /* 3	*/ {2, Sign_Negative | Constness_NotConst, 0x0}, /* x@N@V */
    /* 4	*/ {2, Sign_NoIdea, 0x0}, /* x */
    /* 5	*/ {2, Value_Logical, 0x0}, /* x@L */
    /* 6	*/ {3, Sign_NoIdea, 0x0}, /* y */
    /* 7	*/ {3, 0, 0x0}, /* y */
    /* 8	*/ {3, Value_Logical, 0x0}, /* y@L */
    /* 9	*/ {3, 0, 0x8}, /* y */
    /* 10	*/ {3, Value_OddInt, 0x0}, /* y@O */
    /* 11	*/ {3, Value_NonInteger, 0x0}, /* y@F */
    /* 12	*/ {3, Value_EvenInt, 0x0}, /* y@E */
    /* 13	*/ {3, Sign_Positive, 0x0}, /* y@P */
    /* 14	*/ {0, Sign_Negative | Constness_Const, 0x0}, /* %@N */
    /* 15	*/ {0, Constness_Const, 0x0}, /* % */
    /* 16	*/ {0, Sign_Positive | Constness_Const, 0x0}, /* %@P */
    /* 17	*/ {0, Value_EvenInt | Constness_Const, 0x0}, /* %@E */
    /* 18	*/ {0, Constness_Const, 0x1}, /* % */
    /* 19	*/ {0, Value_IsInteger | Sign_Positive | Constness_Const, 0x0}, /* %@I@P */
    /* 20	*/ {0, Oneness_NotOne | Constness_Const, 0x1}, /* %@M */
    /* 21	*/ {0, Oneness_NotOne | Constness_Const, 0x0}, /* %@M */
    /* 22	*/ {0, Oneness_One | Constness_Const, 0x0}, /* %@1 */
    /* 23	*/ {0, Value_Logical | Constness_Const, 0x0}, /* %@L */
    /* 24	*/ {1, Constness_Const, 0x0}, /* & */
    /* 25	*/ {1, Value_EvenInt | Constness_Const, 0x0}, /* &@E */
    /* 26	*/ {1, Oneness_NotOne | Constness_Const, 0x0}, /* &@M */
    /* 27	*/ {1, Value_IsInteger | Constness_Const, 0x0}, /* &@I */
    /* 28	*/ {1, Sign_Positive | Constness_Const, 0x0}, /* &@P */
    /* 29	*/ {1, Sign_Negative | Constness_Const, 0x0}, /* &@N */
    /* 30	*/ {6, 0, 0x0}, /* b */
    /* 31	*/ {4, 0, 0x0}, /* z */
    /* 32	*/ {4, Value_IsInteger, 0x0}, /* z@I */
    /* 33	*/ {4, Constness_Const, 0x0}, /* z@C */
    /* 34	*/ {4, 0, 0x16}, /* z */
    /* 35	*/ {5, 0, 0x0}, /* a */
    /* 36	*/ {5, Constness_Const, 0x0}, /* a@C */
    };

    template<typename Value_t>
    struct plist_n_container
    {
        static const ParamSpec_NumConstant<Value_t> plist_n[20];
    };
    template<typename Value_t>
    const ParamSpec_NumConstant <Value_t> plist_n_container<Value_t>::plist_n[20] =
    {
    /* 37	*/ {Value_t(-2), 0}, /* -2 */
    /* 38	*/ {Value_t(-1), 0}, /* -1 */
    /* 39	*/ {Value_t(-0.5), 0}, /* -0.5 */
    /* 40	*/ {Value_t(-0.25), 0}, /* -0.25 */
    /* 41	*/ {Value_t(0), 0}, /* 0 */
    /* 42	*/ {fp_const_deg_to_rad<Value_t>(), 0}, /* 0.0174532925199 */
    /* 43	*/ {fp_const_einv<Value_t>(), 0}, /* 0.367879441171 */
    /* 44	*/ {fp_const_log10inv<Value_t>(), 0}, /* 0.434294481903 */
    /* 45	*/ {Value_t(0.5), 0}, /* 0.5 */
    /* 46	*/ {fp_const_log2<Value_t>(), 0}, /* 0.69314718056 */
    /* 47	*/ {Value_t(1), 0}, /* 1 */
    /* 48	*/ {fp_const_log2inv<Value_t>(), 0}, /* 1.44269504089 */
    /* 49	*/ {Value_t(2), 0}, /* 2 */
    /* 50	*/ {fp_const_log10<Value_t>(), 0}, /* 2.30258509299 */
    /* 51	*/ {fp_const_e<Value_t>(), 0}, /* 2.71828182846 */
    /* 52	*/ {fp_const_rad_to_deg<Value_t>(), 0}, /* 57.2957795131 */
    /* 53	*/ {-fp_const_pihalf<Value_t>(), Modulo_Radians}, /* -1.57079632679 */
    /* 54	*/ {Value_t(0), Modulo_Radians}, /* 0 */
    /* 55	*/ {fp_const_pihalf<Value_t>(), Modulo_Radians}, /* 1.57079632679 */
    /* 56	*/ {fp_const_pi<Value_t>(), Modulo_Radians}, /* 3.14159265359 */
    };

    const ParamSpec_SubFunction plist_s[517] =
    {
    /* 57	*/ {{1,/*15         */15        , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* -%@C */
    /* 58	*/ {{1,/*398        */398       , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* -POW( MUL( % 0.5 )@C 2 )@C@C */
    /* 59	*/ {{1,/*477        */477       , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* -MIN( % & )@C@C */
    /* 60	*/ {{1,/*15         */15        , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x1}, /* -%@C */
    /* 61	*/ {{1,/*15         */15        , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /%@C */
    /* 62	*/ {{1,/*24         */24        , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /&@C */
    /* 63	*/ {{1,/*465        */465       , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /LOG( % )@C@C */
    /* 64	*/ {{1,/*466        */466       , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /LOG( & )@C@C */
    /* 65	*/ {{1,/*498        */498       , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /SQRT( % )@C@C */
    /* 66	*/ {{2,/*315,320    */327995    , cAdd        ,PositionalParams,0}, 0, 0x0}, /* (cAdd [(cPow [x %@E]) (cPow [y &@E])]) */
    /* 67	*/ {{2,/*148,47     */48276     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cPow [x 2]) -1}) 1}) */
    /* 68	*/ {{2,/*55,254     */260151    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1.57079632679 (cMul %@N <1>)}) */
    /* 69	*/ {{2,/*155,459    */470171    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {LOG( % )@C y}) (cLog [(cMul  <1>)])}) */
    /* 70	*/ {{2,/*166,165    */169126    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {a@C (cPow [x %@E])}) (cMul {z@C (cPow [y &@E])})}) */
    /* 71	*/ {{2,/*290,47     */48418     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x 2]) 1}) */
    /* 72	*/ {{2,/*304,1      */1328      , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [(cAdd {-1 (cPow [x 2])}) 0.5])@D4 x@D4}) */
    /* 73	*/ {{2,/*314,277    */283962    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x y]) MUL( % 0.5 )@C}) */
    /* 74	*/ {{2,/*315,165    */169275    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x %@E]) (cMul {z@C (cPow [y &@E])})}) */
    /* 75	*/ {{2,/*290,38     */39202     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x 2]) -1}) */
    /* 76	*/ {{2,/*316,277    */283964    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x (cMul {y &})]) MUL( % 0.5 )@C}) */
    /* 77	*/ {{2,/*325,277    */283973    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [& y]) MUL( % 0.5 )@C}) */
    /* 78	*/ {{2,/*459,465    */476619    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cLog [(cMul  <1>)]) LOG( % )@C}) */
    /* 79	*/ {{2,/*38,290     */296998    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {-1 (cPow [x 2])}) */
    /* 80	*/ {{2,/*47,0       */47        , cAdd        ,SelectedParams  ,0}, 0, 0x4}, /* (cAdd {1 x}) */
    /* 81	*/ {{2,/*47,158     */161839    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cMul {-1 x})}) */
    /* 82	*/ {{2,/*460,24     */25036     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cLog [x]) &}) */
    /* 83	*/ {{2,/*7,35       */35847     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {y a}) */
    /* 84	*/ {{2,/*24,59      */60440     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {& -MIN( % & )@C@C}) */
    /* 85	*/ {{2,/*31,30      */30751     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {z b}) */
    /* 86	*/ {{2,/*178,179    */183474    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {x SQRT( % )@C}) (cMul {y MUL( 0.5 MUL( & /SQRT( % )@C@C )@C )@C})}) */
    /* 87	*/ {{2,/*246,253    */259318    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul  <1>) (cMul -1 <2>)}) */
    /* 88	*/ {{2,/*263,264    */270599    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul (cPow [x (cAdd {% -MIN( % & )@C@C})]) <1>) (cMul (cPow [x (cAdd {& -MIN( % & )@C@C})]) <2>)}) */
    /* 89	*/ {{2,/*15,59      */60431     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {% -MIN( % & )@C@C}) */
    /* 90	*/ {{2,/*47,253     */259119    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cMul -1 <2>)}) */
    /* 91	*/ {{2,/*290,324    */332066    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x 2]) (cPow [y 2])}) */
    /* 92	*/ {{2,/*0,7        */7168      , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {x y}) */
    /* 93	*/ {{2,/*0,193      */197632    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {x (cMul {-1 y})}) */
    /* 94	*/ {{2,/*0,285      */291840    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {x MUL( % -0.5 )@C}) */
    /* 95	*/ {{2,/*0,277      */283648    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {x MUL( % 0.5 )@C}) */
    /* 96	*/ {{2,/*274,233    */238866    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul % & <1>) (cMul {& (cAdd  <2>)})}) */
    /* 97	*/ {{2,/*286,234    */239902    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {MUL( % & )@C (cMul {& (cAdd  <1>)})}) */
    /* 98	*/ {{2,/*7,31       */31751     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {y z}) */
    /* 99	*/ {{2,/*7,239      */244743    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})}) */
    /* 100	*/ {{2,/*22,375     */384022    , cAdd        ,SelectedParams  ,0}, 0, 0x4}, /* (cAdd {%@1 (cPow [x z])}) */
    /* 101	*/ {{2,/*238,376    */385262    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cPow [x y]) %}) (cPow [x (cAdd {y z})])}) */
    /* 102	*/ {{2,/*38,377     */386086    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {-1 (cPow [x@P z])}) */
    /* 103	*/ {{2,/*38,384     */393254    , cAdd        ,SelectedParams  ,0}, 0, 0x5}, /* (cAdd {-1 (cPow [% x])}) */
    /* 104	*/ {{2,/*38,384     */393254    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {-1 (cPow [% x])}) */
    /* 105	*/ {{2,/*47,377     */386095    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [x@P z])}) */
    /* 106	*/ {{2,/*240,378    */387312    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cPow [& y]) -1}) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})])}) */
    /* 107	*/ {{2,/*230,18     */18662     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {% (cPow [x@P z])})@D1 %@D1}) */
    /* 108	*/ {{2,/*230,60     */61670     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {% (cPow [x@P z])})@D1 -%@C@D1}) */
    /* 109	*/ {{2,/*325,378    */387397    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [& y]) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})])}) */
    /* 110	*/ {{2,/*47,242     */247855    , cAdd        ,SelectedParams  ,0}, 0, 0x1}, /* (cAdd {1 (cMul {(cLog [x]) /%@C})}) */
    /* 111	*/ {{2,/*47,334     */342063    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [(cAdd  <1>) 2])}) */
    /* 112	*/ {{2,/*47,290     */297007    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [x 2])}) */
    /* 113	*/ {{2,/*460,15     */15820     , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cLog [x]) %}) */
    /* 114	*/ {{2,/*47,384     */393263    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [% x])}) */
    /* 115	*/ {{2,/*47,384     */393263    , cAdd        ,SelectedParams  ,0}, 0, 0x5}, /* (cAdd {1 (cPow [% x])}) */
    /* 116	*/ {{2,/*55,158     */161847    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1.57079632679 (cMul {-1 x})}) */
    /* 117	*/ {{2,/*55,252     */258103    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1.57079632679 (cMul -1 <1>)}) */
    /* 118	*/ {{2,/*241,243    */249073    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cAbs [x]) -%@C}) (cMul {x (cAdd  <1>)})}) */
    /* 119	*/ {{2,/*244,243    */249076    , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cAbs [x]) %}) (cMul {x (cAdd  <1>)})}) */
    /* 120	*/ {{0,/*           */0         , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd  <1>) */
    /* 121	*/ {{0,/*           */0         , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd  <2>) */
    /* 122	*/ {{1,/*45         */45        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd 0.5 <1>) */
    /* 123	*/ {{1,/*53         */53        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd -1.57079632679 <1>) */
    /* 124	*/ {{1,/*54         */54        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd 0 <1>) */
    /* 125	*/ {{1,/*55         */55        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd 1.57079632679 <1>) */
    /* 126	*/ {{1,/*56         */56        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd 3.14159265359 <1>) */
    /* 127	*/ {{1,/*26         */26        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd &@M <1>) */
    /* 128	*/ {{1,/*259        */259       , cAdd        ,AnyParams       ,1}, 0, 0x16}, /* (cAdd (cMul (cPow [(cLog [z]) -1]) <2>) <1>) */
    /* 129	*/ {{1,/*272        */272       , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd (cMul %@M <1>) <2>) */
    /* 130	*/ {{1,/*323        */323       , cAdd        ,AnyParams       ,1}, 0, 0x16}, /* (cAdd (cPow [(cLog [z]) -1]) <1>) */
    /* 131	*/ {{1,/*0          */0         , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd x <1>) */
    /* 132	*/ {{1,/*21         */21        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd %@M <1>) */
    /* 133	*/ {{1,/*447        */447       , cAdd        ,AnyParams       ,1}, 0, 0x4}, /* (cAdd (cIf [(cLess [x 0]) %@D1 -%@C@D1]) <1>) */
    /* 134	*/ {{1,/*449        */449       , cAdd        ,AnyParams       ,1}, 0, 0x4}, /* (cAdd (cIf [(cGreater [x 0]) %@D1 -%@C@D1]) <1>) */
    /* 135	*/ {{1,/*0          */0         , cAdd        ,AnyParams       ,1}, 0, 0x4}, /* (cAdd x <1>) */
    /* 136	*/ {{1,/*0          */0         , cAdd        ,AnyParams       ,2}, 0, 0x4}, /* (cAdd x <2>) */
    /* 137	*/ {{1,/*15         */15        , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd % <1>) */
    /* 138	*/ {{1,/*24         */24        , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd & <2>) */
    /* 139	*/ {{2,/*24,57      */58392     , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd & -%@C <2>) */
    /* 140	*/ {{0,/*           */0         , cAdd        ,AnyParams       ,1}, Sign_Positive, 0x0}, /* (cAdd  <1>)@P */
    /* 141	*/ {{2,/*15,24      */24591     , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( % & )@C */
    /* 142	*/ {{2,/*15,33      */33807     , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( % z@C )@C */
    /* 143	*/ {{2,/*15,47      */48143     , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( % 1 )@C */
    /* 144	*/ {{2,/*24,279     */285720    , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( & MUL( -2 % )@C )@C */
    /* 145	*/ {{2,/*24,284     */290840    , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( & MUL( 2 % )@C )@C */
    /* 146	*/ {{2,/*0,298      */305152    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cPow [y %@N])}) */
    /* 147	*/ {{2,/*80,305     */312400    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cAdd {1 x})@D4 (cPow [(cAdd {1 (cMul {-1 x})}) -1])@D4}) */
    /* 148	*/ {{2,/*290,38     */39202     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x 2]) -1}) */
    /* 149	*/ {{2,/*38,120     */122918    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cAdd  <1>)}) */
    /* 150	*/ {{2,/*38,412     */421926    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCeil [(cMul  <1>)])}) */
    /* 151	*/ {{2,/*38,419     */429094    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCos [(cAdd  <1>)])}) */
    /* 152	*/ {{2,/*38,433     */443430    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cFloor [(cMul  <1>)])}) */
    /* 153	*/ {{2,/*394,310    */317834    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {POW( % & )@C (cPow [(cMul  <1>) &])}) */
    /* 154	*/ {{2,/*394,321    */329098    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {POW( % & )@C (cPow [% (cAdd  <1>)])}) */
    /* 155	*/ {{2,/*465,7      */7633      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {LOG( % )@C y}) */
    /* 156	*/ {{2,/*538,7      */7706      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cNot [x]) y}) */
    /* 157	*/ {{2,/*562,7      */7730      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cNotNot [x]) y}) */
    /* 158	*/ {{2,/*38,0       */38        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 x}) */
    /* 159	*/ {{2,/*411,49     */50587     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cAtanh [x]) 2}) */
    /* 160	*/ {{2,/*0,397      */406528    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x POW( a@C /%@C )@C}) */
    /* 161	*/ {{2,/*7,24       */24583     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y &}) */
    /* 162	*/ {{2,/*7,31       */31751     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y z}) */
    /* 163	*/ {{2,/*7,396      */405511    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y POW( z@C /&@C )@C}) */
    /* 164	*/ {{2,/*15,314     */321551    , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {% (cPow [x y])}) */
    /* 165	*/ {{2,/*33,320     */327713    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {z@C (cPow [y &@E])}) */
    /* 166	*/ {{2,/*36,315     */322596    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {a@C (cPow [x %@E])}) */
    /* 167	*/ {{2,/*297,88     */90409     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x MIN( % & )@C]) (cAdd {(cMul (cPow [x (cAdd {% -MIN( % & )@C@C})]) <1>) (cMul (cPow [x (cAdd {& -MIN( % & )@C@C})]) <2>)})}) */
    /* 168	*/ {{2,/*326,327    */335174    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [z (cAdd  <1>)]) (cPow [2.71828182846 (cMul  <2>)])}) */
    /* 169	*/ {{2,/*394,319    */327050    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {POW( % & )@C (cPow [x LOG( % )@C])}) */
    /* 170	*/ {{2,/*38,482     */493606    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSin [(cAdd  <1>)])}) */
    /* 171	*/ {{2,/*38,485     */496678    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSin [(cMul  <1>)])}) */
    /* 172	*/ {{2,/*38,492     */503846    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSinh [(cMul  <1>)])}) */
    /* 173	*/ {{2,/*38,504     */516134    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cTan [(cMul  <1>)])}) */
    /* 174	*/ {{2,/*49,7       */7217      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {2 y}) */
    /* 175	*/ {{2,/*51,326     */333875    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {2.71828182846 (cPow [z (cAdd  <1>)])}) */
    /* 176	*/ {{2,/*0,329      */336896    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cPow [y -1])}) */
    /* 177	*/ {{2,/*38,512     */524326    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cTanh [(cMul  <1>)])}) */
    /* 178	*/ {{2,/*0,498      */509952    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x SQRT( % )@C}) */
    /* 179	*/ {{2,/*7,280      */286727    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y MUL( 0.5 MUL( & /SQRT( % )@C@C )@C )@C}) */
    /* 180	*/ {{2,/*15,87      */89103     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% (cAdd {(cMul  <1>) (cMul -1 <2>)})}) */
    /* 181	*/ {{2,/*15,90      */92175     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% (cAdd {1 (cMul -1 <2>)})}) */
    /* 182	*/ {{2,/*16,290     */296976    , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {%@P (cPow [x 2])}) */
    /* 183	*/ {{2,/*15,317     */324623    , cMul        ,SelectedParams  ,0}, 0, 0x14}, /* (cMul {% (cPow [x (cMul {& y})])}) */
    /* 184	*/ {{2,/*15,325     */332815    , cMul        ,SelectedParams  ,0}, 0, 0x10}, /* (cMul {% (cPow [& y])}) */
    /* 185	*/ {{3,/*24,0,7     */7340056   , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {& x y}) */
    /* 186	*/ {{2,/*324,282    */289092    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [y 2]) MUL( MUL( -0.25 /%@C )@C POW( & 2 )@C )@C}) */
    /* 187	*/ {{2,/*16,91      */93200     , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {%@P (cAdd {(cPow [x 2]) (cPow [y 2])})}) */
    /* 188	*/ {{2,/*15,330     */337935    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% (cPow [(cAdd {x y}) 2])}) */
    /* 189	*/ {{3,/*28,0,7     */7340060   , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {&@P x y}) */
    /* 190	*/ {{3,/*144,0,7    */7340176   , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {ADD( & MUL( -2 % )@C )@C x y}) */
    /* 191	*/ {{2,/*15,331     */338959    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% (cPow [(cAdd {x (cMul {-1 y})}) 2])}) */
    /* 192	*/ {{3,/*29,0,7     */7340061   , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {&@N x y}) */
    /* 193	*/ {{2,/*38,7       */7206      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 y}) */
    /* 194	*/ {{2,/*0,7        */7168      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x y}) */
    /* 195	*/ {{2,/*38,349     */357414    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cPow [(cCos [x]) 2])}) */
    /* 196	*/ {{2,/*38,360     */368678    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cPow [(cSin [x]) 2])}) */
    /* 197	*/ {{2,/*57,362     */370745    , cMul        ,SelectedParams  ,0}, 0, 0x7}, /* (cMul {-%@C (cPow [/&@C x])}) */
    /* 198	*/ {{3,/*145,0,7    */7340177   , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {ADD( & MUL( 2 % )@C )@C x y}) */
    /* 199	*/ {{2,/*365,38     */39277     , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cPow [0.367879441171 x]) -1}) */
    /* 200	*/ {{2,/*414,416    */426398    , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cCos [y])}) */
    /* 201	*/ {{3,/*414,416,38 */40272286  , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cCos [y]) -1}) */
    /* 202	*/ {{2,/*414,479    */490910    , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cSin [y])}) */
    /* 203	*/ {{3,/*414,479,38 */40336798  , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cSin [y]) -1}) */
    /* 204	*/ {{2,/*424,49     */50600     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cCosh [(cLog [(cPow [& x])])]) 2}) */
    /* 205	*/ {{2,/*478,416    */426462    , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cSin [x]) (cCos [y])}) */
    /* 206	*/ {{2,/*478,479    */490974    , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cSin [x]) (cSin [y])}) */
    /* 207	*/ {{2,/*38,362     */370726    , cMul        ,SelectedParams  ,0}, 0, 0x6}, /* (cMul {-1 (cPow [/&@C x])}) */
    /* 208	*/ {{2,/*38,363     */371750    , cMul        ,SelectedParams  ,0}, 0, 0x6}, /* (cMul {-1 (cPow [& x])}) */
    /* 209	*/ {{2,/*38,418     */428070    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCos [(cAdd {x (cMul {-1 y})})])}) */
    /* 210	*/ {{3,/*478,479,38 */40336862  , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cSin [x]) (cSin [y]) -1}) */
    /* 211	*/ {{2,/*490,37     */38378     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cSinh [(cMul {x LOG( & )@C})]) -2}) */
    /* 212	*/ {{2,/*495,49     */50671     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cSinh [(cLog [(cPow [& x])])]) 2}) */
    /* 213	*/ {{3,/*0,465,45   */47662080  , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x LOG( % )@C 0.5}) */
    /* 214	*/ {{2,/*0,466      */477184    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x LOG( & )@C}) */
    /* 215	*/ {{2,/*0,555      */568320    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cAnd  <1>)}) */
    /* 216	*/ {{2,/*15,363     */371727    , cMul        ,SelectedParams  ,0}, 0, 0x7}, /* (cMul {% (cPow [& x])}) */
    /* 217	*/ {{3,/*490,49,15  */15779306  , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cSinh [(cMul {x LOG( & )@C})]) 2 %}) */
    /* 218	*/ {{2,/*15,362     */370703    , cMul        ,SelectedParams  ,0}, 0, 0x7}, /* (cMul {% (cPow [/&@C x])}) */
    /* 219	*/ {{2,/*365,38     */39277     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [0.367879441171 x]) -1}) */
    /* 220	*/ {{2,/*367,38     */39279     , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cPow [2.71828182846 x]) -1}) */
    /* 221	*/ {{3,/*422,49,15  */15779238  , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cCosh [(cMul {x LOG( & )@C})]) 2 %}) */
    /* 222	*/ {{2,/*426,38     */39338     , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cCosh [x]) -1}) */
    /* 223	*/ {{2,/*38,426     */436262    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCosh [x])}) */
    /* 224	*/ {{2,/*38,497     */508966    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSinh [x])}) */
    /* 225	*/ {{2,/*497,38     */39409     , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cSinh [x]) -1}) */
    /* 226	*/ {{2,/*38,290     */296998    , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {-1 (cPow [x 2])}) */
    /* 227	*/ {{2,/*7,35       */35847     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y a}) */
    /* 228	*/ {{2,/*15,0       */15        , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {% x}) */
    /* 229	*/ {{2,/*38,369     */377894    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cPow [(cAdd {x MUL( % -0.5 )@C}) 2])}) */
    /* 230	*/ {{2,/*15,377     */386063    , cMul        ,SelectedParams  ,0}, 0, 0x1}, /* (cMul {% (cPow [x@P z])}) */
    /* 231	*/ {{2,/*15,0       */15        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% x}) */
    /* 232	*/ {{2,/*24,7       */7192      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {& y}) */
    /* 233	*/ {{2,/*24,121     */123928    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {& (cAdd  <2>)}) */
    /* 234	*/ {{2,/*24,120     */122904    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {& (cAdd  <1>)}) */
    /* 235	*/ {{2,/*31,30      */30751     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {z b}) */
    /* 236	*/ {{2,/*57,0       */57        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-%@C x}) */
    /* 237	*/ {{2,/*288,7      */7456      , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {MUL( & 2 )@C y}) */
    /* 238	*/ {{2,/*314,15     */15674     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x y]) %}) */
    /* 239	*/ {{3,/*31,460,64  */67579935  , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {z (cLog [x]) /LOG( & )@C@C}) */
    /* 240	*/ {{2,/*325,38     */39237     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [& y]) -1}) */
    /* 241	*/ {{2,/*400,57     */58768     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cAbs [x]) -%@C}) */
    /* 242	*/ {{2,/*460,61     */62924     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cLog [x]) /%@C}) */
    /* 243	*/ {{2,/*0,120      */122880    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cAdd  <1>)}) */
    /* 244	*/ {{2,/*400,15     */15760     , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cAbs [x]) %}) */
    /* 245	*/ {{3,/*0,45,61    */64009216  , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x 0.5 /%@C}) */
    /* 246	*/ {{0,/*           */0         , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul  <1>) */
    /* 247	*/ {{0,/*           */0         , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul  <2>) */
    /* 248	*/ {{1,/*2          */2         , cMul        ,AnyParams       ,1}, 0, 0x4}, /* (cMul x@P <1>) */
    /* 249	*/ {{1,/*2          */2         , cMul        ,AnyParams       ,2}, 0, 0x4}, /* (cMul x@P <2>) */
    /* 250	*/ {{1,/*3          */3         , cMul        ,AnyParams       ,1}, 0, 0x4}, /* (cMul x@N@V <1>) */
    /* 251	*/ {{1,/*3          */3         , cMul        ,AnyParams       ,2}, 0, 0x4}, /* (cMul x@N@V <2>) */
    /* 252	*/ {{1,/*38         */38        , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul -1 <1>) */
    /* 253	*/ {{1,/*38         */38        , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul -1 <2>) */
    /* 254	*/ {{1,/*14         */14        , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul %@N <1>) */
    /* 255	*/ {{1,/*57         */57        , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul -%@C <1>) */
    /* 256	*/ {{1,/*16         */16        , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul %@P <1>) */
    /* 257	*/ {{2,/*63,460     */471103    , cMul        ,AnyParams       ,1}, 0, 0x1}, /* (cMul /LOG( % )@C@C (cLog [x]) <1>) */
    /* 258	*/ {{1,/*303        */303       , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cPow [% y]) <1>) */
    /* 259	*/ {{1,/*323        */323       , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul (cPow [(cLog [z]) -1]) <2>) */
    /* 260	*/ {{2,/*323,460    */471363    , cMul        ,AnyParams       ,1}, 0, 0x16}, /* (cMul (cPow [(cLog [z]) -1]) (cLog [x]) <1>) */
    /* 261	*/ {{1,/*293        */293       , cMul        ,AnyParams       ,1}, 0, 0x4}, /* (cMul (cPow [x %@I@P]) <1>) */
    /* 262	*/ {{1,/*294        */294       , cMul        ,AnyParams       ,2}, 0, 0x4}, /* (cMul (cPow [x &@I]) <2>) */
    /* 263	*/ {{1,/*295        */295       , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cPow [x (cAdd {% -MIN( % & )@C@C})]) <1>) */
    /* 264	*/ {{1,/*296        */296       , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul (cPow [x (cAdd {& -MIN( % & )@C@C})]) <2>) */
    /* 265	*/ {{1,/*400        */400       , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cAbs [x]) <1>) */
    /* 266	*/ {{1,/*0          */0         , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul x <1>) */
    /* 267	*/ {{1,/*460        */460       , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cLog [x]) <1>) */
    /* 268	*/ {{1,/*465        */465       , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul LOG( % )@C <1>) */
    /* 269	*/ {{1,/*16         */16        , cMul        ,AnyParams       ,1}, 0, 0x1}, /* (cMul %@P <1>) */
    /* 270	*/ {{1,/*57         */57        , cMul        ,AnyParams       ,2}, 0, 0x1}, /* (cMul -%@C <2>) */
    /* 271	*/ {{1,/*0          */0         , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul x <2>) */
    /* 272	*/ {{1,/*21         */21        , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul %@M <1>) */
    /* 273	*/ {{1,/*15         */15        , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul % <1>) */
    /* 274	*/ {{2,/*15,24      */24591     , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul % & <1>) */
    /* 275	*/ {{1,/*24         */24        , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul & <2>) */
    /* 276	*/ {{1,/*517        */517       , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul DIV( & % )@C <2>) */
    /* 277	*/ {{2,/*15,45      */46095     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( % 0.5 )@C */
    /* 278	*/ {{2,/*24,45      */46104     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( & 0.5 )@C */
    /* 279	*/ {{2,/*37,15      */15397     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( -2 % )@C */
    /* 280	*/ {{2,/*45,281     */287789    , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( 0.5 MUL( & /SQRT( % )@C@C )@C )@C */
    /* 281	*/ {{2,/*24,65      */66584     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( & /SQRT( % )@C@C )@C */
    /* 282	*/ {{2,/*283,395    */404763    , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( MUL( -0.25 /%@C )@C POW( & 2 )@C )@C */
    /* 283	*/ {{2,/*40,61      */62504     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( -0.25 /%@C )@C */
    /* 284	*/ {{2,/*49,15      */15409     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( 2 % )@C */
    /* 285	*/ {{2,/*15,39      */39951     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( % -0.5 )@C */
    /* 286	*/ {{2,/*15,24      */24591     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( % & )@C */
    /* 287	*/ {{2,/*15,33      */33807     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( % z@C )@C */
    /* 288	*/ {{2,/*24,49      */50200     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( & 2 )@C */
    /* 289	*/ {{2,/*45,61      */62509     , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( 0.5 /%@C )@C */
    /* 290	*/ {{2,/*0,49       */50176     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x 2]) */
    /* 291	*/ {{2,/*0,174      */178176    , cPow        ,PositionalParams,0}, 0, 0x12}, /* (cPow [x (cMul {2 y})]) */
    /* 292	*/ {{2,/*0,277      */283648    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x MUL( % 0.5 )@C]) */
    /* 293	*/ {{2,/*0,19       */19456     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %@I@P]) */
    /* 294	*/ {{2,/*0,27       */27648     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x &@I]) */
    /* 295	*/ {{2,/*0,89       */91136     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {% -MIN( % & )@C@C})]) */
    /* 296	*/ {{2,/*0,84       */86016     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {& -MIN( % & )@C@C})]) */
    /* 297	*/ {{2,/*0,477      */488448    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x MIN( % & )@C]) */
    /* 298	*/ {{2,/*6,14       */14342     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y %@N]) */
    /* 299	*/ {{2,/*7,57       */58375     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y -%@C]) */
    /* 300	*/ {{2,/*67,45      */46147     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {(cMul {(cPow [x 2]) -1}) 1}) 0.5]) */
    /* 301	*/ {{2,/*71,45      */46151     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [x 2]) 1}) 0.5]) */
    /* 302	*/ {{2,/*7,278      */284679    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y MUL( & 0.5 )@C]) */
    /* 303	*/ {{2,/*15,7       */7183      , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [% y]) */
    /* 304	*/ {{2,/*79,45      */46159     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {-1 (cPow [x 2])}) 0.5]) */
    /* 305	*/ {{2,/*81,38      */38993     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {1 (cMul {-1 x})}) -1]) */
    /* 306	*/ {{2,/*86,49      */50262     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cMul {x SQRT( % )@C}) (cMul {y MUL( 0.5 MUL( & /SQRT( % )@C@C )@C )@C})}) 2]) */
    /* 307	*/ {{2,/*73,49      */50249     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [x y]) MUL( % 0.5 )@C}) 2]) */
    /* 308	*/ {{2,/*160,277    */283808    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cMul {x POW( a@C /%@C )@C}) MUL( % 0.5 )@C]) */
    /* 309	*/ {{2,/*163,278    */284835    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cMul {y POW( z@C /&@C )@C}) MUL( & 0.5 )@C]) */
    /* 310	*/ {{2,/*246,24     */24822     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cMul  <1>) &]) */
    /* 311	*/ {{2,/*0,10       */10240     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x y@O]) */
    /* 312	*/ {{2,/*0,11       */11264     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x y@F]) */
    /* 313	*/ {{2,/*2,7        */7170      , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x@P y]) */
    /* 314	*/ {{2,/*0,7        */7168      , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x y]) */
    /* 315	*/ {{2,/*0,17       */17408     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %@E]) */
    /* 316	*/ {{2,/*0,161      */164864    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cMul {y &})]) */
    /* 317	*/ {{2,/*0,232      */237568    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cMul {& y})]) */
    /* 318	*/ {{2,/*0,237      */242688    , cPow        ,PositionalParams,0}, 0, 0x14}, /* (cPow [x (cMul {MUL( & 2 )@C y})]) */
    /* 319	*/ {{2,/*0,465      */476160    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x LOG( % )@C]) */
    /* 320	*/ {{2,/*7,25       */25607     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y &@E]) */
    /* 321	*/ {{2,/*15,120     */122895    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [% (cAdd  <1>)]) */
    /* 322	*/ {{2,/*76,49      */50252     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [x (cMul {y &})]) MUL( % 0.5 )@C}) 2]) */
    /* 323	*/ {{2,/*462,38     */39374     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog [z]) -1]) */
    /* 324	*/ {{2,/*7,49       */50183     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y 2]) */
    /* 325	*/ {{2,/*24,7       */7192      , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [& y]) */
    /* 326	*/ {{2,/*31,120     */122911    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [z (cAdd  <1>)]) */
    /* 327	*/ {{2,/*51,247     */252979    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 (cMul  <2>)]) */
    /* 328	*/ {{2,/*75,45      */46155     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [x 2]) -1}) 0.5]) */
    /* 329	*/ {{2,/*7,38       */38919     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y -1]) */
    /* 330	*/ {{2,/*92,49      */50268     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {x y}) 2]) */
    /* 331	*/ {{2,/*93,49      */50269     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {x (cMul {-1 y})}) 2]) */
    /* 332	*/ {{2,/*77,49      */50253     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [& y]) MUL( % 0.5 )@C}) 2]) */
    /* 333	*/ {{2,/*111,45     */46191     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) 0.5]) */
    /* 334	*/ {{2,/*120,49     */50296     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd  <1>) 2]) */
    /* 335	*/ {{2,/*395,7      */7563      , cPow        ,PositionalParams,0}, 0, 0x10}, /* (cPow [POW( & 2 )@C y]) */
    /* 336	*/ {{2,/*43,407     */416811    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [0.367879441171 (cAsinh [(cAdd  <1>)])]) */
    /* 337	*/ {{2,/*51,407     */416819    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 (cAsinh [(cAdd  <1>)])]) */
    /* 338	*/ {{2,/*111,39     */40047     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) -0.5]) */
    /* 339	*/ {{2,/*112,45     */46192     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {1 (cPow [x 2])}) 0.5]) */
    /* 340	*/ {{2,/*51,406     */415795    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 (cAsinh [x])]) */
    /* 341	*/ {{2,/*112,39     */40048     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {1 (cPow [x 2])}) -0.5]) */
    /* 342	*/ {{2,/*43,406     */415787    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [0.367879441171 (cAsinh [x])]) */
    /* 343	*/ {{2,/*104,38     */39016     , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cAdd {-1 (cPow [% x])}) -1]) */
    /* 344	*/ {{2,/*414,38     */39326     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCos [x]) -1]) */
    /* 345	*/ {{2,/*414,38     */39326     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCos [x]) -1]) */
    /* 346	*/ {{2,/*420,38     */39332     , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cCos [(cMul {-%@C x})]) -1]) */
    /* 347	*/ {{2,/*421,38     */39333     , cPow        ,PositionalParams,0}, 0, 0x1}, /* (cPow [(cCos [(cMul -%@C <1>)]) -1]) */
    /* 348	*/ {{2,/*414,49     */50590     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCos [x]) 2]) */
    /* 349	*/ {{2,/*414,49     */50590     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCos [x]) 2]) */
    /* 350	*/ {{2,/*426,38     */39338     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCosh [x]) -1]) */
    /* 351	*/ {{2,/*426,38     */39338     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCosh [x]) -1]) */
    /* 352	*/ {{2,/*423,38     */39335     , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cCosh [(cMul {-%@C x})]) -1]) */
    /* 353	*/ {{2,/*426,15     */15786     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCosh [x]) %]) */
    /* 354	*/ {{2,/*426,143    */146858    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCosh [x]) ADD( % 1 )@C]) */
    /* 355	*/ {{2,/*460,38     */39372     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog [x]) -1]) */
    /* 356	*/ {{2,/*467,38     */39379     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog10 [x]) -1]) */
    /* 357	*/ {{2,/*468,38     */39380     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog2 [x]) -1]) */
    /* 358	*/ {{2,/*478,38     */39390     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cSin [x]) -1]) */
    /* 359	*/ {{2,/*478,49     */50654     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cSin [x]) 2]) */
    /* 360	*/ {{2,/*478,49     */50654     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cSin [x]) 2]) */
    /* 361	*/ {{2,/*24,0       */24        , cPow        ,PositionalParams,0}, 0, 0x6}, /* (cPow [& x]) */
    /* 362	*/ {{2,/*62,0       */62        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [/&@C x]) */
    /* 363	*/ {{2,/*24,0       */24        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [& x]) */
    /* 364	*/ {{2,/*62,0       */62        , cPow        ,PositionalParams,0}, 0, 0x6}, /* (cPow [/&@C x]) */
    /* 365	*/ {{2,/*43,0       */43        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [0.367879441171 x]) */
    /* 366	*/ {{2,/*43,0       */43        , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [0.367879441171 x]) */
    /* 367	*/ {{2,/*51,0       */51        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 x]) */
    /* 368	*/ {{2,/*51,0       */51        , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [2.71828182846 x]) */
    /* 369	*/ {{2,/*94,49      */50270     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {x MUL( % -0.5 )@C}) 2]) */
    /* 370	*/ {{2,/*0,49       */50176     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [x 2]) */
    /* 371	*/ {{2,/*95,49      */50271     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {x MUL( % 0.5 )@C}) 2]) */
    /* 372	*/ {{2,/*247,38     */39159     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cMul  <2>) -1]) */
    /* 373	*/ {{2,/*271,38     */39183     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cMul x <2>) -1]) */
    /* 374	*/ {{2,/*0,7        */7168      , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [x y]) */
    /* 375	*/ {{2,/*0,31       */31744     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x z]) */
    /* 376	*/ {{2,/*0,98       */100352    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {y z})]) */
    /* 377	*/ {{2,/*2,31       */31746     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x@P z]) */
    /* 378	*/ {{2,/*24,99      */101400    , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})]) */
    /* 379	*/ {{2,/*497,38     */39409     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cSinh [x]) -1]) */
    /* 380	*/ {{2,/*499,38     */39411     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cTan [x]) -1]) */
    /* 381	*/ {{2,/*499,38     */39411     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cTan [x]) -1]) */
    /* 382	*/ {{2,/*508,38     */39420     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cTanh [x]) -1]) */
    /* 383	*/ {{2,/*508,38     */39420     , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cTanh [x]) -1]) */
    /* 384	*/ {{2,/*15,0       */15        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [% x]) */
    /* 385	*/ {{2,/*114,38     */39026     , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cAdd {1 (cPow [% x])}) -1]) */
    /* 386	*/ {{2,/*510,38     */39422     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cTanh [(cMul {x LOG( % )@C 0.5})]) -1]) */
    /* 387	*/ {{2,/*0,16       */16384     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %@P]) */
    /* 388	*/ {{2,/*389,61     */62853     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cPow [x %]) /%@C]) */
    /* 389	*/ {{2,/*0,15       */15360     , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %]) */
    /* 390	*/ {{2,/*15,0       */15        , cPow        ,PositionalParams,0}, 0, 0x1}, /* (cPow [% x]) */
    /* 391	*/ {{2,/*16,0       */16        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [%@P x]) */
    /* 392	*/ {{2,/*15,7       */7183      , cPow        ,PositionalParams,0}, 0, 0x1}, /* (cPow [% y]) */
    /* 393	*/ {{2,/*4,7        */7172      , cPow        ,PositionalParams,0}, Sign_Positive, 0x0}, /* (cPow [x y])@P */
    /* 394	*/ {{2,/*15,24      */24591     , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( % & )@C */
    /* 395	*/ {{2,/*24,49      */50200     , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( & 2 )@C */
    /* 396	*/ {{2,/*33,62      */63521     , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( z@C /&@C )@C */
    /* 397	*/ {{2,/*36,61      */62500     , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( a@C /%@C )@C */
    /* 398	*/ {{2,/*277,49     */50453     , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( MUL( % 0.5 )@C 2 )@C */
    /* 399	*/ {{2,/*24,61      */62488     , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( & /%@C )@C */
    /* 400	*/ {{1,/*0          */0         , cAbs        ,PositionalParams,0}, 0, 0x0}, /* (cAbs [x]) */
    /* 401	*/ {{1,/*7          */7         , cAbs        ,PositionalParams,0}, 0, 0x0}, /* (cAbs [y]) */
    /* 402	*/ {{1,/*194        */194       , cAbs        ,PositionalParams,0}, 0, 0x0}, /* (cAbs [(cMul {x y})]) */
    /* 403	*/ {{1,/*0          */0         , cAcos       ,PositionalParams,0}, 0, 0x0}, /* (cAcos [x]) */
    /* 404	*/ {{1,/*0          */0         , cAcosh      ,PositionalParams,0}, 0, 0x0}, /* (cAcosh [x]) */
    /* 405	*/ {{1,/*0          */0         , cAsin       ,PositionalParams,0}, 0, 0x0}, /* (cAsin [x]) */
    /* 406	*/ {{1,/*0          */0         , cAsinh      ,PositionalParams,0}, 0, 0x0}, /* (cAsinh [x]) */
    /* 407	*/ {{1,/*120        */120       , cAsinh      ,PositionalParams,0}, 0, 0x0}, /* (cAsinh [(cAdd  <1>)]) */
    /* 408	*/ {{1,/*0          */0         , cAtan       ,PositionalParams,0}, 0, 0x0}, /* (cAtan [x]) */
    /* 409	*/ {{2,/*0,299      */306176    , cAtan2      ,PositionalParams,0}, 0, 0x0}, /* (cAtan2 [x (cPow [y -%@C])]) */
    /* 410	*/ {{2,/*0,7        */7168      , cAtan2      ,PositionalParams,0}, 0, 0x0}, /* (cAtan2 [x y]) */
    /* 411	*/ {{1,/*0          */0         , cAtanh      ,PositionalParams,0}, 0, 0x0}, /* (cAtanh [x]) */
    /* 412	*/ {{1,/*246        */246       , cCeil       ,PositionalParams,0}, 0, 0x0}, /* (cCeil [(cMul  <1>)]) */
    /* 413	*/ {{1,/*0          */0         , cCeil       ,PositionalParams,0}, 0, 0x4}, /* (cCeil [x]) */
    /* 414	*/ {{1,/*0          */0         , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [x]) */
    /* 415	*/ {{1,/*0          */0         , cCos        ,PositionalParams,0}, 0, 0x4}, /* (cCos [x]) */
    /* 416	*/ {{1,/*7          */7         , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [y]) */
    /* 417	*/ {{1,/*92         */92        , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cAdd {x y})]) */
    /* 418	*/ {{1,/*93         */93        , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cAdd {x (cMul {-1 y})})]) */
    /* 419	*/ {{1,/*120        */120       , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cAdd  <1>)]) */
    /* 420	*/ {{1,/*236        */236       , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cMul {-%@C x})]) */
    /* 421	*/ {{1,/*255        */255       , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cMul -%@C <1>)]) */
    /* 422	*/ {{1,/*214        */214       , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [(cMul {x LOG( & )@C})]) */
    /* 423	*/ {{1,/*236        */236       , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [(cMul {-%@C x})]) */
    /* 424	*/ {{1,/*464        */464       , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [(cLog [(cPow [& x])])]) */
    /* 425	*/ {{1,/*0          */0         , cCosh       ,PositionalParams,0}, 0, 0x4}, /* (cCosh [x]) */
    /* 426	*/ {{1,/*0          */0         , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [x]) */
    /* 427	*/ {{1,/*0          */0         , cExp        ,PositionalParams,0}, 0, 0x0}, /* (cExp [x]) */
    /* 428	*/ {{1,/*7          */7         , cExp        ,PositionalParams,0}, 0, 0x0}, /* (cExp [y]) */
    /* 429	*/ {{1,/*92         */92        , cExp        ,PositionalParams,0}, 0, 0x0}, /* (cExp [(cAdd {x y})]) */
    /* 430	*/ {{1,/*0          */0         , cExp2       ,PositionalParams,0}, 0, 0x0}, /* (cExp2 [x]) */
    /* 431	*/ {{1,/*7          */7         , cExp2       ,PositionalParams,0}, 0, 0x0}, /* (cExp2 [y]) */
    /* 432	*/ {{1,/*92         */92        , cExp2       ,PositionalParams,0}, 0, 0x0}, /* (cExp2 [(cAdd {x y})]) */
    /* 433	*/ {{1,/*246        */246       , cFloor      ,PositionalParams,0}, 0, 0x0}, /* (cFloor [(cMul  <1>)]) */
    /* 434	*/ {{1,/*0          */0         , cFloor      ,PositionalParams,0}, 0, 0x4}, /* (cFloor [x]) */
    /* 435	*/ {{2,/*292,302    */309540    , cHypot      ,PositionalParams,0}, 0, 0x0}, /* (cHypot [(cPow [x MUL( % 0.5 )@C]) (cPow [y MUL( & 0.5 )@C])]) */
    /* 436	*/ {{2,/*292,309    */316708    , cHypot      ,PositionalParams,0}, 0, 0x0}, /* (cHypot [(cPow [x MUL( % 0.5 )@C]) (cPow [(cMul {y POW( z@C /&@C )@C}) MUL( & 0.5 )@C])]) */
    /* 437	*/ {{2,/*308,309    */316724    , cHypot      ,PositionalParams,0}, 0, 0x0}, /* (cHypot [(cPow [(cMul {x POW( a@C /%@C )@C}) MUL( % 0.5 )@C]) (cPow [(cMul {y POW( z@C /&@C )@C}) MUL( & 0.5 )@C])]) */
    /* 438	*/ {{3,/*0,7,31     */32513024  , cIf         ,PositionalParams,0}, 0, 0x4}, /* (cIf [x y z]) */
    /* 439	*/ {{3,/*0,24,33    */34627584  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x & z@C]) */
    /* 440	*/ {{3,/*0,35,30    */31493120  , cIf         ,PositionalParams,0}, 0, 0x4}, /* (cIf [x a b]) */
    /* 441	*/ {{3,/*0,83,85    */89213952  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cAdd {y a}) (cAdd {z b})]) */
    /* 442	*/ {{3,/*0,141,142  */149042176 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x ADD( % & )@C ADD( % z@C )@C]) */
    /* 443	*/ {{3,/*0,227,235  */246647808 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cMul {y a}) (cMul {z b})]) */
    /* 444	*/ {{3,/*0,286,287  */301234176 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x MUL( % & )@C MUL( % z@C )@C]) */
    /* 445	*/ {{3,/*0,470,471  */494360576 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cMax [y a]) (cMax [z b])]) */
    /* 446	*/ {{3,/*0,474,475  */498558976 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cMin [y a]) (cMin [z b])]) */
    /* 447	*/ {{3,/*528,18,60  */62933520  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [(cLess [x 0]) %@D1 -%@C@D1]) */
    /* 448	*/ {{3,/*528,18,60  */62933520  , cIf         ,PositionalParams,0}, 0, 0x4}, /* (cIf [(cLess [x 0]) %@D1 -%@C@D1]) */
    /* 449	*/ {{3,/*534,18,60  */62933526  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [(cGreater [x 0]) %@D1 -%@C@D1]) */
    /* 450	*/ {{3,/*534,18,60  */62933526  , cIf         ,PositionalParams,0}, 0, 0x4}, /* (cIf [(cGreater [x 0]) %@D1 -%@C@D1]) */
    /* 451	*/ {{3,/*0,540,23   */24670208  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cNot [y]) %@L]) */
    /* 452	*/ {{3,/*0,551,552  */579378176 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cAnd {y a}) (cAnd {z b})]) */
    /* 453	*/ {{3,/*0,7,547    */573578240 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x y (cNot [%])]) */
    /* 454	*/ {{3,/*0,7,31     */32513024  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x y z]) */
    /* 455	*/ {{3,/*0,23,540   */566254592 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x %@L (cNot [y])]) */
    /* 456	*/ {{3,/*0,547,7    */7900160   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cNot [%]) y]) */
    /* 457	*/ {{3,/*0,558,561  */588822528 , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cOr {y a}) (cOr {z b})]) */
    /* 458	*/ {{1,/*120        */120       , cInt        ,PositionalParams,0}, 0, 0x0}, /* (cInt [(cAdd  <1>)]) */
    /* 459	*/ {{1,/*246        */246       , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [(cMul  <1>)]) */
    /* 460	*/ {{1,/*0          */0         , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [x]) */
    /* 461	*/ {{1,/*7          */7         , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [y]) */
    /* 462	*/ {{1,/*31         */31        , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [z]) */
    /* 463	*/ {{1,/*194        */194       , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [(cMul {x y})]) */
    /* 464	*/ {{1,/*363        */363       , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [(cPow [& x])]) */
    /* 465	*/ {{1,/*15         */15        , cLog        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* LOG( % )@C */
    /* 466	*/ {{1,/*24         */24        , cLog        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* LOG( & )@C */
    /* 467	*/ {{1,/*0          */0         , cLog10      ,PositionalParams,0}, 0, 0x0}, /* (cLog10 [x]) */
    /* 468	*/ {{1,/*0          */0         , cLog2       ,PositionalParams,0}, 0, 0x0}, /* (cLog2 [x]) */
    /* 469	*/ {{2,/*0,7        */7168      , cMax        ,PositionalParams,0}, 0, 0x0}, /* (cMax [x y]) */
    /* 470	*/ {{2,/*7,35       */35847     , cMax        ,PositionalParams,0}, 0, 0x0}, /* (cMax [y a]) */
    /* 471	*/ {{2,/*31,30      */30751     , cMax        ,PositionalParams,0}, 0, 0x0}, /* (cMax [z b]) */
    /* 472	*/ {{1,/*0          */0         , cMax        ,AnyParams       ,1}, 0, 0x4}, /* (cMax x <1>) */
    /* 473	*/ {{2,/*0,7        */7168      , cMin        ,PositionalParams,0}, 0, 0x0}, /* (cMin [x y]) */
    /* 474	*/ {{2,/*7,35       */35847     , cMin        ,PositionalParams,0}, 0, 0x0}, /* (cMin [y a]) */
    /* 475	*/ {{2,/*31,30      */30751     , cMin        ,PositionalParams,0}, 0, 0x0}, /* (cMin [z b]) */
    /* 476	*/ {{1,/*0          */0         , cMin        ,AnyParams       ,1}, 0, 0x4}, /* (cMin x <1>) */
    /* 477	*/ {{2,/*15,24      */24591     , cMin        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MIN( % & )@C */
    /* 478	*/ {{1,/*0          */0         , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [x]) */
    /* 479	*/ {{1,/*7          */7         , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [y]) */
    /* 480	*/ {{1,/*92         */92        , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cAdd {x y})]) */
    /* 481	*/ {{1,/*93         */93        , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cAdd {x (cMul {-1 y})})]) */
    /* 482	*/ {{1,/*120        */120       , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cAdd  <1>)]) */
    /* 483	*/ {{1,/*149        */149       , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul {-1 (cAdd  <1>)})]) */
    /* 484	*/ {{1,/*231        */231       , cSin        ,PositionalParams,0}, 0, 0x5}, /* (cSin [(cMul {% x})]) */
    /* 485	*/ {{1,/*246        */246       , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul  <1>)]) */
    /* 486	*/ {{1,/*255        */255       , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul -%@C <1>)]) */
    /* 487	*/ {{1,/*254        */254       , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul %@N <1>)]) */
    /* 488	*/ {{1,/*0          */0         , cSin        ,PositionalParams,0}, 0, 0x4}, /* (cSin [x]) */
    /* 489	*/ {{1,/*273        */273       , cSin        ,PositionalParams,0}, 0, 0x1}, /* (cSin [(cMul % <1>)]) */
    /* 490	*/ {{1,/*214        */214       , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul {x LOG( & )@C})]) */
    /* 491	*/ {{1,/*231        */231       , cSinh       ,PositionalParams,0}, 0, 0x5}, /* (cSinh [(cMul {% x})]) */
    /* 492	*/ {{1,/*246        */246       , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul  <1>)]) */
    /* 493	*/ {{1,/*254        */254       , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul %@N <1>)]) */
    /* 494	*/ {{1,/*255        */255       , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul -%@C <1>)]) */
    /* 495	*/ {{1,/*464        */464       , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cLog [(cPow [& x])])]) */
    /* 496	*/ {{1,/*0          */0         , cSinh       ,PositionalParams,0}, 0, 0x4}, /* (cSinh [x]) */
    /* 497	*/ {{1,/*0          */0         , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [x]) */
    /* 498	*/ {{1,/*15         */15        , cSqrt       ,GroupFunction   ,0}, Constness_Const, 0x0}, /* SQRT( % )@C */
    /* 499	*/ {{1,/*0          */0         , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [x]) */
    /* 500	*/ {{1,/*0          */0         , cTan        ,PositionalParams,0}, 0, 0x4}, /* (cTan [x]) */
    /* 501	*/ {{1,/*116        */116       , cTan        ,PositionalParams,0}, 0, 0x4}, /* (cTan [(cAdd {1.57079632679 (cMul {-1 x})})]) */
    /* 502	*/ {{1,/*117        */117       , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cAdd {1.57079632679 (cMul -1 <1>)})]) */
    /* 503	*/ {{1,/*231        */231       , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul {% x})]) */
    /* 504	*/ {{1,/*246        */246       , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul  <1>)]) */
    /* 505	*/ {{1,/*273        */273       , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul % <1>)]) */
    /* 506	*/ {{1,/*254        */254       , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul %@N <1>)]) */
    /* 507	*/ {{1,/*255        */255       , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul -%@C <1>)]) */
    /* 508	*/ {{1,/*0          */0         , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [x]) */
    /* 509	*/ {{1,/*0          */0         , cTanh       ,PositionalParams,0}, 0, 0x4}, /* (cTanh [x]) */
    /* 510	*/ {{1,/*213        */213       , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul {x LOG( % )@C 0.5})]) */
    /* 511	*/ {{1,/*231        */231       , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul {% x})]) */
    /* 512	*/ {{1,/*246        */246       , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul  <1>)]) */
    /* 513	*/ {{1,/*254        */254       , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul %@N <1>)]) */
    /* 514	*/ {{1,/*255        */255       , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul -%@C <1>)]) */
    /* 515	*/ {{1,/*0          */0         , cTrunc      ,PositionalParams,0}, 0, 0x0}, /* (cTrunc [x]) */
    /* 516	*/ {{2,/*24,15      */15384     , cSub        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* SUB( & % )@C */
    /* 517	*/ {{2,/*24,15      */15384     , cDiv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* DIV( & % )@C */
    /* 518	*/ {{2,/*466,465    */476626    , cDiv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* DIV( LOG( & )@C LOG( % )@C )@C */
    /* 519	*/ {{2,/*57,120     */122937    , cEqual      ,PositionalParams,0}, 0, 0x0}, /* (cEqual [-%@C (cAdd  <1>)]) */
    /* 520	*/ {{2,/*0,7        */7168      , cEqual      ,PositionalParams,0}, 0, 0x12}, /* (cEqual [x y]) */
    /* 521	*/ {{2,/*0,7        */7168      , cEqual      ,PositionalParams,0}, 0, 0x0}, /* (cEqual [x y]) */
    /* 522	*/ {{2,/*0,31       */31744     , cEqual      ,PositionalParams,0}, 0, 0x20}, /* (cEqual [x z]) */
    /* 523	*/ {{2,/*7,31       */31751     , cEqual      ,PositionalParams,0}, 0, 0x24}, /* (cEqual [y z]) */
    /* 524	*/ {{2,/*7,31       */31751     , cEqual      ,PositionalParams,0}, 0, 0x0}, /* (cEqual [y z]) */
    /* 525	*/ {{2,/*57,120     */122937    , cNEqual     ,PositionalParams,0}, 0, 0x0}, /* (cNEqual [-%@C (cAdd  <1>)]) */
    /* 526	*/ {{2,/*0,7        */7168      , cLess       ,PositionalParams,0}, 0, 0x12}, /* (cLess [x y]) */
    /* 527	*/ {{2,/*0,41       */41984     , cLess       ,PositionalParams,0}, 0, 0x4}, /* (cLess [x 0]) */
    /* 528	*/ {{2,/*0,41       */41984     , cLess       ,PositionalParams,0}, 0, 0x0}, /* (cLess [x 0]) */
    /* 529	*/ {{2,/*7,0        */7         , cLess       ,PositionalParams,0}, 0, 0x0}, /* (cLess [y x]) */
    /* 530	*/ {{2,/*0,7        */7168      , cLessOrEq   ,PositionalParams,0}, 0, 0x0}, /* (cLessOrEq [x y]) */
    /* 531	*/ {{2,/*246,289    */296182    , cLessOrEq   ,PositionalParams,0}, 0, 0x0}, /* (cLessOrEq [(cMul  <1>) MUL( 0.5 /%@C )@C]) */
    /* 532	*/ {{2,/*0,7        */7168      , cGreater    ,PositionalParams,0}, 0, 0x12}, /* (cGreater [x y]) */
    /* 533	*/ {{2,/*0,41       */41984     , cGreater    ,PositionalParams,0}, 0, 0x4}, /* (cGreater [x 0]) */
    /* 534	*/ {{2,/*0,41       */41984     , cGreater    ,PositionalParams,0}, 0, 0x0}, /* (cGreater [x 0]) */
    /* 535	*/ {{2,/*7,0        */7         , cGreater    ,PositionalParams,0}, 0, 0x0}, /* (cGreater [y x]) */
    /* 536	*/ {{2,/*0,7        */7168      , cGreaterOrEq,PositionalParams,0}, 0, 0x0}, /* (cGreaterOrEq [x y]) */
    /* 537	*/ {{2,/*246,289    */296182    , cGreaterOrEq,PositionalParams,0}, 0, 0x0}, /* (cGreaterOrEq [(cMul  <1>) MUL( 0.5 /%@C )@C]) */
    /* 538	*/ {{1,/*0          */0         , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [x]) */
    /* 539	*/ {{1,/*245        */245       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cMul {x 0.5 /%@C})]) */
    /* 540	*/ {{1,/*7          */7         , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [y]) */
    /* 541	*/ {{1,/*550        */550       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cAnd {x y})]) */
    /* 542	*/ {{1,/*553        */553       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cAnd {z (cIf [x y (cNot [%])])})]) */
    /* 543	*/ {{1,/*554        */554       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cAnd {z (cIf [x (cNot [%]) y])})]) */
    /* 544	*/ {{1,/*556        */556       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cOr {x y})]) */
    /* 545	*/ {{1,/*31         */31        , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [z]) */
    /* 546	*/ {{1,/*559        */559       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cOr {z (cIf [x y (cNot [%])])})]) */
    /* 547	*/ {{1,/*15         */15        , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [%]) */
    /* 548	*/ {{1,/*560        */560       , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cOr {z (cIf [x (cNot [%]) y])})]) */
    /* 549	*/ {{2,/*538,7      */7706      , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {(cNot [x]) y}) */
    /* 550	*/ {{2,/*0,7        */7168      , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {x y}) */
    /* 551	*/ {{2,/*7,35       */35847     , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {y a}) */
    /* 552	*/ {{2,/*31,30      */30751     , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {z b}) */
    /* 553	*/ {{2,/*31,453     */463903    , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {z (cIf [x y (cNot [%])])}) */
    /* 554	*/ {{2,/*31,456     */466975    , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {z (cIf [x (cNot [%]) y])}) */
    /* 555	*/ {{0,/*           */0         , cAnd        ,AnyParams       ,1}, 0, 0x0}, /* (cAnd  <1>) */
    /* 556	*/ {{2,/*0,7        */7168      , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {x y}) */
    /* 557	*/ {{2,/*538,7      */7706      , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {(cNot [x]) y}) */
    /* 558	*/ {{2,/*7,35       */35847     , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {y a}) */
    /* 559	*/ {{2,/*31,453     */463903    , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {z (cIf [x y (cNot [%])])}) */
    /* 560	*/ {{2,/*31,456     */466975    , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {z (cIf [x (cNot [%]) y])}) */
    /* 561	*/ {{2,/*31,30      */30751     , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {z b}) */
    /* 562	*/ {{1,/*0          */0         , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [x]) */
    /* 563	*/ {{1,/*92         */92        , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cAdd {x y})]) */
    /* 564	*/ {{1,/*131        */131       , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cAdd x <1>)]) */
    /* 565	*/ {{1,/*245        */245       , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cMul {x 0.5 /%@C})]) */
    /* 566	*/ {{1,/*215        */215       , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cMul {x (cAnd  <1>)})]) */
    /* 567	*/ {{1,/*246        */246       , cDeg        ,PositionalParams,0}, 0, 0x0}, /* (cDeg [(cMul  <1>)]) */
    /* 568	*/ {{1,/*246        */246       , cRad        ,PositionalParams,0}, 0, 0x0}, /* (cRad [(cMul  <1>)]) */
    /* 569	*/ {{2,/*0,7        */7168      , cAbsAnd     ,SelectedParams  ,0}, 0, 0x0}, /* (cAbsAnd {x y}) */
    /* 570	*/ {{2,/*0,7        */7168      , cAbsOr      ,SelectedParams  ,0}, 0, 0x0}, /* (cAbsOr {x y}) */
    /* 571	*/ {{1,/*0          */0         , cAbsNot     ,PositionalParams,0}, 0, 0x0}, /* (cAbsNot [x]) */
    /* 572	*/ {{1,/*0          */0         , cAbsNotNot  ,PositionalParams,0}, 0, 0x0}, /* (cAbsNotNot [x]) */
    /* 573	*/ {{3,/*0,7,31     */32513024  , cAbsIf      ,PositionalParams,0}, 0, 0x0}, /* (cAbsIf [x y z]) */
    };

}
namespace FPoptimizer_Grammar
{
    const Rule grammar_rules[262] =
    {
        /* 0:	@R @L (cAbs [x])
         *	->	x
         */		 {ProduceNewTree, 17, 1,/*0          */0         , {1,/*0          */0         , cAbs        ,PositionalParams,0}},
        /* 1:	@R @F (cAtan [(cMul {x (cPow [y %@N])})])
         *	->	(cAtan2 [x (cPow [y -%@C])])
         */		 {ProduceNewTree, 18, 1,/*409        */409       , {1,/*146        */146       , cAtan       ,PositionalParams,0}},
        /* 2:	@R @F (cAtan2 [(cPow [(cAdd {(cMul {(cPow [x 2]) -1}) 1}) 0.5])@D4 x@D4])
         *	->	(cAcos [x])
         */		 {ProduceNewTree, 18, 1,/*403        */403       , {2,/*300,1      */1324      , cAtan2      ,PositionalParams,0}},
        /* 3:	@R @F (cAtan2 [x@D4 (cPow [(cAdd {(cMul {(cPow [x 2]) -1}) 1}) 0.5])@D4])
         *	->	(cAsin [x])
         */		 {ProduceNewTree, 18, 1,/*405        */405       , {2,/*1,300      */307201    , cAtan2      ,PositionalParams,0}},
        /* 4:	@R @F (cAtan2 [(cMul x@P <1>)@D4 (cMul x@P <2>)@D4])
         *	:	(cMul  <1>) (cMul  <2>)
         */		 {ReplaceParams , 18, 2,/*246,247    */253174    , {2,/*248,249    */255224    , cAtan2      ,PositionalParams,0}},
        /* 5:	@R @F (cAtan2 [(cMul x@N@V <1>)@D4 (cMul x@N@V <2>)@D4])
         *	:	(cMul -1 <1>) (cMul -1 <2>)
         */		 {ReplaceParams , 18, 2,/*252,253    */259324    , {2,/*250,251    */257274    , cAtan2      ,PositionalParams,0}},
        /* 6:	@R @F (cCeil [(cMul -1 <1>)])
         *	->	(cMul {-1 (cFloor [(cMul  <1>)])})
         */		 {ProduceNewTree, 18, 1,/*152        */152       , {1,/*252        */252       , cCeil       ,PositionalParams,0}},
        /* 7:	@F (cCos [(cAdd {1.57079632679 (cMul %@N <1>)})])
         *	->	(cSin [(cMul -%@C <1>)])
         */		 {ProduceNewTree, 2, 1,/*486        */486       , {1,/*68         */68        , cCos        ,PositionalParams,0}},
        /* 8:	@F (cCos [(cAdd -1.57079632679 <1>)])
         *	->	(cSin [(cAdd  <1>)])
         */		 {ProduceNewTree, 2, 1,/*482        */482       , {1,/*123        */123       , cCos        ,PositionalParams,0}},
        /* 9:	@F (cCos [(cAdd 1.57079632679 <1>)])
         *	->	(cSin [(cMul {-1 (cAdd  <1>)})])
         */		 {ProduceNewTree, 2, 1,/*483        */483       , {1,/*125        */125       , cCos        ,PositionalParams,0}},
        /* 10:	@F (cCos [(cAdd 3.14159265359 <1>)])
         *	->	(cMul {-1 (cCos [(cAdd  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*151        */151       , {1,/*126        */126       , cCos        ,PositionalParams,0}},
        /* 11:	@F (cCos [(cAdd 0 <1>)])
         *	->	(cCos [(cAdd  <1>)])
         */		 {ProduceNewTree, 2, 1,/*419        */419       , {1,/*124        */124       , cCos        ,PositionalParams,0}},
        /* 12:	@F (cCos [(cAcos [x])])
         *	->	x
         */		 {ProduceNewTree, 2, 1,/*0          */0         , {1,/*403        */403       , cCos        ,PositionalParams,0}},
        /* 13:	@F (cCos [(cMul -1 <1>)])
         *	:	(cMul  <1>)
         */		 {ReplaceParams , 2, 1,/*246        */246       , {1,/*252        */252       , cCos        ,PositionalParams,0}},
        /* 14:	@R @F (cCos [(cAbs [x])])
         *	:	x
         */		 {ReplaceParams , 18, 1,/*0          */0         , {1,/*400        */400       , cCos        ,PositionalParams,0}},
        /* 15:	@F (cCosh [(cAsinh [x])])
         *	->	(cPow [(cAdd {(cPow [x 2]) 1}) 0.5])
         */		 {ProduceNewTree, 2, 1,/*301        */301       , {1,/*406        */406       , cCosh       ,PositionalParams,0}},
        /* 16:	@F (cCosh [(cMul -1 <1>)])
         *	:	(cMul  <1>)
         */		 {ReplaceParams , 2, 1,/*246        */246       , {1,/*252        */252       , cCosh       ,PositionalParams,0}},
        /* 17:	@R @F (cCosh [(cAbs [x])])
         *	:	x
         */		 {ReplaceParams , 18, 1,/*0          */0         , {1,/*400        */400       , cCosh       ,PositionalParams,0}},
        /* 18:	@F (cFloor [(cAdd 0.5 <1>)])
         *	->	(cInt [(cAdd  <1>)])
         */		 {ProduceNewTree, 2, 1,/*458        */458       , {1,/*122        */122       , cFloor      ,PositionalParams,0}},
        /* 19:	@R @F (cFloor [(cMul -1 <1>)])
         *	->	(cMul {-1 (cCeil [(cMul  <1>)])})
         */		 {ProduceNewTree, 18, 1,/*150        */150       , {1,/*252        */252       , cFloor      ,PositionalParams,0}},
        /* 20:	(cIf [x 0 y])
         *	->	(cMul {(cNot [x]) y})
         */		 {ProduceNewTree, 0, 1,/*156        */156       , {3,/*0,41,7     */7382016   , cIf         ,PositionalParams,0}},
        /* 21:	(cIf [x 0 y@L])
         *	->	(cAnd {(cNot [x]) y})
         */		 {ProduceNewTree, 0, 1,/*549        */549       , {3,/*0,41,8     */8430592   , cIf         ,PositionalParams,0}},
        /* 22:	(cIf [x 1 y@L])
         *	->	(cOr {x y})
         */		 {ProduceNewTree, 0, 1,/*556        */556       , {3,/*0,47,8     */8436736   , cIf         ,PositionalParams,0}},
        /* 23:	(cIf [x y 0])
         *	->	(cMul {(cNotNot [x]) y})
         */		 {ProduceNewTree, 0, 1,/*157        */157       , {3,/*0,7,41     */42998784  , cIf         ,PositionalParams,0}},
        /* 24:	(cIf [x y@L 0])
         *	->	(cAnd {x y})
         */		 {ProduceNewTree, 0, 1,/*550        */550       , {3,/*0,8,41     */42999808  , cIf         ,PositionalParams,0}},
        /* 25:	(cIf [x 1 0])
         *	->	(cNotNot [x])
         */		 {ProduceNewTree, 0, 1,/*562        */562       , {3,/*0,47,41    */43039744  , cIf         ,PositionalParams,0}},
        /* 26:	(cIf [x y@L 1])
         *	->	(cOr {(cNot [x]) y})
         */		 {ProduceNewTree, 0, 1,/*557        */557       , {3,/*0,8,47     */49291264  , cIf         ,PositionalParams,0}},
        /* 27:	(cIf [x 0 1])
         *	->	(cNot [x])
         */		 {ProduceNewTree, 0, 1,/*538        */538       , {3,/*0,41,47    */49325056  , cIf         ,PositionalParams,0}},
        /* 28:	(cIf [(cLess [x y])@D12 y@D8 x@D4])
         *	->	(cMax [x y])
         */		 {ProduceNewTree, 0, 1,/*469        */469       , {3,/*526,9,1    */1058318   , cIf         ,PositionalParams,0}},
        /* 29:	(cIf [(cGreater [x y])@D12 y@D8 x@D4])
         *	->	(cMin [x y])
         */		 {ProduceNewTree, 0, 1,/*473        */473       , {3,/*532,9,1    */1058324   , cIf         ,PositionalParams,0}},
        /* 30:	(cIf [(cLess [x y])@D12 x@D4 y@D8])
         *	->	(cMin [x y])
         */		 {ProduceNewTree, 0, 1,/*473        */473       , {3,/*526,1,9    */9438734   , cIf         ,PositionalParams,0}},
        /* 31:	(cIf [(cGreater [x y])@D12 x@D4 y@D8])
         *	->	(cMax [x y])
         */		 {ProduceNewTree, 0, 1,/*469        */469       , {3,/*532,1,9    */9438740   , cIf         ,PositionalParams,0}},
        /* 32:	(cIf [(cLessOrEq [x y]) z a])
         *	:	(cLess [y x]) a z
         */		 {ReplaceParams , 0, 3,/*529,35,31  */32542225  , {3,/*530,31,35  */36732434  , cIf         ,PositionalParams,0}},
        /* 33:	(cIf [(cGreaterOrEq [x y]) z a])
         *	:	(cGreater [y x]) a z
         */		 {ReplaceParams , 0, 3,/*535,35,31  */32542231  , {3,/*536,31,35  */36732440  , cIf         ,PositionalParams,0}},
        /* 34:	@R (cIf [x@P y z])
         *	->	(cAbsIf [x y z])
         */		 {ProduceNewTree, 16, 1,/*573        */573       , {3,/*2,7,31     */32513026  , cIf         ,PositionalParams,0}},
        /* 35:	@R (cIf [(cLess [x 0])@D4 (cCeil [x])@D4 (cFloor [x])@D4])
         *	->	(cTrunc [x])
         */		 {ProduceNewTree, 16, 1,/*515        */515       , {3,/*527,413,434*/455505423 , cIf         ,PositionalParams,0}},
        /* 36:	@R (cIf [(cGreater [x 0])@D4 (cFloor [x])@D4 (cCeil [x])@D4])
         *	->	(cTrunc [x])
         */		 {ProduceNewTree, 16, 1,/*515        */515       , {3,/*533,434,413*/433506837 , cIf         ,PositionalParams,0}},
        /* 37:	@F (cLog [(cMul %@P <1>)])
         *	->	(cAdd {(cLog [(cMul  <1>)]) LOG( % )@C})
         */		 {ProduceNewTree, 2, 1,/*78         */78        , {1,/*256        */256       , cLog        ,PositionalParams,0}},
        /* 38:	@F (cLog [(cMul (cPow [% y]) <1>)])
         *	->	(cAdd {(cMul {LOG( % )@C y}) (cLog [(cMul  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*69         */69        , {1,/*258        */258       , cLog        ,PositionalParams,0}},
        /* 39:	@F (cLog [(cAdd {(cPow [(cAdd {-1 (cPow [x 2])}) 0.5])@D4 x@D4})])
         *	->	(cAcosh [x])
         */		 {ProduceNewTree, 2, 1,/*404        */404       , {1,/*72         */72        , cLog        ,PositionalParams,0}},
        /* 40:	@F (cLog [(cMul {(cAdd {1 x})@D4 (cPow [(cAdd {1 (cMul {-1 x})}) -1])@D4})])
         *	->	(cMul {(cAtanh [x]) 2})
         */		 {ProduceNewTree, 2, 1,/*159        */159       , {1,/*147        */147       , cLog        ,PositionalParams,0}},
        /* 41:	(cMax x@D4 (cMin x <1>)@D4)
         *	:	x
         */		 {ReplaceParams , 0, 1,/*0          */0         , {2,/*1,476      */487425    , cMax        ,AnyParams       ,0}},
        /* 42:	@R (cMax (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cMax [y a]) (cMax [z b])])
         */		 {ReplaceParams , 16, 1,/*445        */445       , {2,/*438,440    */450998    , cMax        ,AnyParams       ,0}},
        /* 43:	(cMin x@D4 (cMax x <1>)@D4)
         *	:	x
         */		 {ReplaceParams , 0, 1,/*0          */0         , {2,/*1,472      */483329    , cMin        ,AnyParams       ,0}},
        /* 44:	@R (cMin (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cMin [y a]) (cMin [z b])])
         */		 {ReplaceParams , 16, 1,/*446        */446       , {2,/*438,440    */450998    , cMin        ,AnyParams       ,0}},
        /* 45:	(cPow [(cMul %@P <1>) &])
         *	->	(cMul {POW( % & )@C (cPow [(cMul  <1>) &])})
         */		 {ProduceNewTree, 0, 1,/*153        */153       , {2,/*256,24     */24832     , cPow        ,PositionalParams,0}},
        /* 46:	(cPow [(cMul %@N <1>) &@E])
         *	->	(cMul {POW( % & )@C (cPow [(cMul  <1>) &])})
         */		 {ProduceNewTree, 0, 1,/*153        */153       , {2,/*254,25     */25854     , cPow        ,PositionalParams,0}},
        /* 47:	(cPow [% (cAdd &@M <1>)])
         *	->	(cMul {POW( % & )@C (cPow [% (cAdd  <1>)])})
         */		 {ProduceNewTree, 0, 1,/*154        */154       , {2,/*15,127     */130063    , cPow        ,PositionalParams,0}},
        /* 48:	(cPow [(cPow [x y@O]) z])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , 0, 2,/*0,162      */165888    , {2,/*311,31     */32055     , cPow        ,PositionalParams,0}},
        /* 49:	(cPow [(cPow [x y@F]) z])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , 0, 2,/*0,162      */165888    , {2,/*312,31     */32056     , cPow        ,PositionalParams,0}},
        /* 50:	(cPow [(cPow [x@P y]) z])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , 0, 2,/*0,162      */165888    , {2,/*313,31     */32057     , cPow        ,PositionalParams,0}},
        /* 51:	(cPow [(cPow [x y])@P z])
         *	:	(cAbs [x]) (cMul {y z})
         */		 {ReplaceParams , 0, 2,/*400,162    */166288    , {2,/*393,31     */32137     , cPow        ,PositionalParams,0}},
        /* 52:	(cPow [(cPow [x y]) z@I])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , 0, 2,/*0,162      */165888    , {2,/*314,32     */33082     , cPow        ,PositionalParams,0}},
        /* 53:	(cPow [(cAbs [x]) y@E])
         *	:	x y
         */		 {ReplaceParams , 0, 2,/*0,7        */7168      , {2,/*400,12     */12688     , cPow        ,PositionalParams,0}},
        /* 54:	(cPow [(cMul (cAbs [x]) <1>) y@E])
         *	:	(cMul x <1>) y
         */		 {ReplaceParams , 0, 2,/*266,7      */7434      , {2,/*265,12     */12553     , cPow        ,PositionalParams,0}},
        /* 55:	@F (cPow [(cAdd [(cPow [x %@E]) (cPow [y &@E])]) 0.5])
         *	->	(cHypot [(cPow [x MUL( % 0.5 )@C]) (cPow [y MUL( & 0.5 )@C])])
         */		 {ProduceNewTree, 2, 1,/*435        */435       , {2,/*66,45      */46146     , cPow        ,PositionalParams,0}},
        /* 56:	@F (cPow [(cAdd {(cPow [x %@E]) (cMul {z@C (cPow [y &@E])})}) 0.5])
         *	->	(cHypot [(cPow [x MUL( % 0.5 )@C]) (cPow [(cMul {y POW( z@C /&@C )@C}) MUL( & 0.5 )@C])])
         */		 {ProduceNewTree, 2, 1,/*436        */436       , {2,/*74,45      */46154     , cPow        ,PositionalParams,0}},
        /* 57:	@F (cPow [(cAdd {(cMul {a@C (cPow [x %@E])}) (cMul {z@C (cPow [y &@E])})}) 0.5])
         *	->	(cHypot [(cPow [(cMul {x POW( a@C /%@C )@C}) MUL( % 0.5 )@C]) (cPow [(cMul {y POW( z@C /&@C )@C}) MUL( & 0.5 )@C])])
         */		 {ProduceNewTree, 2, 1,/*437        */437       , {2,/*70,45      */46150     , cPow        ,PositionalParams,0}},
        /* 58:	@F (cPow [% (cAdd {(cLog [x]) &})])
         *	->	(cMul {POW( % & )@C (cPow [x LOG( % )@C])})
         */		 {ProduceNewTree, 2, 1,/*169        */169       , {2,/*15,82      */83983     , cPow        ,PositionalParams,0}},
        /* 59:	@F (cPow [z@D16 (cAdd (cMul (cPow [(cLog [z]) -1]) <2>) <1>)@D16])
         *	->	(cMul {(cPow [z (cAdd  <1>)]) (cPow [2.71828182846 (cMul  <2>)])})
         */		 {ProduceNewTree, 2, 1,/*168        */168       , {2,/*34,128     */131106    , cPow        ,PositionalParams,0}},
        /* 60:	@F (cPow [z@D16 (cAdd (cPow [(cLog [z]) -1]) <1>)@D16])
         *	->	(cMul {2.71828182846 (cPow [z (cAdd  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*175        */175       , {2,/*34,130     */133154    , cPow        ,PositionalParams,0}},
        /* 61:	@F (cPow [% (cLog [x])])
         *	:	x LOG( % )@C
         */		 {ReplaceParams , 2, 2,/*0,465      */476160    , {2,/*15,460     */471055    , cPow        ,PositionalParams,0}},
        /* 62:	@F (cPow [% (cMul (cLog [x]) <1>)])
         *	:	x (cMul LOG( % )@C <1>)
         */		 {ReplaceParams , 2, 2,/*0,268      */274432    , {2,/*15,267     */273423    , cPow        ,PositionalParams,0}},
        /* 63:	@F (cPow [z@D16 (cMul (cPow [(cLog [z]) -1]) (cLog [x]) <1>)@D16])
         *	:	x (cMul  <1>)
         */		 {ReplaceParams , 2, 2,/*0,246      */251904    , {2,/*34,260     */266274    , cPow        ,PositionalParams,0}},
        /* 64:	@F (cPow [%@D1 (cMul /LOG( % )@C@C (cLog [x]) <1>)@D1])
         *	:	x (cMul  <1>)
         */		 {ReplaceParams , 2, 2,/*0,246      */251904    , {2,/*18,257     */263186    , cPow        ,PositionalParams,0}},
        /* 65:	@F (cSin [(cMul -1 <1>)])
         *	->	(cMul {-1 (cSin [(cMul  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*171        */171       , {1,/*252        */252       , cSin        ,PositionalParams,0}},
        /* 66:	@F (cSin [(cAdd {1.57079632679 (cMul %@N <1>)})])
         *	->	(cCos [(cMul -%@C <1>)])
         */		 {ProduceNewTree, 2, 1,/*421        */421       , {1,/*68         */68        , cSin        ,PositionalParams,0}},
        /* 67:	@F (cSin [(cAdd -1.57079632679 <1>)])
         *	->	(cMul {-1 (cCos [(cAdd  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*151        */151       , {1,/*123        */123       , cSin        ,PositionalParams,0}},
        /* 68:	@F (cSin [(cAdd 1.57079632679 <1>)])
         *	->	(cCos [(cAdd  <1>)])
         */		 {ProduceNewTree, 2, 1,/*419        */419       , {1,/*125        */125       , cSin        ,PositionalParams,0}},
        /* 69:	@F (cSin [(cAdd 3.14159265359 <1>)])
         *	->	(cMul {-1 (cSin [(cAdd  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*170        */170       , {1,/*126        */126       , cSin        ,PositionalParams,0}},
        /* 70:	@F (cSin [(cAdd 0 <1>)])
         *	->	(cSin [(cAdd  <1>)])
         */		 {ProduceNewTree, 2, 1,/*482        */482       , {1,/*124        */124       , cSin        ,PositionalParams,0}},
        /* 71:	@F (cSin [(cAsin [x])])
         *	->	x
         */		 {ProduceNewTree, 2, 1,/*0          */0         , {1,/*405        */405       , cSin        ,PositionalParams,0}},
        /* 72:	@F (cSinh [(cMul -1 <1>)])
         *	->	(cMul {-1 (cSinh [(cMul  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*172        */172       , {1,/*252        */252       , cSinh       ,PositionalParams,0}},
        /* 73:	@F (cSinh [(cAcosh [x])])
         *	->	(cPow [(cAdd {(cPow [x 2]) -1}) 0.5])
         */		 {ProduceNewTree, 2, 1,/*328        */328       , {1,/*404        */404       , cSinh       ,PositionalParams,0}},
        /* 74:	@F (cTan [(cMul -1 <1>)])
         *	->	(cMul {-1 (cTan [(cMul  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*173        */173       , {1,/*252        */252       , cTan        ,PositionalParams,0}},
        /* 75:	@F (cTan [(cAtan [x])])
         *	->	x
         */		 {ProduceNewTree, 2, 1,/*0          */0         , {1,/*408        */408       , cTan        ,PositionalParams,0}},
        /* 76:	@F (cTan [(cAtan2 [x y])])
         *	->	(cMul {x (cPow [y -1])})
         */		 {ProduceNewTree, 2, 1,/*176        */176       , {1,/*410        */410       , cTan        ,PositionalParams,0}},
        /* 77:	@F (cTanh [(cMul -1 <1>)])
         *	->	(cMul {-1 (cTanh [(cMul  <1>)])})
         */		 {ProduceNewTree, 2, 1,/*177        */177       , {1,/*252        */252       , cTanh       ,PositionalParams,0}},
        /* 78:	(cAdd % (cIf [x & z@C]))
         *	:	(cIf [x ADD( % & )@C ADD( % z@C )@C])
         */		 {ReplaceParams , 0, 1,/*442        */442       , {2,/*15,439     */449551    , cAdd        ,AnyParams       ,0}},
        /* 79:	(cAdd (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cAdd {y a}) (cAdd {z b})])
         */		 {ReplaceParams , 0, 1,/*441        */441       , {2,/*438,440    */450998    , cAdd        ,AnyParams       ,0}},
        /* 80:	(cAdd (cMul (cPow [x %@I@P]) <1>)@D4 (cMul (cPow [x &@I]) <2>)@D4)
         *	:	(cMul {(cPow [x MIN( % & )@C]) (cAdd {(cMul (cPow [x (cAdd {% -MIN( % & )@C@C})]) <1>) (cMul (cPow [x (cAdd {& -MIN( % & )@C@C})]) <2>)})})
         */		 {ReplaceParams , 0, 1,/*167        */167       , {2,/*261,262    */268549    , cAdd        ,AnyParams       ,0}},
        /* 81:	(cAdd (cMul %@P <1>)@D1 (cMul -%@C <2>)@D1)
         *	:	(cMul {% (cAdd {(cMul  <1>) (cMul -1 <2>)})})
         */		 {ReplaceParams , 0, 1,/*180        */180       , {2,/*269,270    */276749    , cAdd        ,AnyParams       ,0}},
        /* 82:	(cAdd %@M@D1 (cMul -%@C <2>)@D1)
         *	:	(cMul {% (cAdd {1 (cMul -1 <2>)})})
         */		 {ReplaceParams , 0, 1,/*181        */181       , {2,/*20,270     */276500    , cAdd        ,AnyParams       ,0}},
        /* 83:	(cAdd (cMul {%@P (cPow [x 2])})@D4 (cMul {& x y})@D4)
         *	:	(cPow [(cAdd {(cMul {x SQRT( % )@C}) (cMul {y MUL( 0.5 MUL( & /SQRT( % )@C@C )@C )@C})}) 2]) (cMul {(cPow [y 2]) MUL( MUL( -0.25 /%@C )@C POW( & 2 )@C )@C})
         */		 {ReplaceParams , 0, 2,/*306,186    */190770    , {2,/*182,185    */189622    , cAdd        ,AnyParams       ,0}},
        /* 84:	(cAdd (cMul {%@P (cAdd {(cPow [x 2]) (cPow [y 2])})})@D12 (cMul {&@P x y})@D12)
         *	:	(cMul {% (cPow [(cAdd {x y}) 2])}) (cMul {ADD( & MUL( -2 % )@C )@C x y})
         */		 {ReplaceParams , 0, 2,/*188,190    */194748    , {2,/*187,189    */193723    , cAdd        ,AnyParams       ,0}},
        /* 85:	(cAdd (cMul {%@P (cAdd {(cPow [x 2]) (cPow [y 2])})})@D12 (cMul {&@N x y})@D12)
         *	:	(cMul {% (cPow [(cAdd {x (cMul {-1 y})}) 2])}) (cMul {ADD( & MUL( 2 % )@C )@C x y})
         */		 {ReplaceParams , 0, 2,/*191,198    */202943    , {2,/*187,192    */196795    , cAdd        ,AnyParams       ,0}},
        /* 86:	(cAdd (cMul {% (cPow [x y])})@D12 (cPow [x (cMul {2 y})])@D12)
         *	:	(cPow [(cAdd {(cPow [x y]) MUL( % 0.5 )@C}) 2]) -POW( MUL( % 0.5 )@C 2 )@C@C
         */		 {ReplaceParams , 0, 2,/*307,58     */59699     , {2,/*164,291    */298148    , cAdd        ,AnyParams       ,0}},
        /* 87:	(cAdd (cMul {% (cPow [x (cMul {& y})])})@D14 (cPow [x (cMul {MUL( & 2 )@C y})])@D14)
         *	:	(cPow [(cAdd {(cPow [x (cMul {y &})]) MUL( % 0.5 )@C}) 2]) -POW( MUL( % 0.5 )@C 2 )@C@C
         */		 {ReplaceParams , 0, 2,/*322,58     */59714     , {2,/*183,318    */325815    , cAdd        ,AnyParams       ,0}},
        /* 88:	(cAdd (cMul {% (cPow [& y])})@D10 (cPow [POW( & 2 )@C y])@D10)
         *	:	(cPow [(cAdd {(cPow [& y]) MUL( % 0.5 )@C}) 2]) -POW( MUL( % 0.5 )@C 2 )@C@C
         */		 {ReplaceParams , 0, 2,/*332,58     */59724     , {2,/*184,335    */343224    , cAdd        ,AnyParams       ,0}},
        /* 89:	@F (cAdd (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) 0.5]) <1>)
         *	->	(cPow [2.71828182846 (cAsinh [(cAdd  <1>)])])
         */		 {ProduceNewTree, 2, 1,/*337        */337       , {1,/*333        */333       , cAdd        ,AnyParams       ,1}},
        /* 90:	@F (cAdd (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) -0.5]) <1>)
         *	->	(cPow [0.367879441171 (cAsinh [(cAdd  <1>)])])
         */		 {ProduceNewTree, 2, 1,/*336        */336       , {1,/*338        */338       , cAdd        ,AnyParams       ,1}},
        /* 91:	@F (cAdd (cPow [(cAdd {1 (cPow [x 2])}) 0.5])@D4 x@D4)
         *	:	(cPow [2.71828182846 (cAsinh [x])])
         */		 {ReplaceParams , 2, 1,/*340        */340       , {2,/*339,1      */1363      , cAdd        ,AnyParams       ,0}},
        /* 92:	@F (cAdd (cPow [(cAdd {1 (cPow [x 2])}) -0.5])@D4 x@D4)
         *	:	(cPow [0.367879441171 (cAsinh [x])])
         */		 {ReplaceParams , 2, 1,/*342        */342       , {2,/*341,1      */1365      , cAdd        ,AnyParams       ,0}},
        /* 93:	@F (cAdd (cLog [x]) (cLog [y]))
         *	:	(cLog [(cMul {x y})])
         */		 {ReplaceParams , 2, 1,/*463        */463       , {2,/*460,461    */472524    , cAdd        ,AnyParams       ,0}},
        /* 94:	@F (cAdd (cPow [(cSin [x]) 2])@D4 (cPow [(cCos [x]) 2])@D4)
         *	:	1
         */		 {ReplaceParams , 2, 1,/*47         */47        , {2,/*359,348    */356711    , cAdd        ,AnyParams       ,0}},
        /* 95:	@F (cAdd 1 (cMul {-1 (cPow [(cSin [x]) 2])}))
         *	:	(cPow [(cCos [x]) 2])
         */		 {ReplaceParams , 2, 1,/*349        */349       , {2,/*47,196     */200751    , cAdd        ,AnyParams       ,0}},
        /* 96:	@F (cAdd 1 (cMul {-1 (cPow [(cCos [x]) 2])}))
         *	:	(cPow [(cSin [x]) 2])
         */		 {ReplaceParams , 2, 1,/*360        */360       , {2,/*47,195     */199727    , cAdd        ,AnyParams       ,0}},
        /* 97:	@F (cAdd (cMul {(cSin [x]) (cCos [y])})@D12 (cMul {(cCos [x]) (cSin [y])})@D12)
         *	:	(cSin [(cAdd {x y})])
         */		 {ReplaceParams , 2, 1,/*480        */480       , {2,/*205,202    */207053    , cAdd        ,AnyParams       ,0}},
        /* 98:	@F (cAdd (cMul {(cSin [x]) (cCos [y])})@D12 (cMul {(cCos [x]) (cSin [y]) -1})@D12)
         *	:	(cSin [(cAdd {x (cMul {-1 y})})])
         */		 {ReplaceParams , 2, 1,/*481        */481       , {2,/*205,203    */208077    , cAdd        ,AnyParams       ,0}},
        /* 99:	@F (cAdd (cMul {(cCos [x]) (cCos [y])})@D12 (cMul {(cSin [x]) (cSin [y])})@D12)
         *	:	(cCos [(cAdd {x y})])
         */		 {ReplaceParams , 2, 1,/*417        */417       , {2,/*200,206    */211144    , cAdd        ,AnyParams       ,0}},
        /* 100:	@F (cAdd (cMul {(cCos [x]) (cCos [y]) -1})@D12 (cMul {(cSin [x]) (cSin [y])})@D12)
         *	:	(cMul {-1 (cCos [(cAdd {x (cMul {-1 y})})])})
         */		 {ReplaceParams , 2, 1,/*209        */209       , {2,/*201,206    */211145    , cAdd        ,AnyParams       ,0}},
        /* 101:	@F (cAdd (cMul {(cCos [x]) (cCos [y])})@D12 (cMul {(cSin [x]) (cSin [y]) -1})@D12)
         *	:	(cCos [(cAdd {x (cMul {-1 y})})])
         */		 {ReplaceParams , 2, 1,/*418        */418       , {2,/*200,210    */215240    , cAdd        ,AnyParams       ,0}},
        /* 102:	@F (cAdd (cPow [& x])@D6 (cMul {-1 (cPow [/&@C x])})@D6)
         *	:	(cMul {(cSinh [(cLog [(cPow [& x])])]) 2})
         */		 {ReplaceParams , 2, 1,/*212        */212       , {2,/*361,207    */212329    , cAdd        ,AnyParams       ,0}},
        /* 103:	@F (cAdd (cPow [& x])@D6 (cPow [/&@C x])@D6)
         *	:	(cMul {(cCosh [(cLog [(cPow [& x])])]) 2})
         */		 {ReplaceParams , 2, 1,/*204        */204       , {2,/*361,364    */373097    , cAdd        ,AnyParams       ,0}},
        /* 104:	@F (cAdd (cMul {-1 (cPow [& x])})@D6 (cPow [/&@C x])@D6)
         *	:	(cMul {(cSinh [(cMul {x LOG( & )@C})]) -2})
         */		 {ReplaceParams , 2, 1,/*211        */211       , {2,/*208,364    */372944    , cAdd        ,AnyParams       ,0}},
        /* 105:	@F (cAdd (cMul {% (cPow [& x])})@D7 (cMul {-%@C (cPow [/&@C x])})@D7)
         *	:	(cMul {(cSinh [(cMul {x LOG( & )@C})]) 2 %})
         */		 {ReplaceParams , 2, 1,/*217        */217       , {2,/*216,197    */201944    , cAdd        ,AnyParams       ,0}},
        /* 106:	@F (cAdd (cMul {% (cPow [& x])})@D7 (cMul {% (cPow [/&@C x])})@D7)
         *	:	(cMul {(cCosh [(cMul {x LOG( & )@C})]) 2 %})
         */		 {ReplaceParams , 2, 1,/*221        */221       , {2,/*216,218    */223448    , cAdd        ,AnyParams       ,0}},
        /* 107:	@F (cAdd (cCosh [x])@D4 (cSinh [x])@D4)
         *	:	(cPow [2.71828182846 x])
         */		 {ReplaceParams , 2, 1,/*367        */367       , {2,/*425,496    */508329    , cAdd        ,AnyParams       ,0}},
        /* 108:	@F (cAdd (cMul {(cCosh [x]) -1})@D4 (cSinh [x])@D4)
         *	:	(cMul {(cPow [0.367879441171 x]) -1})
         */		 {ReplaceParams , 2, 1,/*219        */219       , {2,/*222,496    */508126    , cAdd        ,AnyParams       ,0}},
        /* 109:	@F (cAdd (cCosh [x])@D4 (cMul {(cPow [2.71828182846 x]) -1})@D4)
         *	:	(cMul {-1 (cSinh [x])})
         */		 {ReplaceParams , 2, 1,/*224        */224       , {2,/*425,220    */225705    , cAdd        ,AnyParams       ,0}},
        /* 110:	@F (cAdd (cSinh [x])@D4 (cMul {(cPow [2.71828182846 x]) -1})@D4)
         *	:	(cMul {-1 (cCosh [x])})
         */		 {ReplaceParams , 2, 1,/*223        */223       , {2,/*496,220    */225776    , cAdd        ,AnyParams       ,0}},
        /* 111:	@F (cAdd (cCosh [x])@D4 (cMul {(cSinh [x]) -1})@D4)
         *	:	(cPow [0.367879441171 x])
         */		 {ReplaceParams , 2, 1,/*365        */365       , {2,/*425,225    */230825    , cAdd        ,AnyParams       ,0}},
        /* 112:	@F (cAdd (cMul {(cSinh [x]) -1})@D4 (cPow [2.71828182846 x])@D4)
         *	:	(cCosh [x])
         */		 {ReplaceParams , 2, 1,/*426        */426       , {2,/*225,368    */377057    , cAdd        ,AnyParams       ,0}},
        /* 113:	@F (cAdd (cMul {(cCosh [x]) -1})@D4 (cPow [2.71828182846 x])@D4)
         *	:	(cSinh [x])
         */		 {ReplaceParams , 2, 1,/*497        */497       , {2,/*222,368    */377054    , cAdd        ,AnyParams       ,0}},
        /* 114:	@F (cAdd (cCosh [x])@D4 (cMul {(cPow [0.367879441171 x]) -1})@D4)
         *	:	(cSinh [x])
         */		 {ReplaceParams , 2, 1,/*497        */497       , {2,/*425,199    */204201    , cAdd        ,AnyParams       ,0}},
        /* 115:	@F (cAdd (cSinh [x])@D4 (cPow [0.367879441171 x])@D4)
         *	:	(cCosh [x])
         */		 {ReplaceParams , 2, 1,/*426        */426       , {2,/*496,366    */375280    , cAdd        ,AnyParams       ,0}},
        /* 116:	@F (cAdd (cMul {(cCosh [x]) -1})@D4 (cPow [0.367879441171 x])@D4)
         *	:	(cMul {-1 (cSinh [x])})
         */		 {ReplaceParams , 2, 1,/*224        */224       , {2,/*222,366    */375006    , cAdd        ,AnyParams       ,0}},
        /* 117:	@F (cAdd (cMul {-1 (cPow [x 2])})@D4 (cMul {% x})@D4)
         *	:	(cMul {-1 (cPow [(cAdd {x MUL( % -0.5 )@C}) 2])}) POW( MUL( % 0.5 )@C 2 )@C
         */		 {ReplaceParams , 2, 2,/*229,398    */407781    , {2,/*226,228    */233698    , cAdd        ,AnyParams       ,0}},
        /* 118:	@F (cAdd (cPow [x 2])@D4 (cMul {% x})@D4)
         *	:	(cPow [(cAdd {x MUL( % 0.5 )@C}) 2]) -POW( MUL( % 0.5 )@C 2 )@C@C
         */		 {ReplaceParams , 2, 2,/*371,58     */59763     , {2,/*370,228    */233842    , cAdd        ,AnyParams       ,0}},
        /* 119:	(cMul (cPow [(cMul x <2>) -1])@D4 x@D4)
         *	:	(cPow [(cMul  <2>) -1])
         */		 {ReplaceParams , 0, 1,/*372        */372       , {2,/*373,1      */1397      , cMul        ,AnyParams       ,0}},
        /* 120:	(cMul (cAdd (cMul %@M <1>) <2>) &)
         *	:	(cAdd {(cMul % & <1>) (cMul {& (cAdd  <2>)})})
         */		 {ReplaceParams , 0, 1,/*96         */96        , {2,/*129,24     */24705     , cMul        ,AnyParams       ,0}},
        /* 121:	(cMul (cAdd %@M <1>) &)
         *	:	(cAdd {MUL( % & )@C (cMul {& (cAdd  <1>)})})
         */		 {ReplaceParams , 0, 1,/*97         */97        , {2,/*132,24     */24708     , cMul        ,AnyParams       ,0}},
        /* 122:	(cMul % (cIf [x & z@C]))
         *	:	(cIf [x MUL( % & )@C MUL( % z@C )@C])
         */		 {ReplaceParams , 0, 1,/*444        */444       , {2,/*15,439     */449551    , cMul        ,AnyParams       ,0}},
        /* 123:	(cMul (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cMul {y a}) (cMul {z b})])
         */		 {ReplaceParams , 0, 1,/*443        */443       , {2,/*438,440    */450998    , cMul        ,AnyParams       ,0}},
        /* 124:	(cMul (cPow [x y])@D4 (cAdd {%@1 (cPow [x z])})@D4)
         *	:	(cAdd {(cMul {(cPow [x y]) %}) (cPow [x (cAdd {y z})])})
         */		 {ReplaceParams , 0, 1,/*101        */101       , {2,/*374,100    */102774    , cMul        ,AnyParams       ,0}},
        /* 125:	(cMul (cPow [& y]) (cAdd {1 (cPow [x@P z])}))
         *	:	(cAdd {(cPow [& y]) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})])})
         */		 {ReplaceParams , 0, 1,/*109        */109       , {2,/*325,105    */107845    , cMul        ,AnyParams       ,0}},
        /* 126:	(cMul (cPow [& y]) (cAdd {-1 (cPow [x@P z])}))
         *	:	(cAdd {(cMul {(cPow [& y]) -1}) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})])})
         */		 {ReplaceParams , 0, 1,/*106        */106       , {2,/*325,102    */104773    , cMul        ,AnyParams       ,0}},
        /* 127:	(cMul (cPow [& y]) (cAdd {(cMul {% (cPow [x@P z])})@D1 %@D1}))
         *	:	% (cAdd {(cPow [& y]) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})])})
         */		 {ReplaceParams , 0, 2,/*15,109     */111631    , {2,/*325,107    */109893    , cMul        ,AnyParams       ,0}},
        /* 128:	(cMul (cPow [& y]) (cAdd {(cMul {% (cPow [x@P z])})@D1 -%@C@D1}))
         *	:	% (cAdd {(cMul {(cPow [& y]) -1}) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )@C@C})})])})
         */		 {ReplaceParams , 0, 2,/*15,106     */108559    , {2,/*325,108    */110917    , cMul        ,AnyParams       ,0}},
        /* 129:	@F (cMul {%@D1 (cAdd {1 (cMul {(cLog [x]) /%@C})})@D1})
         *	->	(cAdd {(cLog [x]) %})
         */		 {ProduceNewTree, 2, 1,/*113        */113       , {2,/*18,110     */112658    , cMul        ,SelectedParams  ,0}},
        /* 130:	@F (cMul 57.2957795131 <1>)
         *	->	(cDeg [(cMul  <1>)])
         */		 {ProduceNewTree, 2, 1,/*567        */567       , {1,/*52         */52        , cMul        ,AnyParams       ,1}},
        /* 131:	@F (cMul 0.0174532925199 <1>)
         *	->	(cRad [(cMul  <1>)])
         */		 {ProduceNewTree, 2, 1,/*568        */568       , {1,/*42         */42        , cMul        ,AnyParams       ,1}},
        /* 132:	@F (cMul (cLog [x]) 0.434294481903)
         *	:	(cLog10 [x])
         */		 {ReplaceParams , 2, 1,/*467        */467       , {2,/*460,44     */45516     , cMul        ,AnyParams       ,0}},
        /* 133:	@F (cMul (cPow [(cLog [x]) -1]) 2.30258509299)
         *	:	(cPow [(cLog10 [x]) -1])
         */		 {ReplaceParams , 2, 1,/*356        */356       , {2,/*355,50     */51555     , cMul        ,AnyParams       ,0}},
        /* 134:	@F (cMul (cLog [x]) 1.44269504089)
         *	:	(cLog2 [x])
         */		 {ReplaceParams , 2, 1,/*468        */468       , {2,/*460,48     */49612     , cMul        ,AnyParams       ,0}},
        /* 135:	@F (cMul (cPow [(cLog [x]) -1]) 0.69314718056)
         *	:	(cPow [(cLog2 [x]) -1])
         */		 {ReplaceParams , 2, 1,/*357        */357       , {2,/*355,46     */47459     , cMul        ,AnyParams       ,0}},
        /* 136:	@F (cMul (cExp [x]) (cExp [y]))
         *	:	(cExp [(cAdd {x y})])
         */		 {ReplaceParams , 2, 1,/*429        */429       , {2,/*427,428    */438699    , cMul        ,AnyParams       ,0}},
        /* 137:	@F (cMul (cExp2 [x]) (cExp2 [y]))
         *	:	(cExp2 [(cAdd {x y})])
         */		 {ReplaceParams , 2, 1,/*432        */432       , {2,/*430,431    */441774    , cMul        ,AnyParams       ,0}},
        /* 138:	@F (cMul -1 (cSin [(cMul %@N <1>)]))
         *	:	(cSin [(cMul -%@C <1>)])
         */		 {ReplaceParams , 2, 1,/*486        */486       , {2,/*38,487     */498726    , cMul        ,AnyParams       ,0}},
        /* 139:	@F (cMul -1 (cSinh [(cMul %@N <1>)]))
         *	:	(cSinh [(cMul -%@C <1>)])
         */		 {ReplaceParams , 2, 1,/*494        */494       , {2,/*38,493     */504870    , cMul        ,AnyParams       ,0}},
        /* 140:	@F (cMul (cPow [(cSinh [x]) -1])@D4 (cCosh [x])@D4)
         *	:	(cPow [(cTanh [x]) -1])
         */		 {ReplaceParams , 2, 1,/*382        */382       , {2,/*379,425    */435579    , cMul        ,AnyParams       ,0}},
        /* 141:	@F (cMul (cTanh [x])@D4 (cCosh [x])@D4)
         *	:	(cSinh [x])
         */		 {ReplaceParams , 2, 1,/*497        */497       , {2,/*509,425    */435709    , cMul        ,AnyParams       ,0}},
        /* 142:	@F (cMul (cPow [(cTanh [x]) -1])@D4 (cSinh [x])@D4)
         *	:	(cCosh [x])
         */		 {ReplaceParams , 2, 1,/*426        */426       , {2,/*383,496    */508287    , cMul        ,AnyParams       ,0}},
        /* 143:	@F (cMul (cPow [(cTan [x]) -1])@D4 (cSin [x])@D4)
         *	:	(cCos [x])
         */		 {ReplaceParams , 2, 1,/*414        */414       , {2,/*380,488    */500092    , cMul        ,AnyParams       ,0}},
        /* 144:	@F (cMul (cSin [x])@D4 (cPow [(cCos [x]) -1])@D4)
         *	:	(cTan [x])
         */		 {ReplaceParams , 2, 1,/*499        */499       , {2,/*488,344    */352744    , cMul        ,AnyParams       ,0}},
        /* 145:	@F (cMul (cTan [x])@D4 (cPow [(cSin [x]) -1])@D4)
         *	:	(cPow [(cCos [x]) -1])
         */		 {ReplaceParams , 2, 1,/*345        */345       , {2,/*500,358    */367092    , cMul        ,AnyParams       ,0}},
        /* 146:	@F (cMul (cPow [(cSin [x]) -1])@D4 (cCos [x])@D4)
         *	:	(cPow [(cTan [x]) -1])
         */		 {ReplaceParams , 2, 1,/*381        */381       , {2,/*358,415    */425318    , cMul        ,AnyParams       ,0}},
        /* 147:	@F (cMul (cTan [x])@D4 (cCos [x])@D4)
         *	:	(cSin [x])
         */		 {ReplaceParams , 2, 1,/*478        */478       , {2,/*500,415    */425460    , cMul        ,AnyParams       ,0}},
        /* 148:	@F (cMul (cTan [(cAdd {1.57079632679 (cMul {-1 x})})])@D4 (cTan [x])@D4)
         *	:	1
         */		 {ReplaceParams , 2, 1,/*47         */47        , {2,/*501,500    */512501    , cMul        ,AnyParams       ,0}},
        /* 149:	@F (cMul (cSin [(cMul % <1>)])@D1 (cPow [(cCos [(cMul -%@C <1>)]) -1])@D1)
         *	:	(cTan [(cMul % <1>)])
         */		 {ReplaceParams , 2, 1,/*505        */505       , {2,/*489,347    */355817    , cMul        ,AnyParams       ,0}},
        /* 150:	@F (cMul (cTan [(cAdd {1.57079632679 (cMul -1 <1>)})]) (cTan [(cMul  <1>)]))
         *	:	1
         */		 {ReplaceParams , 2, 1,/*47         */47        , {2,/*502,504    */516598    , cMul        ,AnyParams       ,0}},
        /* 151:	@F (cMul -1 (cTan [(cMul %@N <1>)]))
         *	:	(cTan [(cMul -%@C <1>)])
         */		 {ReplaceParams , 2, 1,/*507        */507       , {2,/*38,506     */518182    , cMul        ,AnyParams       ,0}},
        /* 152:	@F (cMul (cSinh [x])@D4 (cPow [(cCosh [x]) -1])@D4)
         *	:	(cTanh [x])
         */		 {ReplaceParams , 2, 1,/*508        */508       , {2,/*496,350    */358896    , cMul        ,AnyParams       ,0}},
        /* 153:	@F (cMul (cTanh [x])@D4 (cPow [(cSinh [x]) -1])@D4)
         *	:	(cPow [(cCosh [x]) -1])
         */		 {ReplaceParams , 2, 1,/*351        */351       , {2,/*509,379    */388605    , cMul        ,AnyParams       ,0}},
        /* 154:	@F (cMul (cSinh [(cMul {% x})])@D5 (cPow [(cCosh [(cMul {-%@C x})]) -1])@D5)
         *	:	(cTanh [(cMul {% x})])
         */		 {ReplaceParams , 2, 1,/*511        */511       , {2,/*491,352    */360939    , cMul        ,AnyParams       ,0}},
        /* 155:	@F (cMul (cSin [(cMul {% x})])@D5 (cPow [(cCos [(cMul {-%@C x})]) -1])@D5)
         *	:	(cTan [(cMul {% x})])
         */		 {ReplaceParams , 2, 1,/*503        */503       , {2,/*484,346    */354788    , cMul        ,AnyParams       ,0}},
        /* 156:	@F (cMul -1 (cTanh [(cMul %@N <1>)]))
         *	:	(cTanh [(cMul -%@C <1>)])
         */		 {ReplaceParams , 2, 1,/*514        */514       , {2,/*38,513     */525350    , cMul        ,AnyParams       ,0}},
        /* 157:	@F (cMul (cAdd {-1 (cPow [% x])})@D5 (cPow [(cAdd {1 (cPow [% x])}) -1])@D5)
         *	:	(cTanh [(cMul {x LOG( % )@C 0.5})])
         */		 {ReplaceParams , 2, 1,/*510        */510       , {2,/*103,385    */394343    , cMul        ,AnyParams       ,0}},
        /* 158:	@F (cMul (cAdd {1 (cPow [% x])})@D5 (cPow [(cAdd {-1 (cPow [% x])}) -1])@D5)
         *	:	(cPow [(cTanh [(cMul {x LOG( % )@C 0.5})]) -1])
         */		 {ReplaceParams , 2, 1,/*386        */386       , {2,/*115,343    */351347    , cMul        ,AnyParams       ,0}},
        /* 159:	@F (cMul (cSinh [x])@D4 (cPow [(cCosh [x]) %])@D4)
         *	:	(cTanh [x]) (cPow [(cCosh [x]) ADD( % 1 )@C])
         */		 {ReplaceParams , 2, 2,/*508,354    */363004    , {2,/*496,353    */361968    , cMul        ,AnyParams       ,0}},
        /* 160:	@R (cMul (cAdd (cIf [(cLess [x 0]) %@D1 -%@C@D1]) <1>)@D4 x@D4)
         *	:	(cAdd {(cMul {(cAbs [x]) -%@C}) (cMul {x (cAdd  <1>)})})
         */		 {ReplaceParams , 16, 1,/*118        */118       , {2,/*133,1      */1157      , cMul        ,AnyParams       ,0}},
        /* 161:	@R (cMul (cAdd (cIf [(cGreater [x 0]) %@D1 -%@C@D1]) <1>)@D4 x@D4)
         *	:	(cAdd {(cMul {(cAbs [x]) %}) (cMul {x (cAdd  <1>)})})
         */		 {ReplaceParams , 16, 1,/*119        */119       , {2,/*134,1      */1158      , cMul        ,AnyParams       ,0}},
        /* 162:	@R (cMul (cAbs [x]) (cAbs [y]))
         *	:	(cAbs [(cMul {x y})])
         */		 {ReplaceParams , 16, 1,/*402        */402       , {2,/*400,401    */411024    , cMul        ,AnyParams       ,0}},
        /* 163:	@R (cMul (cIf [(cLess [x 0]) %@D1 -%@C@D1])@D4 x@D4)
         *	:	(cAbs [x]) -%@C
         */		 {ReplaceParams , 16, 2,/*400,57     */58768     , {2,/*448,1      */1472      , cMul        ,AnyParams       ,0}},
        /* 164:	@R (cMul (cIf [(cGreater [x 0]) %@D1 -%@C@D1])@D4 x@D4)
         *	:	(cAbs [x]) %
         */		 {ReplaceParams , 16, 2,/*400,15     */15760     , {2,/*450,1      */1474      , cMul        ,AnyParams       ,0}},
        /* 165:	@R @L (cMul (cAbs [x]))
         *	:	x
         */		 {ReplaceParams , 17, 1,/*0          */0         , {1,/*400        */400       , cMul        ,AnyParams       ,0}},
        /* 166:	@R @L (cMul %@N)
         *	:	-%@C
         */		 {ReplaceParams , 17, 1,/*57         */57        , {1,/*14         */14        , cMul        ,AnyParams       ,0}},
        /* 167:	@I (cEqual [0 x])
         *	->	(cNot [x])
         */		 {ProduceNewTree, 4, 1,/*538        */538       , {2,/*41,0       */41        , cEqual      ,PositionalParams,0}},
        /* 168:	@I (cEqual [1 x@L])
         *	->	x
         */		 {ProduceNewTree, 4, 1,/*0          */0         , {2,/*47,5       */5167      , cEqual      ,PositionalParams,0}},
        /* 169:	@R (cEqual [0 (cAbs [x])])
         *	:	x 0
         */		 {ReplaceParams , 16, 2,/*0,41       */41984     , {2,/*41,400     */409641    , cEqual      ,PositionalParams,0}},
        /* 170:	@R (cEqual [(cAdd % <1>) &])
         *	:	(cAdd  <1>) SUB( & % )@C
         */		 {ReplaceParams , 16, 2,/*120,516    */528504    , {2,/*137,24     */24713     , cEqual      ,PositionalParams,0}},
        /* 171:	@R (cEqual [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -%@C <2>)
         */		 {ReplaceParams , 16, 2,/*120,139    */142456    , {2,/*137,138    */141449    , cEqual      ,PositionalParams,0}},
        /* 172:	@R (cEqual [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , 16, 2,/*120,121    */124024    , {2,/*135,136    */139399    , cEqual      ,PositionalParams,0}},
        /* 173:	@R @F (cEqual [(cMul % <1>) &])
         *	:	(cMul  <1>) DIV( & % )@C
         */		 {ReplaceParams , 18, 2,/*246,517    */529654    , {2,/*273,24     */24849     , cEqual      ,PositionalParams,0}},
        /* 174:	@R @F (cEqual [(cPow [x %@P]) &])
         *	:	(cPow [(cPow [x %]) /%@C]) POW( & /%@C )@C
         */		 {ReplaceParams , 18, 2,/*388,399    */408964    , {2,/*387,24     */24963     , cEqual      ,PositionalParams,0}},
        /* 175:	@R @F (cEqual [(cMul % <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul DIV( & % )@C <2>)
         */		 {ReplaceParams , 18, 2,/*246,276    */282870    , {2,/*273,275    */281873    , cEqual      ,PositionalParams,0}},
        /* 176:	@R @F (cEqual [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , 18, 2,/*0,7        */7168      , {2,/*390,392    */401798    , cEqual      ,PositionalParams,0}},
        /* 177:	@R @F (cEqual [&@P (cPow [%@P x])])
         *	:	DIV( LOG( & )@C LOG( % )@C )@C x
         */		 {ReplaceParams , 18, 2,/*518,0      */518       , {2,/*28,391     */400412    , cEqual      ,PositionalParams,0}},
        /* 178:	@I (cNEqual [0 x])
         *	->	(cNotNot [x])
         */		 {ProduceNewTree, 4, 1,/*562        */562       , {2,/*41,0       */41        , cNEqual     ,PositionalParams,0}},
        /* 179:	@I (cNEqual [1 x@L])
         *	->	(cNot [x])
         */		 {ProduceNewTree, 4, 1,/*538        */538       , {2,/*47,5       */5167      , cNEqual     ,PositionalParams,0}},
        /* 180:	@R (cNEqual [0 (cAbs [x])])
         *	:	x 0
         */		 {ReplaceParams , 16, 2,/*0,41       */41984     , {2,/*41,400     */409641    , cNEqual     ,PositionalParams,0}},
        /* 181:	@R (cNEqual [(cAdd % <1>) &])
         *	:	(cAdd  <1>) SUB( & % )@C
         */		 {ReplaceParams , 16, 2,/*120,516    */528504    , {2,/*137,24     */24713     , cNEqual     ,PositionalParams,0}},
        /* 182:	@R (cNEqual [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -%@C <2>)
         */		 {ReplaceParams , 16, 2,/*120,139    */142456    , {2,/*137,138    */141449    , cNEqual     ,PositionalParams,0}},
        /* 183:	@R (cNEqual [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , 16, 2,/*120,121    */124024    , {2,/*135,136    */139399    , cNEqual     ,PositionalParams,0}},
        /* 184:	@R @F (cNEqual [(cMul % <1>) &])
         *	:	(cMul  <1>) DIV( & % )@C
         */		 {ReplaceParams , 18, 2,/*246,517    */529654    , {2,/*273,24     */24849     , cNEqual     ,PositionalParams,0}},
        /* 185:	@R @F (cNEqual [(cPow [x %@P]) &])
         *	:	(cPow [(cPow [x %]) /%@C]) POW( & /%@C )@C
         */		 {ReplaceParams , 18, 2,/*388,399    */408964    , {2,/*387,24     */24963     , cNEqual     ,PositionalParams,0}},
        /* 186:	@R @F (cNEqual [(cMul % <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul DIV( & % )@C <2>)
         */		 {ReplaceParams , 18, 2,/*246,276    */282870    , {2,/*273,275    */281873    , cNEqual     ,PositionalParams,0}},
        /* 187:	@R @F (cNEqual [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , 18, 2,/*0,7        */7168      , {2,/*390,392    */401798    , cNEqual     ,PositionalParams,0}},
        /* 188:	@R @F (cNEqual [&@P (cPow [%@P x])])
         *	:	DIV( LOG( & )@C LOG( % )@C )@C x
         */		 {ReplaceParams , 18, 2,/*518,0      */518       , {2,/*28,391     */400412    , cNEqual     ,PositionalParams,0}},
        /* 189:	@R (cLess [(cAdd % <1>) &])
         *	:	(cAdd  <1>) SUB( & % )@C
         */		 {ReplaceParams , 16, 2,/*120,516    */528504    , {2,/*137,24     */24713     , cLess       ,PositionalParams,0}},
        /* 190:	@R (cLess [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -%@C <2>)
         */		 {ReplaceParams , 16, 2,/*120,139    */142456    , {2,/*137,138    */141449    , cLess       ,PositionalParams,0}},
        /* 191:	@R (cLess [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , 16, 2,/*120,121    */124024    , {2,/*135,136    */139399    , cLess       ,PositionalParams,0}},
        /* 192:	@R @F (cLess [x 0.5])
         *	->	(cAbsNot [x])
         */		 {ProduceNewTree, 18, 1,/*571        */571       , {2,/*0,45       */46080     , cLess       ,PositionalParams,0}},
        /* 193:	@R @F (cLess [(cMul %@P <1>) &])
         *	:	(cMul  <1>) DIV( & % )@C
         */		 {ReplaceParams , 18, 2,/*246,517    */529654    , {2,/*256,24     */24832     , cLess       ,PositionalParams,0}},
        /* 194:	@R @F (cLess [(cMul %@N <1>) &])
         *	:	DIV( & % )@C (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*517,246    */252421    , {2,/*254,24     */24830     , cLess       ,PositionalParams,0}},
        /* 195:	@R @F (cLess [(cPow [x %@P]) &])
         *	:	(cPow [(cPow [x %]) /%@C]) POW( & /%@C )@C
         */		 {ReplaceParams , 18, 2,/*388,399    */408964    , {2,/*387,24     */24963     , cLess       ,PositionalParams,0}},
        /* 196:	@R @F (cLess [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul DIV( & % )@C <2>)
         */		 {ReplaceParams , 18, 2,/*246,276    */282870    , {2,/*256,275    */281856    , cLess       ,PositionalParams,0}},
        /* 197:	@R @F (cLess [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul DIV( & % )@C <2>) (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*276,246    */252180    , {2,/*254,275    */281854    , cLess       ,PositionalParams,0}},
        /* 198:	@R @F (cLess [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , 18, 2,/*0,7        */7168      , {2,/*390,392    */401798    , cLess       ,PositionalParams,0}},
        /* 199:	@R @F (cLess [&@P (cPow [%@P x])])
         *	:	DIV( LOG( & )@C LOG( % )@C )@C x
         */		 {ReplaceParams , 18, 2,/*518,0      */518       , {2,/*28,391     */400412    , cLess       ,PositionalParams,0}},
        /* 200:	@R @I (cLess [0 (cAbs [x])])
         *	->	(cNotNot [x])
         */		 {ProduceNewTree, 20, 1,/*562        */562       , {2,/*41,400     */409641    , cLess       ,PositionalParams,0}},
        /* 201:	@R (cLessOrEq [(cAdd % <1>) &])
         *	:	(cAdd  <1>) SUB( & % )@C
         */		 {ReplaceParams , 16, 2,/*120,516    */528504    , {2,/*137,24     */24713     , cLessOrEq   ,PositionalParams,0}},
        /* 202:	@R (cLessOrEq [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -%@C <2>)
         */		 {ReplaceParams , 16, 2,/*120,139    */142456    , {2,/*137,138    */141449    , cLessOrEq   ,PositionalParams,0}},
        /* 203:	@R (cLessOrEq [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , 16, 2,/*120,121    */124024    , {2,/*135,136    */139399    , cLessOrEq   ,PositionalParams,0}},
        /* 204:	@R @F (cLessOrEq [% (cAbs [x])])
         *	->	(cNotNot [(cMul {x 0.5 /%@C})])
         */		 {ProduceNewTree, 18, 1,/*565        */565       , {2,/*15,400     */409615    , cLessOrEq   ,PositionalParams,0}},
        /* 205:	@R @F (cLessOrEq [(cMul %@P <1>) &])
         *	:	(cMul  <1>) DIV( & % )@C
         */		 {ReplaceParams , 18, 2,/*246,517    */529654    , {2,/*256,24     */24832     , cLessOrEq   ,PositionalParams,0}},
        /* 206:	@R @F (cLessOrEq [(cMul %@N <1>) &])
         *	:	DIV( & % )@C (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*517,246    */252421    , {2,/*254,24     */24830     , cLessOrEq   ,PositionalParams,0}},
        /* 207:	@R @F (cLessOrEq [(cPow [x %@P]) &])
         *	:	(cPow [(cPow [x %]) /%@C]) POW( & /%@C )@C
         */		 {ReplaceParams , 18, 2,/*388,399    */408964    , {2,/*387,24     */24963     , cLessOrEq   ,PositionalParams,0}},
        /* 208:	@R @F (cLessOrEq [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul DIV( & % )@C <2>)
         */		 {ReplaceParams , 18, 2,/*246,276    */282870    , {2,/*256,275    */281856    , cLessOrEq   ,PositionalParams,0}},
        /* 209:	@R @F (cLessOrEq [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul DIV( & % )@C <2>) (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*276,246    */252180    , {2,/*254,275    */281854    , cLessOrEq   ,PositionalParams,0}},
        /* 210:	@R @F (cLessOrEq [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , 18, 2,/*0,7        */7168      , {2,/*390,392    */401798    , cLessOrEq   ,PositionalParams,0}},
        /* 211:	@R @F (cLessOrEq [&@P (cPow [%@P x])])
         *	:	DIV( LOG( & )@C LOG( % )@C )@C x
         */		 {ReplaceParams , 18, 2,/*518,0      */518       , {2,/*28,391     */400412    , cLessOrEq   ,PositionalParams,0}},
        /* 212:	@R @I (cLessOrEq [1 (cAbs [x])])
         *	->	(cNotNot [x])
         */		 {ProduceNewTree, 20, 1,/*562        */562       , {2,/*47,400     */409647    , cLessOrEq   ,PositionalParams,0}},
        /* 213:	@R (cGreater [(cAdd % <1>) &])
         *	:	(cAdd  <1>) SUB( & % )@C
         */		 {ReplaceParams , 16, 2,/*120,516    */528504    , {2,/*137,24     */24713     , cGreater    ,PositionalParams,0}},
        /* 214:	@R (cGreater [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -%@C <2>)
         */		 {ReplaceParams , 16, 2,/*120,139    */142456    , {2,/*137,138    */141449    , cGreater    ,PositionalParams,0}},
        /* 215:	@R (cGreater [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , 16, 2,/*120,121    */124024    , {2,/*135,136    */139399    , cGreater    ,PositionalParams,0}},
        /* 216:	@R @F (cGreater [% (cAbs [x])])
         *	->	(cNot [(cMul {x 0.5 /%@C})])
         */		 {ProduceNewTree, 18, 1,/*539        */539       , {2,/*15,400     */409615    , cGreater    ,PositionalParams,0}},
        /* 217:	@R @F (cGreater [(cMul %@P <1>) &])
         *	:	(cMul  <1>) DIV( & % )@C
         */		 {ReplaceParams , 18, 2,/*246,517    */529654    , {2,/*256,24     */24832     , cGreater    ,PositionalParams,0}},
        /* 218:	@R @F (cGreater [(cMul %@N <1>) &])
         *	:	DIV( & % )@C (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*517,246    */252421    , {2,/*254,24     */24830     , cGreater    ,PositionalParams,0}},
        /* 219:	@R @F (cGreater [(cPow [x %@P]) &])
         *	:	(cPow [(cPow [x %]) /%@C]) POW( & /%@C )@C
         */		 {ReplaceParams , 18, 2,/*388,399    */408964    , {2,/*387,24     */24963     , cGreater    ,PositionalParams,0}},
        /* 220:	@R @F (cGreater [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul DIV( & % )@C <2>)
         */		 {ReplaceParams , 18, 2,/*246,276    */282870    , {2,/*256,275    */281856    , cGreater    ,PositionalParams,0}},
        /* 221:	@R @F (cGreater [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul DIV( & % )@C <2>) (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*276,246    */252180    , {2,/*254,275    */281854    , cGreater    ,PositionalParams,0}},
        /* 222:	@R @F (cGreater [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , 18, 2,/*0,7        */7168      , {2,/*390,392    */401798    , cGreater    ,PositionalParams,0}},
        /* 223:	@R @F (cGreater [&@P (cPow [%@P x])])
         *	:	DIV( LOG( & )@C LOG( % )@C )@C x
         */		 {ReplaceParams , 18, 2,/*518,0      */518       , {2,/*28,391     */400412    , cGreater    ,PositionalParams,0}},
        /* 224:	@R @I (cGreater [1 (cAbs [x])])
         *	->	(cNot [x])
         */		 {ProduceNewTree, 20, 1,/*538        */538       , {2,/*47,400     */409647    , cGreater    ,PositionalParams,0}},
        /* 225:	@R (cGreaterOrEq [(cAdd % <1>) &])
         *	:	(cAdd  <1>) SUB( & % )@C
         */		 {ReplaceParams , 16, 2,/*120,516    */528504    , {2,/*137,24     */24713     , cGreaterOrEq,PositionalParams,0}},
        /* 226:	@R (cGreaterOrEq [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -%@C <2>)
         */		 {ReplaceParams , 16, 2,/*120,139    */142456    , {2,/*137,138    */141449    , cGreaterOrEq,PositionalParams,0}},
        /* 227:	@R (cGreaterOrEq [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , 16, 2,/*120,121    */124024    , {2,/*135,136    */139399    , cGreaterOrEq,PositionalParams,0}},
        /* 228:	@R @F (cGreaterOrEq [x 0.5])
         *	->	(cAbsNotNot [x])
         */		 {ProduceNewTree, 18, 1,/*572        */572       , {2,/*0,45       */46080     , cGreaterOrEq,PositionalParams,0}},
        /* 229:	@R @F (cGreaterOrEq [(cMul %@P <1>) &])
         *	:	(cMul  <1>) DIV( & % )@C
         */		 {ReplaceParams , 18, 2,/*246,517    */529654    , {2,/*256,24     */24832     , cGreaterOrEq,PositionalParams,0}},
        /* 230:	@R @F (cGreaterOrEq [(cMul %@N <1>) &])
         *	:	DIV( & % )@C (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*517,246    */252421    , {2,/*254,24     */24830     , cGreaterOrEq,PositionalParams,0}},
        /* 231:	@R @F (cGreaterOrEq [(cPow [x %@P]) &])
         *	:	(cPow [(cPow [x %]) /%@C]) POW( & /%@C )@C
         */		 {ReplaceParams , 18, 2,/*388,399    */408964    , {2,/*387,24     */24963     , cGreaterOrEq,PositionalParams,0}},
        /* 232:	@R @F (cGreaterOrEq [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul DIV( & % )@C <2>)
         */		 {ReplaceParams , 18, 2,/*246,276    */282870    , {2,/*256,275    */281856    , cGreaterOrEq,PositionalParams,0}},
        /* 233:	@R @F (cGreaterOrEq [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul DIV( & % )@C <2>) (cMul  <1>)
         */		 {ReplaceParams , 18, 2,/*276,246    */252180    , {2,/*254,275    */281854    , cGreaterOrEq,PositionalParams,0}},
        /* 234:	@R @F (cGreaterOrEq [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , 18, 2,/*0,7        */7168      , {2,/*390,392    */401798    , cGreaterOrEq,PositionalParams,0}},
        /* 235:	@R @F (cGreaterOrEq [&@P (cPow [%@P x])])
         *	:	DIV( LOG( & )@C LOG( % )@C )@C x
         */		 {ReplaceParams , 18, 2,/*518,0      */518       , {2,/*28,391     */400412    , cGreaterOrEq,PositionalParams,0}},
        /* 236:	@R @I (cGreaterOrEq [0 (cAbs [x])])
         *	->	(cNot [x])
         */		 {ProduceNewTree, 20, 1,/*538        */538       , {2,/*41,400     */409641    , cGreaterOrEq,PositionalParams,0}},
        /* 237:	@I (cNot [(cAdd % <1>)])
         *	->	(cEqual [-%@C (cAdd  <1>)])
         */		 {ProduceNewTree, 4, 1,/*519        */519       , {1,/*137        */137       , cNot        ,PositionalParams,0}},
        /* 238:	@R (cNot [x@P])
         *	->	(cAbsNot [x])
         */		 {ProduceNewTree, 16, 1,/*571        */571       , {1,/*2          */2         , cNot        ,PositionalParams,0}},
        /* 239:	(cAnd (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cAnd {y a}) (cAnd {z b})])
         */		 {ReplaceParams , 0, 1,/*452        */452       , {2,/*438,440    */450998    , cAnd        ,AnyParams       ,0}},
        /* 240:	(cAnd (cEqual [x y])@D12 (cEqual [y z])@D24 (cEqual [x z])@D20)
         *	:	(cEqual [x y]) (cEqual [y z])
         */		 {ReplaceParams , 0, 2,/*521,524    */537097    , {3,/*520,523,522*/547892744 , cAnd        ,AnyParams       ,0}},
        /* 241:	@R (cAnd x@L <1>)
         *	->	(cNotNot [(cMul {x (cAnd  <1>)})])
         */		 {ProduceNewTree, 16, 1,/*566        */566       , {1,/*5          */5         , cAnd        ,AnyParams       ,1}},
        /* 242:	@R (cAnd x@P y@P)
         *	:	(cAbsAnd {x y})
         */		 {ReplaceParams , 16, 1,/*569        */569       , {2,/*2,13       */13314     , cAnd        ,AnyParams       ,0}},
        /* 243:	@R (cAnd (cNot [x]) (cNot [y]))
         *	:	(cNot [(cOr {x y})])
         */		 {ReplaceParams , 16, 1,/*544        */544       , {2,/*538,540    */553498    , cAnd        ,AnyParams       ,0}},
        /* 244:	@R (cAnd (cNot [z]) (cIf [x (cNot [y]) %@L]))
         *	:	(cNot [(cOr {z (cIf [x y (cNot [%])])})])
         */		 {ReplaceParams , 16, 1,/*546        */546       , {2,/*545,451    */462369    , cAnd        ,AnyParams       ,0}},
        /* 245:	@R (cAnd (cNot [z]) (cIf [x %@L (cNot [y])]))
         *	:	(cNot [(cOr {z (cIf [x (cNot [%]) y])})])
         */		 {ReplaceParams , 16, 1,/*548        */548       , {2,/*545,455    */466465    , cAnd        ,AnyParams       ,0}},
        /* 246:	(cOr (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cOr {y a}) (cOr {z b})])
         */		 {ReplaceParams , 0, 1,/*457        */457       , {2,/*438,440    */450998    , cOr         ,AnyParams       ,0}},
        /* 247:	@R (cOr x@P y@P)
         *	:	(cAbsOr {x y})
         */		 {ReplaceParams , 16, 1,/*570        */570       , {2,/*2,13       */13314     , cOr         ,AnyParams       ,0}},
        /* 248:	@R (cOr x@L y@L)
         *	:	(cNotNot [(cAdd {x y})])
         */		 {ReplaceParams , 16, 1,/*563        */563       , {2,/*5,8        */8197      , cOr         ,AnyParams       ,0}},
        /* 249:	@R (cOr (cNot [x]) (cNot [y]))
         *	:	(cNot [(cAnd {x y})])
         */		 {ReplaceParams , 16, 1,/*541        */541       , {2,/*538,540    */553498    , cOr         ,AnyParams       ,0}},
        /* 250:	@R (cOr (cNot [z]) (cIf [x (cNot [y]) %@L]))
         *	:	(cNot [(cAnd {z (cIf [x y (cNot [%])])})])
         */		 {ReplaceParams , 16, 1,/*542        */542       , {2,/*545,451    */462369    , cOr         ,AnyParams       ,0}},
        /* 251:	@R (cOr (cNot [z]) (cIf [x %@L (cNot [y])]))
         *	:	(cNot [(cAnd {z (cIf [x (cNot [%]) y])})])
         */		 {ReplaceParams , 16, 1,/*543        */543       , {2,/*545,455    */466465    , cOr         ,AnyParams       ,0}},
        /* 252:	@R (cOr x@L (cAdd  <1>)@P)
         *	:	(cNotNot [(cAdd x <1>)])
         */		 {ReplaceParams , 16, 1,/*564        */564       , {2,/*5,140      */143365    , cOr         ,AnyParams       ,0}},
        /* 253:	@I (cNotNot [(cAdd % <1>)])
         *	->	(cNEqual [-%@C (cAdd  <1>)])
         */		 {ProduceNewTree, 4, 1,/*525        */525       , {1,/*137        */137       , cNotNot     ,PositionalParams,0}},
        /* 254:	@R (cNotNot [x@P])
         *	->	(cAbsNotNot [x])
         */		 {ProduceNewTree, 16, 1,/*572        */572       , {1,/*2          */2         , cNotNot     ,PositionalParams,0}},
        /* 255:	@R @L (cNotNot [x])
         *	->	x
         */		 {ProduceNewTree, 17, 1,/*0          */0         , {1,/*0          */0         , cNotNot     ,PositionalParams,0}},
        /* 256:	@R @F (cAbsNotNot (cMul %@P <1>))
         *	->	(cGreaterOrEq [(cMul  <1>) MUL( 0.5 /%@C )@C])
         */		 {ProduceNewTree, 18, 1,/*537        */537       , {1,/*256        */256       , cAbsNotNot  ,AnyParams       ,0}},
        /* 257:	@R @F (cAbsNotNot (cMul %@N <1>))
         *	->	(cLessOrEq [(cMul  <1>) MUL( 0.5 /%@C )@C])
         */		 {ProduceNewTree, 18, 1,/*531        */531       , {1,/*254        */254       , cAbsNotNot  ,AnyParams       ,0}},
        /* 258:	(cAbsIf [x 1 0])
         *	->	(cAbsNotNot [x])
         */		 {ProduceNewTree, 0, 1,/*572        */572       , {3,/*0,47,41    */43039744  , cAbsIf      ,PositionalParams,0}},
        /* 259:	(cAbsIf [x 0 1])
         *	->	(cAbsNot [x])
         */		 {ProduceNewTree, 0, 1,/*571        */571       , {3,/*0,41,47    */49325056  , cAbsIf      ,PositionalParams,0}},
        /* 260:	@R (cAbsIf [(cNotNot [x]) y z])
         *	->	(cIf [x y z])
         */		 {ProduceNewTree, 16, 1,/*454        */454       , {3,/*562,7,31   */32513586  , cAbsIf      ,PositionalParams,0}},
        /* 261:	@R (cAbsIf [(cLessOrEq [x y]) z a])
         *	:	(cLess [y x]) a z
         */		 {ReplaceParams , 16, 3,/*529,35,31  */32542225  , {3,/*530,31,35  */36732434  , cAbsIf      ,PositionalParams,0}},
    };

    struct grammar_optimize_abslogical_type
    {
        unsigned c;
        unsigned short l[9];
    };
    extern "C"
    {
        grammar_optimize_abslogical_type grammar_optimize_abslogical =
        {
            9,
            { 34,192,228,238,242,247,254,260,261
    }   };  }
    struct grammar_optimize_ignore_if_sideeffects_type
    {
        unsigned c;
        unsigned short l[59];
    };
    extern "C"
    {
        grammar_optimize_ignore_if_sideeffects_type grammar_optimize_ignore_if_sideeffects =
        {
            59,
            { 0,20,21,22,23,24,25,26,27,28,
              29,30,31,32,33,35,36,41,42,43,
              44,78,79,122,123,160,161,163,164,165,
              166,167,168,169,178,179,180,200,204,212,
              216,224,236,237,239,240,243,244,245,246,
              249,250,251,253,255,256,257,258,259
    }   };  }
    struct grammar_optimize_nonshortcut_logical_evaluation_type
    {
        unsigned c;
        unsigned short l[56];
    };
    extern "C"
    {
        grammar_optimize_nonshortcut_logical_evaluation_type grammar_optimize_nonshortcut_logical_evaluation =
        {
            56,
            { 0,25,27,28,29,30,31,32,33,35,
              36,41,42,43,44,78,79,122,123,160,
              161,163,164,165,166,167,168,169,178,179,
              180,200,204,212,216,224,236,237,239,240,
              241,243,244,245,246,248,249,250,251,252,
              253,255,256,257,258,259
    }   };  }
    struct grammar_optimize_recreate_type
    {
        unsigned c;
        unsigned short l[22];
    };
    extern "C"
    {
        grammar_optimize_recreate_type grammar_optimize_recreate =
        {
            22,
            { 18,55,56,57,80,81,82,83,84,85,
              117,118,120,121,130,131,132,133,134,135,
              136,137
    }   };  }
    struct grammar_optimize_round1_type
    {
        unsigned c;
        unsigned short l[125];
    };
    extern "C"
    {
        grammar_optimize_round1_type grammar_optimize_round1 =
        {
            125,
            { 0,1,2,3,4,5,6,7,8,9,
              10,11,12,13,14,19,25,27,28,29,
              30,31,32,33,35,36,37,38,41,42,
              43,44,45,46,47,48,49,50,51,52,
              53,54,58,59,60,61,62,63,64,65,
              66,67,68,69,70,71,78,79,80,81,
              82,83,84,85,86,87,88,93,94,95,
              96,97,98,99,100,101,117,118,119,120,
              121,122,123,124,125,126,127,128,129,138,
              160,161,162,163,164,165,166,167,168,169,
              178,179,180,200,204,212,216,224,236,237,
              239,240,243,244,245,246,249,250,251,253,
              255,256,257,258,259
    }   };  }
    struct grammar_optimize_round2_type
    {
        unsigned c;
        unsigned short l[103];
    };
    extern "C"
    {
        grammar_optimize_round2_type grammar_optimize_round2 =
        {
            103,
            { 0,15,16,17,25,27,28,29,30,31,
              32,33,35,36,39,40,41,42,43,44,
              45,46,47,48,49,50,51,52,53,54,
              59,60,72,73,78,79,86,87,88,89,
              90,91,92,102,103,104,105,106,107,108,
              109,110,111,112,113,114,115,116,119,122,
              123,124,125,126,127,128,139,159,160,161,
              162,163,164,165,166,167,168,169,178,179,
              180,200,204,212,216,224,236,237,239,240,
              243,244,245,246,249,250,251,253,255,256,
              257,258,259
    }   };  }
    struct grammar_optimize_round3_type
    {
        unsigned c;
        unsigned short l[79];
    };
    extern "C"
    {
        grammar_optimize_round3_type grammar_optimize_round3 =
        {
            79,
            { 74,75,76,77,140,141,142,143,144,145,
              146,147,148,149,150,151,152,153,154,155,
              156,157,158,170,171,172,173,174,175,176,
              177,181,182,183,184,185,186,187,188,189,
              190,191,193,194,195,196,197,198,199,201,
              202,203,205,206,207,208,209,210,211,213,
              214,215,217,218,219,220,221,222,223,225,
              226,227,229,230,231,232,233,234,235
    }   };  }
    struct grammar_optimize_round4_type
    {
        unsigned c;
        unsigned short l[12];
    };
    extern "C"
    {
        grammar_optimize_round4_type grammar_optimize_round4 =
        {
            12,
            { 18,55,56,57,130,131,132,133,134,135,
              136,137
    }   };  }
    struct grammar_optimize_shortcut_logical_evaluation_type
    {
        unsigned c;
        unsigned short l[53];
    };
    extern "C"
    {
        grammar_optimize_shortcut_logical_evaluation_type grammar_optimize_shortcut_logical_evaluation =
        {
            53,
            { 0,25,27,28,29,30,31,32,33,35,
              36,41,42,43,44,78,79,122,123,160,
              161,163,164,165,166,167,168,169,178,179,
              180,200,204,212,216,224,236,237,239,240,
              243,244,245,246,249,250,251,253,255,256,
              257,258,259
    }   };  }
}
namespace FPoptimizer_Grammar
{
    template<typename Value_t>
    ParamSpec ParamSpec_Extract(unsigned paramlist, unsigned index)
    {
        index = (paramlist >> (index * 10)) & 1023 /* % (1 << 10) */;
        if(index >= 57)
            return ParamSpec(SubFunction,(const void*)&plist_s[index-57]);
        if(index >= 37)
            return ParamSpec(NumConstant,(const void*)&plist_n_container<Value_t>::plist_n[index-37]);
        return ParamSpec(ParamHolder,(const void*)&plist_p[index]);
    }
}
/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_Grammar
{
#define FP_INSTANTIATE(type) \
    template ParamSpec ParamSpec_Extract<type>(unsigned paramlist, unsigned index);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#line 1 "fpoptimizer/optimize.cc"
// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

// line removed for fpoptimizer.cc: #include "grammar.hh"
// line removed for fpoptimizer.cc: #include "consts.hh"
// line removed for fpoptimizer.cc: #include "opcodename.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"

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
            if(found.specs.isNull())
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

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_Optimize
{
#define FP_INSTANTIATE(type) \
    template void ApplyGrammars(FPoptimizer_CodeTree::CodeTree<type>& tree);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/optimize_match.cc"
// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <cmath>

#include <memory> /* for auto_ptr */

// line removed for fpoptimizer.cc: #include "grammar.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"
// line removed for fpoptimizer.cc: #include "rangeestimation.hh"
// line removed for fpoptimizer.cc: #include "consts.hh"

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_Grammar;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace
{
    /* Test the given constraints to a given CodeTree */
    template<typename Value_t>
    bool TestImmedConstraints(unsigned bitmask, const CodeTree<Value_t>& tree)
    {
        switch(bitmask & ValueMask)
        {
            case Value_AnyNum: case ValueMask: break;
            case Value_EvenInt:
                if(GetEvennessInfo(tree) != IsAlways)
                    return false;
                break;
            case Value_OddInt:
                if(GetEvennessInfo(tree) != IsNever)
                    return false;
                break;
            case Value_IsInteger:
                if(GetIntegerInfo(tree) != IsAlways) return false;
                break;
            case Value_NonInteger:
                if(GetIntegerInfo(tree) != IsNever) return false;
                break;
            case Value_Logical:
                if(!IsLogicalValue(tree)) return false;
                break;
        }
        switch(bitmask & SignMask)
        {
            case Sign_AnySign: /*case SignMask:*/ break;
            case Sign_Positive:
                if(GetPositivityInfo(tree) != IsAlways) return false;
                break;
            case Sign_Negative:
                if(GetPositivityInfo(tree) != IsNever) return false;
                break;
            case Sign_NoIdea:
                if(GetPositivityInfo(tree) != Unknown) return false;
                break;
        }
        switch(bitmask & OnenessMask)
        {
            case Oneness_Any: case OnenessMask: break;
            case Oneness_One:
                if(!tree.IsImmed()) return false;
                if(!fp_equal(fp_abs(tree.GetImmed()), Value_t(1))) return false;
                break;
            case Oneness_NotOne:
                if(!tree.IsImmed()) return false;
                if(fp_equal(fp_abs(tree.GetImmed()), Value_t(1))) return false;
                break;
        }
        switch(bitmask & ConstnessMask)
        {
            case Constness_Any: /*case ConstnessMask:*/ break;
            case Constness_Const:
                if(!tree.IsImmed()) return false;
                break;
            case Constness_NotConst:
                if(tree.IsImmed()) return false;
                break;
        }
        return true;
    }

    template<unsigned extent, unsigned nbits, typename item_type=unsigned int>
    struct nbitmap
    {
    private:
        static const unsigned bits_in_char = 8;
        static const unsigned per_item = (sizeof(item_type)*bits_in_char)/nbits;
        item_type data[(extent+per_item-1) / per_item];
    public:
        void inc(unsigned index, int by=1)
        {
            data[pos(index)] += by * item_type(1 << shift(index));
        }
        inline void dec(unsigned index) { inc(index, -1); }
        int get(unsigned index) const { return (data[pos(index)] >> shift(index)) & mask(); }

        static inline unsigned pos(unsigned index) { return index/per_item; }
        static inline unsigned shift(unsigned index) { return nbits * (index%per_item); }
        static inline unsigned mask() { return (1 << nbits)-1; }
        static inline unsigned mask(unsigned index) { return mask() << shift(index); }
    };

    struct Needs
    {
        int SubTrees     : 8; // This many subtrees
        int Others       : 8; // This many others (namedholder)
        int minimum_need : 8; // At least this many leaves (restholder may require more)
        int Immeds       : 8; // This many immeds

        nbitmap<VarBegin,2> SubTreesDetail; // This many subtrees of each opcode type

        Needs()
        {
            std::memset(this, 0, sizeof(*this));
        }
        Needs(const Needs& b)
        {
            std::memcpy(this, &b, sizeof(b));
        }
        Needs& operator= (const Needs& b)
        {
            std::memcpy(this, &b, sizeof(b));
            return *this;
        }
    };

    template<typename Value_t>
    Needs CreateNeedList_uncached(const ParamSpec_SubFunctionData& params)
    {
        Needs NeedList;

        // Figure out what we need
        for(unsigned a = 0; a < params.param_count; ++a)
        {
            const ParamSpec& parampair = ParamSpec_Extract<Value_t>(params.param_list, a);
            switch(parampair.first)
            {
                case SubFunction:
                {
                    const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                    if(param.data.match_type == GroupFunction)
                        ++NeedList.Immeds;
                    else
                    {
                        ++NeedList.SubTrees;
                        assert( param.data.subfunc_opcode < VarBegin );
                        NeedList.SubTreesDetail.inc(param.data.subfunc_opcode);
                    }
                    ++NeedList.minimum_need;
                    break;
                }
                case NumConstant:
                case ParamHolder:
                    ++NeedList.Others;
                    ++NeedList.minimum_need;
                    break;
            }
        }

        return NeedList;
    }

    template<typename Value_t>
    Needs& CreateNeedList(const ParamSpec_SubFunctionData& params)
    {
        typedef std::map<const ParamSpec_SubFunctionData*, Needs> needlist_cached_t;
        static needlist_cached_t needlist_cached;

        needlist_cached_t::iterator i = needlist_cached.lower_bound(&params);
        if(i != needlist_cached.end() && i->first == &params)
            return i->second;

        return
            needlist_cached.insert(i,
                 std::make_pair(&params, CreateNeedList_uncached<Value_t> (params))
            )->second;
    }
    /* Construct CodeTree from a GroupFunction, hopefully evaluating to a constant value */

    template<typename Value_t>
    CodeTree<Value_t> CalculateGroupFunction(
        const ParamSpec& parampair,
        const MatchInfo<Value_t>& info)
    {
        switch( parampair.first )
        {
            case NumConstant:
            {
                const ParamSpec_NumConstant<Value_t>& param = *(const ParamSpec_NumConstant<Value_t>*) parampair.second;
                return CodeTreeImmed( param.constvalue ); // Note: calculates hash too.
            }
            case ParamHolder:
            {
                const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                return info.GetParamHolderValueIfFound( param.index );
                // If the ParamHolder is not defined, it will simply
                // return an Undefined tree. This is ok.
            }
            case SubFunction:
            {
                const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                /* Synthesize a CodeTree which will take care of
                 * constant-folding our expression. It will also
                 * indicate whether the result is, in fact,
                 * a constant at all. */
                CodeTree<Value_t> result;
                result.SetOpcode( param.data.subfunc_opcode );
                result.GetParams().reserve(param.data.param_count);
                for(unsigned a=0; a<param.data.param_count; ++a)
                {
                    CodeTree<Value_t> tmp(
                        CalculateGroupFunction
                        (ParamSpec_Extract<Value_t> (param.data.param_list, a), info)
                                );
                    result.AddParamMove(tmp);
                }
                result.Rehash(); // This will also call ConstantFolding().
                return result;
            }
        }
        // Issue an un-calculatable tree. (This should be unreachable)
        return CodeTree<Value_t>(); // cNop
    }
}

namespace FPoptimizer_Optimize
{
    /* Test the list of parameters to a given CodeTree */
    /* A helper function which simply checks whether the
     * basic shape of the tree matches what we are expecting
     * i.e. given number of numeric constants, etc.
     */
    template<typename Value_t>
    bool IsLogisticallyPlausibleParamsMatch(
        const ParamSpec_SubFunctionData& params,
        const CodeTree<Value_t>& tree)
    {
        /* First, check if the tree has any chances of matching... */
        /* Figure out what we need. */
        Needs NeedList ( CreateNeedList<Value_t> (params) );

        size_t nparams = tree.GetParamCount();

        if(nparams < size_t(NeedList.minimum_need))
        {
            // Impossible to satisfy
            return false;
        }

        // Figure out what we have (note: we already assume that the opcode of the tree matches!)
        for(size_t a=0; a<nparams; ++a)
        {
            unsigned opcode = tree.GetParam(a).GetOpcode();
            switch(opcode)
            {
                case cImmed:
                    if(NeedList.Immeds > 0) --NeedList.Immeds;
                    else --NeedList.Others;
                    break;
                case VarBegin:
                case cFCall:
                case cPCall:
                    --NeedList.Others;
                    break;
                default:
                    assert( opcode < VarBegin );
                    if(NeedList.SubTrees > 0
                    && NeedList.SubTreesDetail.get(opcode) > 0)
                    {
                        --NeedList.SubTrees;
                        NeedList.SubTreesDetail.dec(opcode);
                    }
                    else --NeedList.Others;
            }
        }

        // Check whether all needs were satisfied
        if(NeedList.Immeds > 0
        || NeedList.SubTrees > 0
        || NeedList.Others > 0)
        {
            // Something came short, impossible to satisfy.
            return false;
        }

        if(params.match_type != AnyParams)
        {
            if(0
            //|| NeedList.Immeds < 0 - already checked
            || NeedList.SubTrees < 0
            || NeedList.Others < 0
            //|| params.count != nparams - already checked
              )
            {
                // Something was too much.
                return false;
            }
        }
        return true;
    }

    /* Test the given parameter to a given CodeTree */
    template<typename Value_t>
    MatchResultType TestParam(
        const ParamSpec& parampair,
        const CodeTree<Value_t>& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo<Value_t>& info)
    {
        /*std::cout << "TestParam(";
        DumpParam(parampair);
        std::cout << ", ";
        DumpTree(tree);
        std::cout << ")\n";*/

        /* What kind of param are we expecting */
        switch( parampair.first )
        {
            case NumConstant: /* A particular numeric value */
            {
                const ParamSpec_NumConstant<Value_t>& param = *(const ParamSpec_NumConstant<Value_t>*) parampair.second;
                if(!tree.IsImmed()) return false;
                Value_t imm = tree.GetImmed();
                switch(param.modulo)
                {
                    case Modulo_None: break;
                    case Modulo_Radians:
                        imm = fp_mod(imm, fp_const_twopi<Value_t>());
                        if(imm < Value_t(0))
                            imm += fp_const_twopi<Value_t>();
                        if(imm > fp_const_pi<Value_t>())
                            imm -= fp_const_twopi<Value_t>();
                        break;
                }
                return fp_equal(imm, param.constvalue);
            }
            case ParamHolder: /* Any arbitrary node */
            {
                const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                if(!TestImmedConstraints(param.constraints, tree)) return false;
                return info.SaveOrTestParamHolder(param.index, tree);
            }
            case SubFunction:
            {
                const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                if(param.data.match_type == GroupFunction)
                { /* A constant value acquired from this formula */
                    if(!TestImmedConstraints(param.constraints, tree)) return false;
                    /* Construct the formula */
                    CodeTree<Value_t> grammar_func = CalculateGroupFunction(parampair, info);
        #ifdef DEBUG_SUBSTITUTIONS
                    DumpHashes(grammar_func);
                    std::cout << *(const void**)&grammar_func.GetImmed();
                    std::cout << "\n";
                    std::cout << *(const void**)&tree.GetImmed();
                    std::cout << "\n";
                    DumpHashes(tree);
                    std::cout << "Comparing ";
                    DumpTree(grammar_func);
                    std::cout << " and ";
                    DumpTree(tree);
                    std::cout << ": ";
                    std::cout << (grammar_func.IsIdenticalTo(tree) ? "true" : "false");
                    std::cout << "\n";
        #endif
                    /* Evaluate it and compare */
                    return grammar_func.IsIdenticalTo(tree);
                }
                else /* A subtree conforming these specs */
                {
                    if(start_at.isNull())
                    {
                        if(!TestImmedConstraints(param.constraints, tree)) return false;
                        if(tree.GetOpcode() != param.data.subfunc_opcode) return false;
                    }
                    return TestParams(param.data, tree, start_at, info, false);
                }
            }
        }
        return false;
    }

    template<typename Value_t>
    struct PositionalParams_Rec
    {
        MatchPositionSpecBaseP start_at; /* child's start_at */
        MatchInfo<Value_t>     info;     /* backup of "info" at start */

        PositionalParams_Rec(): start_at(), info() { }
    };

    template<typename Value_t>
    class MatchPositionSpec_PositionalParams
        : public MatchPositionSpecBase,
          public std::vector<PositionalParams_Rec<Value_t> >
    {
    public:
        explicit MatchPositionSpec_PositionalParams(size_t n)
            : MatchPositionSpecBase(),
              std::vector<PositionalParams_Rec<Value_t> > (n)
              { }
    };

    struct AnyWhere_Rec
    {
        MatchPositionSpecBaseP start_at; /* child's start_at */
        AnyWhere_Rec() : start_at() { }
    };
    class MatchPositionSpec_AnyWhere
        : public MatchPositionSpecBase,
          public std::vector<AnyWhere_Rec>
    {
    public:
        unsigned trypos;   /* which param index to try next */

        explicit MatchPositionSpec_AnyWhere(size_t n)
            : MatchPositionSpecBase(),
              std::vector<AnyWhere_Rec> (n),
              trypos(0)
              { }
    };

    template<typename Value_t>
    MatchResultType TestParam_AnyWhere(
        const ParamSpec& parampair,
        const CodeTree<Value_t>& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo<Value_t>& info,
        std::vector<bool>&  used,
        bool TopLevel)
    {
        FPOPT_autoptr<MatchPositionSpec_AnyWhere> position;
        unsigned a;
        if(!start_at.isNull())
        {
            position = (MatchPositionSpec_AnyWhere*) &*start_at;
            a = position->trypos;
            goto retry_anywhere_2;
        }
        else
        {
            position = new MatchPositionSpec_AnyWhere(tree.GetParamCount());
            a = 0;
        }
        for(; a < tree.GetParamCount(); ++a)
        {
            if(used[a]) continue;

        retry_anywhere:
          { MatchResultType r = TestParam(
                parampair,
                tree.GetParam(a),
                (*position)[a].start_at,
                info);

            (*position)[a].start_at = r.specs;
            if(r.found)
            {
                used[a]               = true; // matched
                if(TopLevel) info.SaveMatchedParamIndex(a);

                position->trypos = a; // in case of backtrack, try a again
                return MatchResultType(true, &*position);
            } }
        retry_anywhere_2:
            if(!(*position)[a].start_at.isNull()) // is there another try?
            {
                goto retry_anywhere;
            }
            // no, move on
        }
        return false;
    }

    template<typename Value_t>
    struct AnyParams_Rec
    {
        MatchPositionSpecBaseP start_at; /* child's start_at */
        MatchInfo<Value_t>     info;     /* backup of "info" at start */
        std::vector<bool>      used;     /* which params are remaining */

        explicit AnyParams_Rec(size_t nparams)
            : start_at(), info(), used(nparams) { }
    };
    template<typename Value_t>
    class MatchPositionSpec_AnyParams
        : public MatchPositionSpecBase,
          public std::vector<AnyParams_Rec<Value_t> >
    {
    public:
        explicit MatchPositionSpec_AnyParams(size_t n, size_t m)
            : MatchPositionSpecBase(),
              std::vector<AnyParams_Rec<Value_t> > (n, AnyParams_Rec<Value_t>(m))
              { }
    };

    /* Test the list of parameters to a given CodeTree */
    template<typename Value_t>
    MatchResultType TestParams(
        const ParamSpec_SubFunctionData& model_tree,
        const CodeTree<Value_t>& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo<Value_t>& info,
        bool TopLevel)
    {
        /* When PositionalParams or SelectedParams, verify that
         * the number of parameters is exactly as expected.
         */
        if(model_tree.match_type != AnyParams)
        {
            if(model_tree.param_count != tree.GetParamCount())
                return false;
        }

        /* Verify that the tree basically conforms the shape we are expecting */
        /* This test is not necessary; it may just save us some work. */
        if(!IsLogisticallyPlausibleParamsMatch(model_tree, tree))
        {
            return false;
        }

        /* Verify each parameter that they are found in the tree as expected. */
        switch(model_tree.match_type)
        {
            case PositionalParams:
            {
                /* Simple: Test all given parameters in succession. */
                FPOPT_autoptr<MatchPositionSpec_PositionalParams<Value_t> > position;
                unsigned a;
                if(!start_at.isNull())
                {
                    position = (MatchPositionSpec_PositionalParams<Value_t> *) &*start_at;
                    a = model_tree.param_count - 1;
                    goto retry_positionalparams_2;
                }
                else
                {
                    position = new MatchPositionSpec_PositionalParams<Value_t> (model_tree.param_count);
                    a = 0;
                }

                for(; a < model_tree.param_count; ++a)
                {
                    (*position)[a].info = info;
                retry_positionalparams:
                  { MatchResultType r = TestParam(
                        ParamSpec_Extract<Value_t>(model_tree.param_list, a),
                        tree.GetParam(a),
                        (*position)[a].start_at,
                        info);

                    (*position)[a].start_at = r.specs;
                    if(r.found)
                    {
                        continue;
                  } }
                retry_positionalparams_2:
                    // doesn't match
                    if(!(*position)[a].start_at.isNull()) // is there another try?
                    {
                        info = (*position)[a].info;
                        goto retry_positionalparams;
                    }
                    // no, backtrack
                    if(a > 0)
                    {
                        --a;
                        goto retry_positionalparams_2;
                    }
                    // cannot backtrack
                    info = (*position)[0].info;
                    return false;
                }
                if(TopLevel)
                    for(unsigned a = 0; a < model_tree.param_count; ++a)
                        info.SaveMatchedParamIndex(a);
                return MatchResultType(true, &*position);
            }
            case SelectedParams:
                // same as AnyParams, except that model_tree.count==tree.GetParamCount()
                //                       and that there are no RestHolders
            case AnyParams:
            {
                /* Ensure that all given parameters are found somewhere, in any order */

                FPOPT_autoptr<MatchPositionSpec_AnyParams<Value_t> > position;
                std::vector<bool> used( tree.GetParamCount() );
                std::vector<unsigned> depcodes( model_tree.param_count );
                std::vector<unsigned> test_order( model_tree.param_count );
                for(unsigned a=0; a<model_tree.param_count; ++a)
                {
                    const ParamSpec parampair = ParamSpec_Extract<Value_t>(model_tree.param_list, a);
                    depcodes[a] = ParamSpec_GetDepCode(parampair);
                }
                { unsigned b=0;
                for(unsigned a=0; a<model_tree.param_count; ++a)
                    if(depcodes[a] != 0)
                        test_order[b++] = a;
                for(unsigned a=0; a<model_tree.param_count; ++a)
                    if(depcodes[a] == 0)
                        test_order[b++] = a;
                }

                unsigned a;
                if(!start_at.isNull())
                {
                    position = (MatchPositionSpec_AnyParams<Value_t>*) &*start_at;
                    if(model_tree.param_count == 0)
                    {
                        a = 0;
                        goto retry_anyparams_4;
                    }
                    a = model_tree.param_count - 1;
                    goto retry_anyparams_2;
                }
                else
                {
                    position = new MatchPositionSpec_AnyParams<Value_t>
                        (model_tree.param_count, tree.GetParamCount());
                    a = 0;
                    if(model_tree.param_count != 0)
                    {
                        (*position)[0].info   = info;
                        (*position)[0].used   = used;
                    }
                }
                // Match all but restholders
                for(; a < model_tree.param_count; ++a)
                {
                    if(a > 0) // this test is not necessary, but it saves from doing
                    {         // duplicate work, because [0] was already saved above.
                        (*position)[a].info   = info;
                        (*position)[a].used   = used;
                    }
                retry_anyparams:
                  { MatchResultType r = TestParam_AnyWhere<Value_t>(
                        ParamSpec_Extract<Value_t>(model_tree.param_list, test_order[a]),
                        tree,
                        (*position)[a].start_at,
                        info,
                        used,
                        TopLevel);
                    (*position)[a].start_at = r.specs;
                    if(r.found)
                    {
                        continue;
                  } }
                retry_anyparams_2:
                    // doesn't match
                    if(!(*position)[a].start_at.isNull()) // is there another try?
                    {
                        info = (*position)[a].info;
                        used = (*position)[a].used;
                        goto retry_anyparams;
                    }
                    // no, backtrack
                retry_anyparams_3:
                    if(a > 0)
                    {
                        --a;
                        goto retry_anyparams_2;
                    }
                    // cannot backtrack
                    info = (*position)[0].info;
                    return false;
                }
            retry_anyparams_4:
                // Capture anything remaining in the restholder
                if(model_tree.restholder_index != 0)
                {
                    //std::vector<bool> used_backup(used);
                    //MatchInfo         info_backup(info);

                    if(!TopLevel
                    || !info.HasRestHolder(model_tree.restholder_index))
                    {
                        std::vector<CodeTree<Value_t> > matches;
                        matches.reserve(tree.GetParamCount());
                        for(unsigned b = 0; b < tree.GetParamCount(); ++b)
                        {
                            if(used[b]) continue; // Ignore subtrees that were already used
                            // Save this tree to this restholder

                            matches.push_back(tree.GetParam(b));
                            used[b] = true;
                            if(TopLevel) info.SaveMatchedParamIndex(b);
                        }
                        if(!info.SaveOrTestRestHolder(model_tree.restholder_index, matches))
                        {
                            // Failure at restholder matching. Backtrack if possible.
                            //used.swap(used_backup);
                            //info.swap(info_backup);
                            goto retry_anyparams_3;
                        }
                        //std::cout << "Saved restholder " << model_tree.restholder_index << "\n";
                    }
                    else
                    {
                        const std::vector<CodeTree<Value_t> >& matches
                            = info.GetRestHolderValues(model_tree.restholder_index);
                        //std::cout << "Testing restholder " << model_tree.restholder_index << std::flush;
                        for(size_t a=0; a<matches.size(); ++a)
                        {
                            bool found = false;
                            for(unsigned b = 0; b < tree.GetParamCount(); ++b)
                            {
                                if(used[b]) continue;
                                if(matches[a].IsIdenticalTo(tree.GetParam(b)))
                                {
                                    used[b] = true;
                                    if(TopLevel) info.SaveMatchedParamIndex(b);
                                    found = true;
                                    break;
                                }
                            }
                            if(!found)
                            {
                                //std::cout << " ... failed\n";
                                // Failure at restholder matching. Backtrack if possible.
                                //used.swap(used_backup);
                                //info.swap(info_backup);
                                goto retry_anyparams_3;
                            }
                        }
                        //std::cout << " ... ok\n";
                    }
                }
                return MatchResultType(true, model_tree.param_count ? &*position : 0);
            }
            case GroupFunction: // never occurs
                break;
        }
        return false; // doesn't match
    }
}


/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_Optimize
{
#define FP_INSTANTIATE(type) \
    template \
    MatchResultType TestParams( \
        const ParamSpec_SubFunctionData& model_tree, \
        const CodeTree<type> & tree, \
        const MatchPositionSpecBaseP& start_at, \
        MatchInfo<type>& info, \
        bool TopLevel); \
    template \
    bool IsLogisticallyPlausibleParamsMatch( \
        const ParamSpec_SubFunctionData& params, \
        const CodeTree<type>& tree);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/optimize_synth.cc"
// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <algorithm>
#include <assert.h>

// line removed for fpoptimizer.cc: #include "optimize.hh"

using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace
{
    /* Synthesize the given grammatic parameter into the codetree */
    template<typename Value_t>
    CodeTree<Value_t> SynthesizeParam(
        const ParamSpec& parampair,
        MatchInfo<Value_t> & info,
        bool inner = true)
    {
        switch( parampair.first )
        {
            case NumConstant:
              { const ParamSpec_NumConstant<Value_t>& param = *(const ParamSpec_NumConstant<Value_t>*) parampair.second;
                return CodeTreeImmed( param.constvalue );
              }
            case ParamHolder:
              { const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                return info.GetParamHolderValue( param.index );
              }
            case SubFunction:
              { const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                CodeTree<Value_t> tree;
                tree.SetOpcode( param.data.subfunc_opcode );
                for(unsigned a=0; a < param.data.param_count; ++a)
                {
                    CodeTree<Value_t> nparam =
                        SynthesizeParam( ParamSpec_Extract<Value_t>(param.data.param_list, a),
                                         info, true );
                    tree.AddParamMove(nparam);
                }
                if(param.data.restholder_index != 0)
                {
                    std::vector<CodeTree<Value_t> > trees
                        ( info.GetRestHolderValues( param.data.restholder_index ) );
                    tree.AddParamsMove(trees);
                    // ^note: this fails if the same restholder is synth'd twice
                    if(tree.GetParamCount() == 1)
                    {
                        /* Convert cMul <1> into <1> when <1> only contains one operand.
                         * This is redundant code; it is also done in ConstantFolding(),
                         * but it might be better for performance to do it here, too.
                         */
                        assert(tree.GetOpcode() == cAdd || tree.GetOpcode() == cMul
                            || tree.GetOpcode() == cMin || tree.GetOpcode() == cMax
                            || tree.GetOpcode() == cAnd || tree.GetOpcode() == cOr
                            || tree.GetOpcode() == cAbsAnd || tree.GetOpcode() == cAbsOr);
                        tree.Become(tree.GetParam(0));
                    }
                    else if(tree.GetParamCount() == 0)
                    {
                        switch(tree.GetOpcode())
                        {
                            case cAdd: case cOr:
                                tree = CodeTreeImmed(Value_t(0));
                                break;
                            case cMul: case cAnd:
                                tree = CodeTreeImmed(Value_t(1));
                            default: break;
                        }
                    }
                }
                if(inner) tree.Rehash();
                return tree;
              }
        }
        return CodeTree<Value_t> ();
    }
}

namespace FPoptimizer_Optimize
{
    template<typename Value_t>
    void SynthesizeRule(
        const Rule& rule,
        CodeTree<Value_t>& tree,
        MatchInfo<Value_t>& info)
    {
        switch(rule.ruletype)
        {
            case ProduceNewTree:
            {
                tree.Become(
                    SynthesizeParam( ParamSpec_Extract<Value_t>(rule.repl_param_list, 0),
                                     info, false ) );
                break;
            }
            case ReplaceParams:
            default:
            {
                /* Delete the matched parameters from the source tree */
                std::vector<unsigned> list = info.GetMatchedParamIndexes();
                std::sort(list.begin(), list.end());
                for(size_t a=list.size(); a-->0; )
                    tree.DelParam( list[a] );

                /* Synthesize the replacement params */
                for(unsigned a=0; a < rule.repl_param_count; ++a)
                {
                    CodeTree<Value_t> nparam
                     = SynthesizeParam( ParamSpec_Extract<Value_t>(rule.repl_param_list, a),
                                        info, true );
                    tree.AddParamMove(nparam);
                }
                break;
            }
        }
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_Optimize
{
#define FP_INSTANTIATE(type) \
    template \
    void SynthesizeRule( \
        const Rule& rule, \
        CodeTree<type>& tree, \
        MatchInfo<type>& info);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/optimize_debug.cc"
// line removed for fpoptimizer.cc: #include "optimize.hh"
#ifdef DEBUG_SUBSTITUTIONS

// line removed for fpoptimizer.cc: #include "grammar.hh"
// line removed for fpoptimizer.cc: #include "opcodename.hh"

#include <sstream>
#include <cstring>

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_Grammar;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace FPoptimizer_Grammar
{
    template<typename Value_t>
    void DumpMatch(const Rule& rule,
                   const CodeTree<Value_t>& tree,
                   const MatchInfo<Value_t>& info,
                   bool DidMatch,
                   std::ostream& o)
    {
        DumpMatch(rule,tree,info,DidMatch?"Found match":"Found mismatch",o);
    }

    template<typename Value_t>
    void DumpMatch(const Rule& rule,
                   const CodeTree<Value_t>& tree,
                   const MatchInfo<Value_t>& info,
                   const char* whydump,
                   std::ostream& o)
    {
        static const char ParamHolderNames[][2] = {"%","&","x","y","z","a","b","c"};

        o << whydump
          << " (rule " << (&rule - grammar_rules) << ")"
          << ":\n"
            "  Pattern    : ";
        { ParamSpec tmp;
          tmp.first = SubFunction;
          ParamSpec_SubFunction tmp2;
          tmp2.data = rule.match_tree;
          tmp.second = (const void*) &tmp2;
          DumpParam<Value_t>(tmp, o);
        }
        o << "\n"
            "  Replacement: ";
        DumpParams<Value_t>(rule.repl_param_list, rule.repl_param_count, o);
        o << "\n";

        o <<
            "  Tree       : ";
        DumpTree(tree, o);
        o << "\n";
        if(!std::strcmp(whydump,"Found match")) DumpHashes(tree, o);

        for(size_t a=0; a<info.paramholder_matches.size(); ++a)
        {
            if(!info.paramholder_matches[a].IsDefined()) continue;
            o << "           " << ParamHolderNames[a] << " = ";
            DumpTree(info.paramholder_matches[a], o);
            o << "\n";
        }

        for(size_t b=0; b<info.restholder_matches.size(); ++b)
        {
            if(!info.restholder_matches[b].first) continue;
            for(size_t a=0; a<info.restholder_matches[b].second.size(); ++a)
            {
                o << "         <" << b << "> = ";
                DumpTree(info.restholder_matches[b].second[a], o);
                o << std::endl;
            }
        }
        o << std::flush;
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_Grammar
{
#define FP_INSTANTIATE(type) \
    template void DumpMatch(const Rule& rule, \
                   const CodeTree<type>& tree, \
                   const MatchInfo<type>& info, \
                   bool DidMatch, \
                   std::ostream& o);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/hash.cc"
#include <list>
#include <algorithm>

// line removed for fpoptimizer.cc: #include "constantfolding.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"
// line removed for fpoptimizer.cc: #include "../lib/crc32.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    template<typename Value_t>
    bool MarkIncompletes(FPoptimizer_CodeTree::CodeTree<Value_t>& tree)
    {
        if(tree.Is_Incompletely_Hashed())
            return true;

        bool needs_rehash = false;
        for(size_t a=0; a<tree.GetParamCount(); ++a)
            needs_rehash |= MarkIncompletes(tree.GetParam(a));
        if(needs_rehash)
            tree.Mark_Incompletely_Hashed();
        return needs_rehash;
    }

    template<typename Value_t>
    void FixIncompletes(FPoptimizer_CodeTree::CodeTree<Value_t>& tree)
    {
        if(tree.Is_Incompletely_Hashed())
        {
            for(size_t a=0; a<tree.GetParamCount(); ++a)
                FixIncompletes(tree.GetParam(a));
            tree.Rehash();
        }
    }
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void CodeTree<Value_t>::Sort()
    {
        data->Sort();
    }

    template<typename Value_t>
    void CodeTree<Value_t>::Rehash(bool constantfolding)
    {
        if(constantfolding)
            ConstantFolding(*this); // also runs Sort()
        else
            Sort();

        data->Recalculate_Hash_NoRecursion();
    }

    template<typename Value_t>
    struct ImmedHashGenerator
    {
        static void MakeHash(
            FUNCTIONPARSERTYPES::fphash_t& NewHash,
            const Value_t& Value)
        {
            /* TODO: For non-POD types, convert the value
             * into a base-62 string (or something) and hash that.
             */
            NewHash.hash1 = 0; // Try to ensure immeds gets always sorted first
          #if 0
            long double value = Value;
            fphash_value_t key = crc32::calc((const unsigned char*)&value, sizeof(value));
            key ^= (key << 24);
          #elif 0
            union
            {
                struct
                {
                    unsigned char filler1[16];
                    Value_t       v;
                    unsigned char filler2[16];
                } buf2;
                struct
                {
                    unsigned char filler3[sizeof(Value_t)+16-sizeof(fphash_value_t)];
                    fphash_value_t key;
                } buf1;
            } data;
            memset(&data, 0, sizeof(data));
            data.buf2.v = Value;
            fphash_value_t key = data.buf1.key;
          #else
            int exponent;
            Value_t fraction = std::frexp(Value, &exponent);
            fphash_value_t key = (unsigned(exponent+0x8000) & 0xFFFF);
            if(fraction < 0)
                { fraction = -fraction; key = key^0xFFFF; }
            else
                key += 0x10000;
            fraction -= Value_t(0.5);
            key <<= 39; // covers bits 39..55 now
            key |= fphash_value_t((fraction+fraction) * Value_t(1u<<31)) << 8;
            // fraction covers bits 8..39 now
          #endif
            /* Key = 56-bit unsigned integer value
             *       that is directly proportional
             *       to the floating point value.
             */
            NewHash.hash1 |= key;
            //crc32_t crc = crc32::calc((const unsigned char*)&Value, sizeof(Value));
            fphash_value_t crc = (key >> 10) | (key << (64-10));
            NewHash.hash2 += ((~fphash_value_t(crc)) * 3) ^ 1234567;
        }
    };

#ifdef FP_SUPPORT_COMPLEX_NUMBERS
    template<typename T>
    struct ImmedHashGenerator< std::complex<T> >
    {
        static void MakeHash(
            FUNCTIONPARSERTYPES::fphash_t& NewHash,
            const std::complex<T>& Value)
        {
            ImmedHashGenerator<T>::MakeHash(NewHash, Value.real());
            FUNCTIONPARSERTYPES::fphash_t temp;
            ImmedHashGenerator<T>::MakeHash(temp, Value.imag());
            NewHash.hash1 ^= temp.hash2;
            NewHash.hash2 ^= temp.hash1;
        }
    };
#endif

#ifdef FP_SUPPORT_LONG_INT_TYPE
    template<>
    struct ImmedHashGenerator<long>
    {
        static void MakeHash(
            FUNCTIONPARSERTYPES::fphash_t& NewHash,
            long Value)
        {
            fphash_value_t key = Value;
            /* Key = 56-bit unsigned integer value
             *       that is directly proportional
             *       to the floating point value.
             */
            NewHash.hash1 |= key;
            //crc32_t crc = crc32::calc((const unsigned char*)&Value, sizeof(Value));
            fphash_value_t crc = (key >> 10) | (key << (64-10));
            NewHash.hash2 += ((~fphash_value_t(crc)) * 3) ^ 1234567;
        }
    };
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
    template<>
    struct ImmedHashGenerator<GmpInt>
    {
        static void MakeHash(
            FUNCTIONPARSERTYPES::fphash_t& NewHash,
            const GmpInt& Value)
        {
            fphash_value_t key = Value.toInt();
            /* Key = 56-bit unsigned integer value
             *       that is directly proportional
             *       to the floating point value.
             */
            NewHash.hash1 |= key;
            //crc32_t crc = crc32::calc((const unsigned char*)&Value, sizeof(Value));
            fphash_value_t crc = (key >> 10) | (key << (64-10));
            NewHash.hash2 += ((~fphash_value_t(crc)) * 3) ^ 1234567;
        }
    };
#endif

    template<typename Value_t>
    void CodeTreeData<Value_t>::Recalculate_Hash_NoRecursion()
    {
        /* Hash structure:
         *     hash1: sorting key (8 bytes, 64 bits)
         *              byte 1: opcode
         *     hash2: unique value
         */
        fphash_t NewHash ( fphash_value_t(Opcode) << 56,
                           Opcode * FPHASH_CONST(0x1131462E270012B) );
        Depth = 1;
        switch(Opcode)
        {
            case cImmed:              // Value
            {
                ImmedHashGenerator<Value_t>::MakeHash(NewHash, Value);
                break; // no params
            }
            case VarBegin:            // Var_or_Funcno
            {
                NewHash.hash1 |= fphash_value_t(Var_or_Funcno) << 48;
                NewHash.hash2 += ((fphash_value_t(Var_or_Funcno)) * 11)
                                   ^ FPHASH_CONST(0x3A83A83A83A83A0);
                break; // no params
            }
            case cFCall: case cPCall: // Var_or_Funcno
            {
                NewHash.hash1 |= fphash_value_t(Var_or_Funcno) << 48;
                NewHash.hash2 += ((~fphash_value_t(Var_or_Funcno)) * 7) ^ 3456789;
                /* passthru */
            }
            default:
            {
                size_t MaxChildDepth = 0;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    if(Params[a].GetDepth() > MaxChildDepth)
                        MaxChildDepth = Params[a].GetDepth();

                    NewHash.hash1 += ((Params[a].GetHash().hash1*(a+1)) >> 12);
                    NewHash.hash2 += Params[a].GetHash().hash1;
                    NewHash.hash2 += (3)*FPHASH_CONST(0x9ABCD801357);
                    NewHash.hash2 *= FPHASH_CONST(0xECADB912345);
                    NewHash.hash2 += (~Params[a].GetHash().hash2) ^ 4567890;
                }
                Depth += MaxChildDepth;
            }
        }
        if(Hash != NewHash)
        {
            Hash = NewHash;
            OptimizedUsing = 0;
        }
    }

    template<typename Value_t>
    void CodeTree<Value_t>::FixIncompleteHashes()
    {
        MarkIncompletes(*this);
        FixIncompletes(*this);
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template void CodeTree<type>::Sort(); \
    template void CodeTree<type>::Rehash(bool); \
    template void CodeTree<type>::FixIncompleteHashes(); \
    template void CodeTreeData<type>::Recalculate_Hash_NoRecursion();
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/makebytecode.cc"
#include <cmath>
#include <list>
#include <cassert>

// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"
// line removed for fpoptimizer.cc: #include "consts.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"
// line removed for fpoptimizer.cc: #include "bytecodesynth.hh"

//#include "grammar.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    using namespace FPoptimizer_CodeTree;

    template<typename Value_t>
    bool AssembleSequence(
                  const CodeTree<Value_t>& tree, long count,
                  const FPoptimizer_ByteCode::SequenceOpCode<Value_t>& sequencing,
                  FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
                  size_t max_bytecode_grow_length);

    /*
    Trigonomic operations are expensive.
    If we need to synthesize a tan(x), we could
    first check if we can synthesize it by utilizing
    some expressions that are already in the stack.
    Such as: 1 / cot(x)
             sin(x) / cos(x)
             sin(x) * sec(x)
             (1/csc(x)) / cos(x)
             sec(x) / csc(x)
    This table lists the instructions for building
    each of these functions from components that
    might already exist in the stack.

    It is up to CSE to make it possible for this
    opportunity to arise as often as possible.
    */
    static const struct SinCosTanDataType
    {
        OPCODE whichopcode;
        OPCODE inverse_opcode;
        enum { nominator,
               denominator,
               inverse_nominator,
               inverse_denominator };
        OPCODE codes[4];
    } SinCosTanData[12] =
    {
        { cTan, cCot, { cSin,cCos, cCsc,cSec } },
        { cCot, cCot, { cCos,cSin, cSec,cCsc } },
        // tan(x) = 1/cot(x)
        //        = sin(x) / cos(x)
        //        = sin(x) * sec(x)
        //        = (1/csc(x)) / cos(x)
        //        = sec(x) / csc(x)

        { cCos, cSec, { cSin,cTan, cCsc,cCot } },
        { cSec, cCos, { cTan,cSin, cCot,cCsc } },
        // cos(x) = 1/sec(x) = sin(x) / tan(x)

        { cSin, cCsc, { cCos,cCot, cSec,cTan } },
        { cCsc, cSin, { cCot,cCos, cTan,cSec } },
        // sin(x) = 1/csc(x) = cos(x)/cot(x)

        { cTanh, cNop, { cSinh,cCosh, cNop,cNop } },
        { cSinh, cNop, { cTanh,cNop,  cNop,cCosh } },
        { cCosh, cNop, { cSinh,cTanh, cNop,cNop } },
        { cNop, cTanh, { cCosh,cSinh, cNop,cNop } },
        { cNop, cSinh, { cNop, cTanh, cCosh,cNop } },
        { cNop, cCosh, { cTanh,cSinh, cNop,cNop } }
    };
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void CodeTree<Value_t>::SynthesizeByteCode(
        std::vector<unsigned>& ByteCode,
        std::vector<Value_t>&   Immed,
        size_t& stacktop_max)
    {
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Making bytecode for:\n";
        DumpTreeWithIndent(*this);
    #endif
        while(RecreateInversionsAndNegations())
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "One change issued, produced:\n";
            DumpTreeWithIndent(*this);
        #endif
            FixIncompleteHashes();

            using namespace FPoptimizer_Optimize;
            using namespace FPoptimizer_Grammar;
            const void* g = (const void*)&grammar_optimize_recreate;
            while(ApplyGrammar(*(const Grammar*)g, *this))
                {   FixIncompleteHashes();
                }
        }
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Actually synthesizing, after recreating inv/neg:\n";
        DumpTreeWithIndent(*this);
    #endif

        FPoptimizer_ByteCode::ByteCodeSynth<Value_t> synth;

        /* Then synthesize the actual expression */
        SynthesizeByteCode(synth, false);
        /* The "false" parameters tells SynthesizeByteCode
         * that at the outermost synthesizing level, it does
         * not matter if leftover temps are left in the stack.
         */
        synth.Pull(ByteCode, Immed, stacktop_max);
    }

    template<typename Value_t>
    void CodeTree<Value_t>::SynthesizeByteCode(
        FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
        bool MustPopTemps) const
    {
        // If the synth can already locate our operand in the stack,
        // never mind synthesizing it again, just dup it.
        if(synth.FindAndDup(*this))
        {
            return;
        }

        for(size_t a=0; a<12; ++a)
        {
            const SinCosTanDataType& data = SinCosTanData[a];
            if(data.whichopcode != cNop)
            {
                if(GetOpcode() != data.whichopcode) continue;

                CodeTree invtree;
                invtree.SetParams(GetParams());
                invtree.SetOpcode( data.inverse_opcode );
                invtree.Rehash(false);
                if(synth.FindAndDup(invtree))
                {
                    synth.AddOperation(cInv,1,1);
                    synth.StackTopIs(*this);
                    return;
                }
            }
            else
            {
                // cNop indicates that there's no dedicated
                // opcode that indicates an inverted function.
                // For example, we have inv(cosh(x)).
                if(GetOpcode() != cInv) continue;
                if(GetParam(0).GetOpcode() != data.inverse_opcode) continue;
                if(synth.FindAndDup(GetParam(0)))
                {
                    synth.AddOperation(cInv,1,1);
                    synth.StackTopIs(*this);
                    return;
                }
            }

            // Check which trees we can find
            size_t   found[4];
            for(size_t b=0; b<4; ++b)
            {
                CodeTree tree;
                if(data.codes[b] == cNop)
                {
                    tree.SetOpcode(cInv);
                    CodeTree subtree;
                    subtree.SetParams(GetParams());
                    subtree.SetOpcode(data.codes[b ^ 2]);
                    subtree.Rehash(false);
                    tree.AddParamMove(subtree);
                }
                else
                {
                    tree.SetParams(GetParams());
                    tree.SetOpcode(data.codes[b]);
                }
                tree.Rehash(false);
                found[b] = synth.FindPos(tree);
            }

            if(found[ data.nominator ]   != ~size_t(0)
            && found[ data.denominator ] != ~size_t(0))
            {
                // tan(x) = sin(x) / cos(x)
                synth.DoDup( found[data.nominator] );
                synth.DoDup( found[data.denominator] );
                synth.AddOperation(cDiv,2,1);
                synth.StackTopIs(*this);
                return;
            }

            if(found[ data.nominator ]           != ~size_t(0)
            && found[ data.inverse_denominator ] != ~size_t(0))
            {
                // tan(x) = sin(x) * sec(x)
                synth.DoDup( found[data.nominator] );
                synth.DoDup( found[data.inverse_denominator] );
                synth.AddOperation(cMul,2,1);
                synth.StackTopIs(*this);
                return;
            }

            if(found[ data.inverse_nominator ]   != ~size_t(0)
            && found[ data.inverse_denominator ] != ~size_t(0))
            {
                // tan(x) = 1 / (csc(x) / sec(x)) = sec(x) / csc(x)
                synth.DoDup( found[data.inverse_nominator] );
                synth.DoDup( found[data.inverse_denominator] );
                synth.AddOperation(cRDiv,2,1);
                synth.StackTopIs(*this);
                return;
            }

            if(found[ data.inverse_nominator ]   != ~size_t(0)
            && found[ data.denominator ]         != ~size_t(0))
            {
                // tan(x) = 1 / csc(x) / cos(x) = 1 / (csc(x) * cos(x))
                synth.DoDup( found[data.inverse_nominator] );
                synth.DoDup( found[data.denominator] );
                synth.AddOperation(cMul,2,1);
                synth.AddOperation(cInv,1,1);
                synth.StackTopIs(*this);
                return;
            }
        }

        size_t n_subexpressions_synthesized = SynthCommonSubExpressions(synth);

        switch(GetOpcode())
        {
            case VarBegin:
                synth.PushVar(GetVar());
                break;
            case cImmed:
                synth.PushImmed(GetImmed());
                break;
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd:
            case cOr:
            case cAbsAnd:
            case cAbsOr:
            {
                if(GetOpcode() == cMul) // Special treatment for cMul sequences
                {
                    // If the paramlist contains an Immed, and that Immed
                    // fits in a long-integer, try to synthesize it
                    // as add-sequences instead.
                    bool did_muli = false;
                    for(size_t a=0; a<GetParamCount(); ++a)
                    {
                        if(GetParam(a).IsImmed() && isLongInteger(GetParam(a).GetImmed()))
                        {
                            long value = makeLongInteger(GetParam(a).GetImmed());

                            CodeTree tmp(*this, typename CodeTree::CloneTag());
                            tmp.DelParam(a);
                            tmp.Rehash();
                            if(AssembleSequence(
                                tmp, value,
                                FPoptimizer_ByteCode::SequenceOpcodes<Value_t>::AddSequence,
                                synth,
                                MAX_MULI_BYTECODE_LENGTH))
                            {
                                did_muli = true;
                                break;
                            }
                        }
                    }
                    if(did_muli)
                        break; // done
                }

                // If any of the params is currently a copy of
                // the stack topmost item, treat it first.
                int n_stacked = 0;
                std::vector<bool> done( GetParamCount() , false );
                CodeTree synthed_tree;
                synthed_tree.SetOpcode(GetOpcode());
                for(;;)
                {
                    bool found = false;
                    for(size_t a=0; a<GetParamCount(); ++a)
                    {
                        if(done[a]) continue;
                        if(synth.IsStackTop(GetParam(a)))
                        {
                            found = true;
                            done[a] = true;
                            GetParam(a).SynthesizeByteCode(synth);
                            synthed_tree.AddParam(GetParam(a));
                            if(++n_stacked > 1)
                            {
                                // Cumulate at the earliest opportunity.
                                synth.AddOperation(GetOpcode(), 2); // stack state: -2+1 = -1
                                synthed_tree.Rehash(false);
                                synth.StackTopIs(synthed_tree);
                                n_stacked = n_stacked - 2 + 1;
                            }
                        }
                    }
                    if(!found) break;
                }

                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    if(done[a]) continue;
                    GetParam(a).SynthesizeByteCode(synth);
                    synthed_tree.AddParam(GetParam(a));
                    if(++n_stacked > 1)
                    {
                        // Cumulate at the earliest opportunity.
                        synth.AddOperation(GetOpcode(), 2); // stack state: -2+1 = -1
                        synthed_tree.Rehash(false);
                        synth.StackTopIs(synthed_tree);
                        n_stacked = n_stacked - 2 + 1;
                    }
                }
                if(n_stacked == 0)
                {
                    // Uh, we got an empty cAdd/cMul/whatever...
                    // Synthesize a default value.
                    // This should never happen.
                    switch(GetOpcode())
                    {
                        case cAdd:
                        case cOr:
                        case cAbsOr:
                            synth.PushImmed(0);
                            break;
                        case cMul:
                        case cAnd:
                        case cAbsAnd:
                            synth.PushImmed(1);
                            break;
                        case cMin:
                        case cMax:
                            //synth.PushImmed(NaN);
                            synth.PushImmed(0);
                            break;
                        default:
                            break;
                    }
                    ++n_stacked;
                }
                assert(n_stacked == 1);
                break;
            }
            case cPow:
            {
                const CodeTree& p0 = GetParam(0);
                const CodeTree& p1 = GetParam(1);

                if(!p1.IsImmed()
                || !isLongInteger(p1.GetImmed())
                || !AssembleSequence( /* Optimize integer exponents */
                        p0, makeLongInteger(p1.GetImmed()),
                        FPoptimizer_ByteCode::SequenceOpcodes<Value_t>::MulSequence,
                        synth,
                        MAX_POWI_BYTECODE_LENGTH)
                  )
                {
                    p0.SynthesizeByteCode(synth);
                    p1.SynthesizeByteCode(synth);
                    synth.AddOperation(GetOpcode(), 2); // Create a vanilla cPow.
                }
                break;
            }
            case cIf:
            case cAbsIf:
            {
                // Assume that the parameter count is 3 as it should.
                typename FPoptimizer_ByteCode::ByteCodeSynth<Value_t>::IfData ifdata;

                GetParam(0).SynthesizeByteCode(synth); // expression

                synth.SynthIfStep1(ifdata, GetOpcode());

                GetParam(1).SynthesizeByteCode(synth); // true branch

                synth.SynthIfStep2(ifdata);

                GetParam(2).SynthesizeByteCode(synth); // false branch

                synth.SynthIfStep3(ifdata);
                break;
            }
            case cFCall:
            case cPCall:
            {
                // Assume that the parameter count is as it should.
                for(size_t a=0; a<GetParamCount(); ++a)
                    GetParam(a).SynthesizeByteCode(synth);
                synth.AddOperation(GetOpcode(), (unsigned) GetParamCount());
                synth.AddOperation(0x80000000u | GetFuncNo(), 0, 0);
                break;
            }
            default:
            {
                // Assume that the parameter count is as it should.
                for(size_t a=0; a<GetParamCount(); ++a)
                    GetParam(a).SynthesizeByteCode(synth);
                synth.AddOperation(GetOpcode(), (unsigned) GetParamCount());
                break;
            }
        }

        // Tell the synthesizer which tree was just produced in the stack
        synth.StackTopIs(*this);

        // If we added subexpressions, peel them off the stack now
        if(MustPopTemps && n_subexpressions_synthesized > 0)
        {
            size_t top = synth.GetStackTop();
            synth.DoPopNMov(top-1-n_subexpressions_synthesized, top-1);
        }
    }
}

namespace
{
    template<typename Value_t>
    bool AssembleSequence(
        const CodeTree<Value_t>& tree, long count,
        const FPoptimizer_ByteCode::SequenceOpCode<Value_t>& sequencing,
        FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
        size_t max_bytecode_grow_length)
    {
        if(count != 0)
        {
            FPoptimizer_ByteCode::ByteCodeSynth<Value_t> backup = synth;

            tree.SynthesizeByteCode(synth);

            // Ignore the size generated by subtree
            size_t bytecodesize_backup = synth.GetByteCodeSize();

            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);

            size_t bytecode_grow_amount = synth.GetByteCodeSize() - bytecodesize_backup;
            if(bytecode_grow_amount > max_bytecode_grow_length)
            {
                synth = backup;
                return false;
            }
            return true;
        }
        else
        {
            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);
            return true;
        }
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template void CodeTree<type>::SynthesizeByteCode( \
        std::vector<unsigned>& ByteCode, \
        std::vector<type>&   Immed, \
        size_t& stacktop_max);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/readbytecode.cc"
#include <cmath>
#include <cassert>

// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"
// line removed for fpoptimizer.cc: #include "opcodename.hh"
// line removed for fpoptimizer.cc: #include "grammar.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

// line removed for fpoptimizer.cc: #include "consts.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    using namespace FPoptimizer_CodeTree;

    #define FactorStack std::vector /* typedef*/

    const struct PowiMuliType
    {
        unsigned opcode_square;
        unsigned opcode_cumulate;
        unsigned opcode_invert;
        unsigned opcode_half;
        unsigned opcode_invhalf;
    } iseq_powi = {cSqr,cMul,cInv,cSqrt,cRSqrt},
      iseq_muli = {~unsigned(0), cAdd,cNeg, ~unsigned(0),~unsigned(0) };

    template<typename Value_t>
    Value_t ParsePowiMuli(
        const PowiMuliType& opcodes,
        const std::vector<unsigned>& ByteCode, size_t& IP,
        size_t limit,
        size_t factor_stack_base,
        FactorStack<Value_t>& stack)
    {
        Value_t result(1);
        while(IP < limit)
        {
            if(ByteCode[IP] == opcodes.opcode_square)
            {
                if(!isInteger(result)) break;
                result *= 2;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_invert)
            {
                result = -result;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_half)
            {
                if(result > Value_t(0) && isEvenInteger(result))
                    break;
                result *= Value_t(0.5);
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_invhalf)
            {
                if(result > Value_t(0) && isEvenInteger(result))
                    break;
                result *= Value_t(-0.5);
                ++IP;
                continue;
            }

            size_t dup_fetch_pos = IP;
            Value_t lhs(1);

            if(ByteCode[IP] == cFetch)
            {
                unsigned index = ByteCode[++IP];
                if(index < factor_stack_base
                || size_t(index-factor_stack_base) >= stack.size())
                {
                    // It wasn't a powi-fetch after all
                    IP = dup_fetch_pos;
                    break;
                }
                lhs = stack[index - factor_stack_base];
                // Note: ^This assumes that cFetch of recentmost
                //        is always converted into cDup.
                goto dup_or_fetch;
            }
            if(ByteCode[IP] == cDup)
            {
                lhs = result;
                goto dup_or_fetch;

            dup_or_fetch:
                stack.push_back(result);
                ++IP;
                Value_t subexponent = ParsePowiMuli
                    (opcodes,
                     ByteCode, IP, limit,
                     factor_stack_base, stack);
                if(IP >= limit || ByteCode[IP] != opcodes.opcode_cumulate)
                {
                    // It wasn't a powi-dup after all
                    IP = dup_fetch_pos;
                    break;
                }
                ++IP; // skip opcode_cumulate
                stack.pop_back();
                result += lhs*subexponent;
                continue;
            }
            break;
        }
        return result;
    }

    template<typename Value_t>
    Value_t ParsePowiSequence
        (const std::vector<unsigned>& ByteCode, size_t& IP,
         size_t limit, size_t factor_stack_base)
    {
        FactorStack<Value_t> stack;
        stack.push_back( Value_t(1) );
        return ParsePowiMuli(iseq_powi, ByteCode, IP, limit, factor_stack_base, stack);
    }

    template<typename Value_t>
    Value_t ParseMuliSequence
        (const std::vector<unsigned>& ByteCode, size_t& IP,
         size_t limit, size_t factor_stack_base)
    {
        FactorStack<Value_t> stack;
        stack.push_back( Value_t(1) );
        return ParsePowiMuli(iseq_muli, ByteCode, IP, limit, factor_stack_base, stack);
    }

    template<typename Value_t>
    class CodeTreeParserData
    {
    public:
        explicit CodeTreeParserData(bool k_powi)
            : stack(), clones(), keep_powi(k_powi) { }

        void Eat(size_t nparams, OPCODE opcode)
        {
            CodeTree<Value_t> newnode;
            newnode.SetOpcode(opcode);

            std::vector<CodeTree<Value_t> > params = Pop(nparams);
            newnode.SetParamsMove(params);

            if(!keep_powi)
            switch(opcode)
            {
                //        asinh: log(x + sqrt(x*x + 1))
                //cAsinh [x] -> cLog (cAdd x (cPow (cAdd (cPow x 2) 1) 0.5))
                // Note: ^ Replacement function refers to x twice

                //        acosh: log(x + sqrt(x*x - 1))
                //cAcosh [x] -> cLog (cAdd x (cPow (cAdd (cPow x 2) -1) 0.5))

                //        atanh: log( (1+x) / (1-x)) / 2
                //cAtanh [x] -> cMul (cLog (cMul (cAdd 1 x) (cPow (cAdd 1 (cMul -1 x)) -1))) 0.5

                //        asin: atan2(x, sqrt(1-x*x))
                //cAsin[x] -> cAtan2 [x (cPow [(cAdd 1 (cMul (cPow [x 2] -1)) 0.5])]

                //        acos: atan2(sqrt(1-x*x), x)
                //cAcos[x] -> cAtan2 [(cPow [(cAdd 1 (cMul (cPow [x 2] -1)) 0.5]) x]

                //     The hyperbolic functions themselves are:
                //        sinh: (exp(x)-exp(-x)) / 2  = exp(-x) * (exp(2*x)-1) / 2
                //cSinh [x] -> cMul 0.5 (cPow [CONSTANT_EI x]) (cAdd [-1 (cPow [CONSTANT_2E x])])

                //        cosh: (exp(x)+exp(-x)) / 2  = exp(-x) * (exp(2*x)+1) / 2
                //        cosh(-x) = cosh(x)
                //cCosh [x] -> cMul 0.5 (cPow [CONSTANT_EI x]) (cAdd [ 1 (cPow [CONSTANT_2E x])])

                //        tanh: sinh/cosh = (exp(2*x)-1) / (exp(2*x)+1)
                //cTanh [x] -> (cMul (cAdd {(cPow [CONSTANT_2E x]) -1}) (cPow [(cAdd {(cPow [CONSTANT_2E x]) 1}) -1]))
                case cTanh:
                {
                    CodeTree<Value_t> sinh, cosh;
                    sinh.SetOpcode(cSinh); sinh.AddParam(newnode.GetParam(0)); sinh.Rehash();
                    cosh.SetOpcode(cCosh); cosh.AddParamMove(newnode.GetParam(0)); cosh.Rehash();
                    CodeTree<Value_t> pow;
                    pow.SetOpcode(cPow);
                    pow.AddParamMove(cosh);
                    pow.AddParam(CodeTreeImmed(Value_t(-1)));
                    pow.Rehash();
                    newnode.SetOpcode(cMul);
                    newnode.SetParamMove(0, sinh);
                    newnode.AddParamMove(pow);
                    break;
                }

                //        tan: sin/cos
                //cTan [x] -> (cMul (cSin [x]) (cPow [(cCos [x]) -1]))
                case cTan:
                {
                    CodeTree<Value_t> sin, cos;
                    sin.SetOpcode(cSin); sin.AddParam(newnode.GetParam(0)); sin.Rehash();
                    cos.SetOpcode(cCos); cos.AddParamMove(newnode.GetParam(0)); cos.Rehash();
                    CodeTree<Value_t> pow;
                    pow.SetOpcode(cPow);
                    pow.AddParamMove(cos);
                    pow.AddParam(CodeTreeImmed(Value_t(-1)));
                    pow.Rehash();
                    newnode.SetOpcode(cMul);
                    newnode.SetParamMove(0, sin);
                    newnode.AddParamMove(pow);
                    break;
                }

                case cPow:
                {
                    const CodeTree<Value_t>& p0 = newnode.GetParam(0);
                    const CodeTree<Value_t>& p1 = newnode.GetParam(1);
                    if(p1.GetOpcode() == cAdd)
                    {
                        // convert x^(a + b) into x^a * x^b just so that
                        // some optimizations can be run on it.
                        // For instance, exp(log(x)*-61.1 + log(z)*-59.1)
                        // won't be changed into exp(log(x*z)*-61.1)*z^2
                        // unless we do this.
                        std::vector<CodeTree<Value_t> > mulgroup(p1.GetParamCount());
                        for(size_t a=0; a<p1.GetParamCount(); ++a)
                        {
                            CodeTree<Value_t> pow;
                            pow.SetOpcode(cPow);
                            pow.AddParam(p0);
                            pow.AddParam(p1.GetParam(a));
                            pow.Rehash();
                            mulgroup[a].swap(pow);
                        }
                        newnode.SetOpcode(cMul);
                        newnode.SetParamsMove(mulgroup);
                    }
                    break;
                }

                // Should we change sin(x) into cos(pi/2-x)
                //               or cos(x) into sin(pi/2-x)?
                //                        note: cos(x-pi/2) = cos(pi/2-x) = sin(x)
                //                        note: sin(x-pi/2) = -sin(pi/2-x) = -cos(x)
                default: break;
            }

            newnode.Rehash(!keep_powi);
        /*
            using namespace FPoptimizer_Grammar;
            bool recurse = false;
            while(ApplyGrammar(pack.glist[0], newnode, recurse)) // intermediate
            { //std::cout << "Rerunning 1\n";
                newnode.FixIncompleteHashes();
                recurse = true;
            }
        */
            FindClone(newnode, false);
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "POP " << nparams << ", " << FP_GetOpcodeName(opcode)
                      << "->" << FP_GetOpcodeName(newnode.GetOpcode())
                      << ": PUSH ";
            DumpTree(newnode);
            std::cout <<std::endl;
            DumpHashes(newnode);
        #endif
            stack.push_back(newnode);
        }

        void EatFunc(size_t nparams, OPCODE opcode, unsigned funcno)
        {
            CodeTree<Value_t> newnode = CodeTreeFuncOp<Value_t> (opcode, funcno);
            std::vector<CodeTree<Value_t> > params = Pop(nparams);
            newnode.SetParamsMove(params);
            newnode.Rehash(false);
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "POP " << nparams << ", PUSH ";
            DumpTree(newnode);
            std::cout << std::endl;
            DumpHashes(newnode);
        #endif
            FindClone(newnode);
            stack.push_back(newnode);
        }

        void AddConst(const Value_t& value)
        {
            CodeTree<Value_t> newnode = CodeTreeImmed(value);
            FindClone(newnode);
            Push(newnode);
        }

        void AddVar(unsigned varno)
        {
            CodeTree<Value_t> newnode = CodeTreeVar<Value_t>(varno);
            FindClone(newnode);
            Push(newnode);
        }

        void SwapLastTwoInStack()
        {
            stack[stack.size()-1].swap( stack[stack.size()-2] );
        }

        void Dup()
        {
            Fetch(stack.size()-1);
        }

        void Fetch(size_t which)
        {
            Push(stack[which]);
        }

        template<typename T>
        void Push(T tree)
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "PUSH ";
            DumpTree(tree);
            std::cout << std::endl;
            DumpHashes(tree);
        #endif
            stack.push_back(tree);
        }

        void PopNMov(size_t target, size_t source)
        {
            stack[target] = stack[source];
            stack.resize(target+1);
        }

        CodeTree<Value_t> PullResult()
        {
            clones.clear();
            CodeTree<Value_t> result(stack.back());
            stack.resize(stack.size()-1);
            return result;
        }
        std::vector<CodeTree<Value_t> > Pop(size_t n_pop)
        {
            std::vector<CodeTree<Value_t> > result(n_pop);
            for(unsigned n=0; n<n_pop; ++n)
                result[n].swap(stack[stack.size()-n_pop+n]);
        #ifdef DEBUG_SUBSTITUTIONS
            for(size_t n=n_pop; n-- > 0; )
            {
                std::cout << "POP ";
                DumpTree(result[n]);
                std::cout << std::endl;
                DumpHashes(result[n]);
            }
        #endif
            stack.resize(stack.size()-n_pop);
            return result;
        }

        size_t GetStackTop() const { return stack.size(); }
    private:
        void FindClone(CodeTree<Value_t> & /*tree*/, bool /*recurse*/ = true)
        {
            // Disabled: Causes problems in optimization when
            // the same subtree is included in logical and non-logical
            // contexts: optimizations applied to the logical one will
            // mess up the non-logical one.
            return;
            /*
            std::multimap<fphash_t, CodeTree>::const_iterator
                i = clones.lower_bound(tree.GetHash());
            for(; i != clones.end() && i->first == tree.GetHash(); ++i)
            {
                if(i->second.IsIdenticalTo(tree))
                    tree.Become(i->second);
            }
            if(recurse)
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                    FindClone(tree.GetParam(a));
            clones.insert(std::make_pair(tree.GetHash(), tree));
            */
        }
    private:
        std::vector<CodeTree<Value_t> > stack;
        std::multimap<fphash_t, CodeTree<Value_t> > clones;

        bool keep_powi;

    private:
        CodeTreeParserData(const CodeTreeParserData&);
        CodeTreeParserData& operator=(const CodeTreeParserData&);
    };

    template<typename Value_t>
    struct IfInfo
    {
        CodeTree<Value_t> condition;
        CodeTree<Value_t> thenbranch;
        size_t endif_location;

        IfInfo(): condition(), thenbranch(), endif_location() { }
    };
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void CodeTree<Value_t>::GenerateFrom(
        const typename FunctionParserBase<Value_t>::Data& fpdata,
        bool keep_powi)
    {
        std::vector<CodeTree<Value_t> > var_trees;
        var_trees.reserve(fpdata.mVariablesAmount);
        for(unsigned n=0; n<fpdata.mVariablesAmount; ++n)
        {
            var_trees.push_back( CodeTreeVar<Value_t> (n+VarBegin) );
        }
        GenerateFrom(fpdata,var_trees,keep_powi);
    }

    template<typename Value_t>
    void CodeTree<Value_t>::GenerateFrom(
        const typename FunctionParserBase<Value_t>::Data& fpdata,
        const std::vector<CodeTree>& var_trees,
        bool keep_powi)
    {
        const std::vector<unsigned>& ByteCode = fpdata.mByteCode;
        const std::vector<Value_t>&  Immed    = fpdata.mImmed;

        /*for(unsigned i=0; i<ByteCode.size(); ++i)
            fprintf(stderr, "by[%u/%u]=%u\n", i, (unsigned)ByteCode.size(), (unsigned) ByteCode[i]);*/

    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "ENTERS GenerateFrom()\n";
    #endif
        CodeTreeParserData<Value_t> sim(keep_powi);
        std::vector<IfInfo<Value_t> > if_stack;

        for(size_t IP=0, DP=0; ; ++IP)
        {
        after_powi:
            while(!if_stack.empty() &&
              (   // Normal If termination rule:
                  if_stack.back().endif_location == IP
                  // This rule matches when cJumps are threaded:
               || (IP < ByteCode.size() && ByteCode[IP] == cJump
                   && if_stack.back().thenbranch.IsDefined())
              ))
            {
                // The "else" of an "if" ends here
                CodeTree elsebranch = sim.PullResult();
                sim.Push(if_stack.back().condition);
                sim.Push(if_stack.back().thenbranch);
                sim.Push(elsebranch);
                sim.Eat(3, cIf);
                if_stack.pop_back();
            }
            if(IP >= ByteCode.size()) break;

            unsigned opcode = ByteCode[IP];
            if((opcode == cSqr || opcode == cDup
             || (opcode == cInv && !IsIntType<Value_t>::result)
             || opcode == cNeg
             || opcode == cSqrt || opcode == cRSqrt
             || opcode == cFetch))
            {
                // Parse a powi sequence
                size_t was_ip = IP;
                Value_t exponent = ParsePowiSequence<Value_t>(
                    ByteCode, IP, if_stack.empty() ? ByteCode.size() : if_stack.back().endif_location,
                    sim.GetStackTop()-1);
                if(exponent != Value_t(1.0))
                {
                    //std::cout << "Found exponent at " << was_ip << ": " << exponent << "\n";
                    sim.AddConst(exponent);
                    sim.Eat(2, cPow);
                    goto after_powi;
                }
                if(opcode == cDup
                || opcode == cFetch
                || opcode == cNeg)
                {
                    Value_t factor = ParseMuliSequence<Value_t>(
                        ByteCode, IP, if_stack.empty() ? ByteCode.size() : if_stack.back().endif_location,
                        sim.GetStackTop()-1);
                    if(factor != Value_t(1.0))
                    {
                        //std::cout << "Found factor at " << was_ip << ": " << factor << "\n";
                        sim.AddConst(factor);
                        sim.Eat(2, cMul);
                        goto after_powi;
                    }
                }
                IP = was_ip;
            }
            if(OPCODE(opcode) >= VarBegin)
            {
                unsigned index = opcode-VarBegin;
                /*std::fprintf(stderr, "IP=%u, opcode=%u, VarBegin=%u, ind=%u\n",
                    (unsigned)IP,(unsigned)opcode,(unsigned)VarBegin, (unsigned)index);*/
                sim.Push(var_trees[index]);
            }
            else
            {
                switch( OPCODE(opcode) )
                {
                    // Specials
                    case cIf:
                    case cAbsIf:
                    {
                        if_stack.resize(if_stack.size() + 1);
                        CodeTree res( sim.PullResult() );
                        if_stack.back().condition.swap( res );
                        if_stack.back().endif_location = ByteCode.size();
                        IP += 2; // dp,sp for elsebranch are irrelevant.
                        continue;
                    }
                    case cJump:
                    {
                        CodeTree res( sim.PullResult() );
                        if_stack.back().thenbranch.swap( res );
                        if_stack.back().endif_location = ByteCode[IP+1]+1;
                        IP += 2;
                        continue;
                    }
                    case cImmed:
                        sim.AddConst(Immed[DP++]);
                        break;
                    case cDup:
                        sim.Dup();
                        break;
                    case cNop:
                        break;
                    case cFCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        assert(funcno < fpdata.mFuncPtrs.size());
                        unsigned params = fpdata.mFuncPtrs[funcno].mParams;
                        sim.EatFunc(params, OPCODE(opcode), funcno);
                        break;
                    }
                    case cPCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        assert(funcno < fpdata.mFuncParsers.size());
                        const FunctionParserBase<Value_t>& p =
                            *fpdata.mFuncParsers[funcno].mParserPtr;
                        unsigned params = fpdata.mFuncParsers[funcno].mParams;

                        /* Inline the procedure call */
                        /* Works because cPCalls can never recurse */
                        std::vector<CodeTree> paramlist = sim.Pop(params);
                        CodeTree pcall_tree;
                        pcall_tree.GenerateFrom(*p.mData, paramlist);
                        sim.Push(pcall_tree);
                        break;
                    }
                    // Unary operators requiring special attention
                    case cInv:  // already handled by powi_opt
                        //sim.Eat(1, cInv);
                        //break;
                        sim.AddConst(1);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cDiv);
                        break;
                    case cNeg: // already handled by powi_opt
                        sim.Eat(1, cNeg);
                        break;
                        sim.AddConst(0);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cSub);
                        break;
                    case cSqr: // already handled by powi_opt
                        //sim.Eat(1, cSqr);
                        //break;
                        sim.AddConst(2);
                        sim.Eat(2, cPow);
                        break;
                    // Unary functions requiring special attention
                    case cSqrt: // already handled by powi_opt
                        sim.AddConst( Value_t(0.5) );
                        sim.Eat(2, cPow);
                        break;
                    case cRSqrt: // already handled by powi_opt
                        sim.AddConst( Value_t(-0.5) );
                        sim.Eat(2, cPow);
                        break;
                    case cCbrt:
                        sim.AddConst(Value_t(1) / Value_t(3));
                        sim.Eat(2, cPow);
                        break;
                    case cDeg:
                        sim.AddConst(fp_const_rad_to_deg<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cRad:
                        sim.AddConst(fp_const_deg_to_rad<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cExp:
                        if(keep_powi) goto default_function_handling;
                        sim.AddConst(fp_const_e<Value_t>());
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cPow);
                        break;
                    case cExp2: // from fpoptimizer
                        if(keep_powi) goto default_function_handling;
                        sim.AddConst(2.0);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cPow);
                        break;
                    case cCot:
                        sim.Eat(1, cTan);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cCsc:
                        sim.Eat(1, cSin);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cSec:
                        sim.Eat(1, cCos);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cInt: // int(x) = floor(x + 0.5)
                    #ifndef __x86_64
                        if(keep_powi) { sim.Eat(1, cInt); break; }
                    #endif
                        sim.AddConst( Value_t(0.5) );
                        sim.Eat(2, cAdd);
                        sim.Eat(1, cFloor);
                        break;
                    case cLog10:
                        sim.Eat(1, cLog);
                        sim.AddConst(fp_const_log10inv<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cLog2:
                        sim.Eat(1, cLog);
                        sim.AddConst(fp_const_log2inv<Value_t>());
                        sim.Eat(2, cMul);
                        break;
                    case cLog2by: // x y     -> log(x)*CONSTANT_L2I*y
                        sim.SwapLastTwoInStack();   // y x
                        sim.Eat(1, cLog);           // y log(x)
                        sim.AddConst(fp_const_log2inv<Value_t>()); // y log(x) CONSTANT_L2I
                        sim.Eat(3, cMul);           // y*log(x)*CONSTANT_L2I
                        break;
                    case cHypot: // x y -> sqrt(x*x + y*y)
                        sim.AddConst(2);
                        sim.Eat(2, cPow); // x y^2
                        sim.SwapLastTwoInStack();
                        sim.AddConst(2);
                        sim.Eat(2, cPow); // y^2 x^2
                        sim.Eat(2, cAdd); // y^2 + x^2
                        sim.AddConst( Value_t(0.5) );
                        sim.Eat(2, cPow); // (y^2 + x^2)^0.5
                        break;
                    case cSinCos:
                        sim.Dup();
                        sim.Eat(1, cSin);
                        sim.SwapLastTwoInStack();
                        sim.Eat(1, cCos);
                        break;
                    case cSinhCosh:
                        sim.Dup();
                        sim.Eat(1, cSinh);
                        sim.SwapLastTwoInStack();
                        sim.Eat(1, cCosh);
                        break;
                    //case cLog:
                    //    sim.Eat(1, cLog2);
                    //    sim.AddConst(fp_const_log2<Value_t>());
                    //    sim.Eat(2, cMul);
                    //    break;
                    // Binary operators requiring special attention
                    case cRSub: // from fpoptimizer
                        sim.SwapLastTwoInStack();
                        // Passthru to cSub
                    case cSub:
                        if(keep_powi) { sim.Eat(2, cSub); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cMul); // -x is x*-1
                        sim.Eat(2, cAdd); // Minus is negative adding
                        break;
                    case cRDiv: // from fpoptimizer
                        sim.SwapLastTwoInStack();
                        // Passthru to cDiv
                    case cDiv:
                        if(keep_powi || IsIntType<Value_t>::result)
                        {
                            sim.Eat(2, cDiv);
                            break;
                        }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow); // 1/x is x^-1
                        sim.Eat(2, cMul); // Divide is inverse multiply
                        break;
                    // Binary operators not requiring special attention
                    case cAdd: case cMul:
                    case cMod: case cPow:
                    case cEqual: case cLess: case cGreater:
                    case cNEqual: case cLessOrEq: case cGreaterOrEq:
                    case cAnd: case cOr:
                    case cAbsAnd: case cAbsOr:
                        sim.Eat(2, OPCODE(opcode));
                        break;
                    // Unary operators not requiring special attention
                    case cNot:
                    case cNotNot: // from fpoptimizer
                    case cAbsNot:
                    case cAbsNotNot:
                        sim.Eat(1, OPCODE(opcode));
                        break;
                    // Special opcodes generated by fpoptimizer itself
                    case cFetch:
                        sim.Fetch(ByteCode[++IP]);
                        break;
                    case cPopNMov:
                    {
                        unsigned stackOffs_target = ByteCode[++IP];
                        unsigned stackOffs_source = ByteCode[++IP];
                        sim.PopNMov(stackOffs_target, stackOffs_source);
                        break;
                    }

                    default:
                    default_function_handling:;
                        unsigned funcno = opcode-cAbs;
                        assert(funcno < FUNC_AMOUNT);
                        const FuncDefinition& func = Functions[funcno];
                        sim.Eat(func.params, OPCODE(opcode));
                        break;
                }
            }
        }
        Become(sim.PullResult());
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Produced tree:\n";
        DumpTreeWithIndent(*this);
    #endif
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template \
    void CodeTree<type>::GenerateFrom( \
        const FunctionParserBase<type>::Data& fpdata, \
        bool keep_powi);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/constantfolding.cc"
#include <algorithm>

// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"
// line removed for fpoptimizer.cc: #include "consts.hh"

#include <assert.h>

// line removed for fpoptimizer.cc: #include "rangeestimation.hh"
// line removed for fpoptimizer.cc: #include "constantfolding.hh"


// line removed for fpoptimizer.cc: #include "logic_boolgroups.hh"
/* ^For ConstantFolding_AndLogic()
 *      ConstantFolding_OrLogic()
 *      ConstantFolding_MulLogicItems()
 *      ConstantFolding_AddLogicItems()
 */

// line removed for fpoptimizer.cc: #include "logic_collections.hh"
/* ^For ConstantFolding_MulGrouping()
 *      ConstantFolding_AddGrouping()
 */

// line removed for fpoptimizer.cc: #include "logic_ifoperations.hh"
/* ^For ConstantFolding_IfOperations()
 */

// line removed for fpoptimizer.cc: #include "logic_powoperations.hh"
/* ^For ConstantFolding_PowOperations()
 */

// line removed for fpoptimizer.cc: #include "logic_comparisons.hh"
/* ^For ConstantFolding_Comparison()
 */

#define FP_MUL_COMBINE_EXPONENTS

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /*************************************/
    /* ADOPTING SAME-TYPE CHILDREN       */
    /*************************************/

    template<typename Value_t>
    static void AdoptChildrenWithSameOpcode(CodeTree<Value_t>& tree)
    {
        /* If the list contains another list of the same kind, assimilate it */
      #ifdef DEBUG_SUBSTITUTIONS
        bool assimilated = false;
      #endif
        for(size_t a=tree.GetParamCount(); a-- > 0; )
            if(tree.GetParam(a).GetOpcode() == tree.GetOpcode())
            {
              #ifdef DEBUG_SUBSTITUTIONS
                if(!assimilated)
                {
                    std::cout << "Before assimilation: "; DumpTree(tree);
                    std::cout << "\n";
                    assimilated = true;
                }
              #endif
                // Assimilate its children and remove it
                tree.AddParamsMove(tree.GetParam(a).GetUniqueRef().GetParams(), a);
            }
      #ifdef DEBUG_SUBSTITUTIONS
        if(assimilated)
        {
            std::cout << "After assimilation:   "; DumpTree(tree);
            std::cout << "\n";
        }
      #endif
    }
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void ConstantFolding(CodeTree<Value_t>& tree)
    {
        tree.Sort(); // Otherwise "0 <= acos(x)" does not get properly optimized
    #ifdef DEBUG_SUBSTITUTIONS
        void* stackptr=0;
        std::cout << "[" << (&stackptr) << "]Runs ConstantFolding for: ";
        DumpTree(tree);
        std::cout << "\n";
        DumpHashes(tree); std::cout << std::flush;
    #endif
        if(false)
        {
    redo:;
            tree.Sort();
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "[" << (&stackptr) << "]Re-runs ConstantFolding: ";
            DumpTree(tree);
            std::cout << "\n";
            DumpHashes(tree);
        #endif
        }

        // Insert here any hardcoded constant-folding optimizations
        // that you want to be done whenever a new subtree is generated.
        /* Not recursive. */

        if(tree.GetOpcode() != cImmed)
        {
            range<Value_t> p = CalculateResultBoundaries(tree);
            if(p.min.known && p.max.known && p.min.val == p.max.val)
            {
                // Replace us with this immed
                tree.ReplaceWithImmed(p.min.val);
                goto do_return;
            }
        }

        if(false)
        {
            ReplaceTreeWithOne:  tree.ReplaceWithImmed(Value_t(1)); goto do_return;
            ReplaceTreeWithZero: tree.ReplaceWithImmed(Value_t(0)); goto do_return;
            ReplaceTreeWithParam0:
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before replace: ";
                std::cout << std::hex
                          << '[' << tree.GetHash().hash1
                          << ',' << tree.GetHash().hash2
                          << ']' << std::dec;
                DumpTree(tree);
                std::cout << "\n";
              #endif
                tree.Become(tree.GetParam(0));
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After replace: ";
                std::cout << std::hex
                          << '[' << tree.GetHash().hash1
                          << ',' << tree.GetHash().hash2
                          << ']' << std::dec;
                DumpTree(tree);
                std::cout << "\n";
              #endif
                goto redo;
        }

        /* Constant folding */
        switch(tree.GetOpcode())
        {
            case cImmed:
                break; // nothing to do
            case VarBegin:
                break; // nothing to do

            case cAnd:
            case cAbsAnd:
            {
                AdoptChildrenWithSameOpcode(tree);
                bool has_nonlogical_values = false;
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                {
                    if(!IsLogicalValue(tree.GetParam(a))) has_nonlogical_values = true;
                    switch(GetLogicalValue(tree.GetParam(a), tree.GetOpcode()==cAbsAnd))
                    {
                        case IsNever: goto ReplaceTreeWithZero;
                        case IsAlways: tree.DelParam(a); break; // x & y & 1 = x & y;  x & 1 = !!x
                        case Unknown: default: ;
                    }
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithOne;
                    case 1: tree.SetOpcode(tree.GetOpcode()==cAnd ? cNotNot : cAbsNotNot); goto redo; // Replace self with the single operand
                    default:
                        if(tree.GetOpcode()==cAnd || !has_nonlogical_values)
                            if(ConstantFolding_AndLogic(tree)) goto redo;
                }
                break;
            }
            case cOr:
            case cAbsOr:
            {
                AdoptChildrenWithSameOpcode(tree);
                bool has_nonlogical_values = false;
                for(size_t a=tree.GetParamCount(); a-- > 0; )
                {
                    if(!IsLogicalValue(tree.GetParam(a))) has_nonlogical_values = true;
                    switch(GetLogicalValue(tree.GetParam(a), tree.GetOpcode()==cAbsOr))
                    {
                        case IsAlways: goto ReplaceTreeWithOne;
                        case IsNever: tree.DelParam(a); break;
                        case Unknown: default: ;
                    }
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: tree.SetOpcode(tree.GetOpcode()==cOr ? cNotNot : cAbsNotNot); goto redo; // Replace self with the single operand
                    default:
                        if(tree.GetOpcode()==cOr || !has_nonlogical_values)
                            if(ConstantFolding_OrLogic(tree)) goto redo;
                }
                break;
            }
            case cNot:
            case cAbsNot:
            {
                unsigned opposite = 0;
                switch(tree.GetParam(0).GetOpcode())
                {
                    case cEqual:       opposite = cNEqual; break;
                    case cNEqual:      opposite = cEqual; break;
                    case cLess:        opposite = cGreaterOrEq; break;
                    case cGreater:     opposite = cLessOrEq; break;
                    case cLessOrEq:    opposite = cGreater; break;
                    case cGreaterOrEq: opposite = cLess; break;
                    case cNotNot:      opposite = cNot; break;
                    case cNot:         opposite = cNotNot; break;
                    case cAbsNot:      opposite = cAbsNotNot; break;
                    case cAbsNotNot:   opposite = cAbsNot; break;
                    default: break;
                }
                if(opposite)
                {
                    tree.SetOpcode(OPCODE(opposite));
                    tree.SetParamsMove(tree.GetParam(0).GetUniqueRef().GetParams());
                    goto redo;
                }

                // If the sub-expression evaluates to approx. zero, yield one.
                // If the sub-expression evaluates to approx. nonzero, yield zero.
                switch(GetLogicalValue(tree.GetParam(0), tree.GetOpcode()==cAbsNot))
                {
                    case IsAlways: goto ReplaceTreeWithZero;
                    case IsNever: goto ReplaceTreeWithOne;
                    case Unknown: default: ;
                }
                if(tree.GetOpcode() == cNot && GetPositivityInfo(tree.GetParam(0)) == IsAlways)
                    tree.SetOpcode(cAbsNot);

                if(tree.GetParam(0).GetOpcode() == cIf
                || tree.GetParam(0).GetOpcode() == cAbsIf)
                {
                    CodeTree<Value_t> iftree = tree.GetParam(0);
                    const CodeTree<Value_t>& ifp1 = iftree.GetParam(1);
                    const CodeTree<Value_t>& ifp2 = iftree.GetParam(2);
                    if(ifp1.GetOpcode() == cNot
                    || ifp1.GetOpcode() == cAbsNot)
                    {
                        // cNot [(cIf [x (cNot[y]) z])] -> cIf [x (cNotNot[y]) (cNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree<Value_t> p1;
                        p1.SetOpcode(ifp1.GetOpcode()==cNot ? cNotNot : cAbsNotNot);
                        p1.AddParam(ifp1.GetParam(0));
                        p1.Rehash();
                        tree.AddParamMove(p1);
                        CodeTree<Value_t> p2;
                        p2.SetOpcode(tree.GetOpcode());
                        p2.AddParam(ifp2);
                        p2.Rehash();
                        tree.AddParamMove(p2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                    if(ifp2.GetOpcode() == cNot
                    || ifp2.GetOpcode() == cAbsNot)
                    {
                        // cNot [(cIf [x y (cNot[z])])] -> cIf [x (cNot[y]) (cNotNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree<Value_t> p1;
                        p1.SetOpcode(tree.GetOpcode());
                        p1.AddParam(ifp1);
                        p1.Rehash();
                        tree.AddParamMove(p1);
                        CodeTree<Value_t> p2;
                        p2.SetOpcode(ifp2.GetOpcode()==cNot ? cNotNot : cAbsNotNot);
                        p2.AddParam(ifp2.GetParam(0));
                        p2.Rehash();
                        tree.AddParamMove(p2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                }
                break;
            }
            case cNotNot:
            case cAbsNotNot:
            {
                // The function of cNotNot is to protect a logical value from
                // changing. If the parameter is already a logical value,
                // then the cNotNot opcode is redundant.
                if(IsLogicalValue(tree.GetParam(0)))
                    goto ReplaceTreeWithParam0;

                // If the sub-expression evaluates to approx. zero, yield zero.
                // If the sub-expression evaluates to approx. nonzero, yield one.
                switch(GetLogicalValue(tree.GetParam(0), tree.GetOpcode()==cAbsNotNot))
                {
                    case IsNever: goto ReplaceTreeWithZero;
                    case IsAlways: goto ReplaceTreeWithOne;
                    case Unknown: default: ;
                }
                if(tree.GetOpcode() == cNotNot && GetPositivityInfo(tree.GetParam(0)) == IsAlways)
                    tree.SetOpcode(cAbsNotNot);

                if(tree.GetParam(0).GetOpcode() == cIf
                || tree.GetParam(0).GetOpcode() == cAbsIf)
                {
                    CodeTree<Value_t> iftree = tree.GetParam(0);
                    const CodeTree<Value_t>& ifp1 = iftree.GetParam(1);
                    const CodeTree<Value_t>& ifp2 = iftree.GetParam(2);
                    if(ifp1.GetOpcode() == cNot
                    || ifp1.GetOpcode() == cAbsNot)
                    {
                        // cNotNot [(cIf [x (cNot[y]) z])] -> cIf [x (cNot[y]) (cNotNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        tree.AddParam(ifp1);
                        CodeTree<Value_t> p2;
                        p2.SetOpcode(tree.GetOpcode());
                        p2.AddParam(ifp2);
                        p2.Rehash();
                        tree.AddParamMove(p2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                    if(ifp2.GetOpcode() == cNot
                    || ifp2.GetOpcode() == cAbsNot)
                    {
                        // cNotNot [(cIf [x y (cNot[z])])] -> cIf [x (cNotNot[y]) (cNot[z])]
                        tree.SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree<Value_t> p1;
                        p1.SetOpcode(tree.GetOpcode());
                        p1.AddParam(ifp1);
                        p1.Rehash();
                        tree.AddParamMove(p1);
                        tree.AddParam(ifp2);
                        tree.SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                }
                break;
            }
            case cIf:
            case cAbsIf:
            {
                if(ConstantFolding_IfOperations(tree))
                    goto redo;
                break;
            }
            case cMul:
            {
            NowWeAreMulGroup: ;
                AdoptChildrenWithSameOpcode(tree);
                // If one sub-expression evalutes to exact zero, yield zero.
                Value_t immed_product = Value_t(1);
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    if(!tree.GetParam(a).IsImmed()) continue;
                    // ^ Only check constant values
                    Value_t immed = tree.GetParam(a).GetImmed();
                    if(immed == Value_t(0) ) goto ReplaceTreeWithZero;
                    immed_product *= immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1 || (n_immeds == 1 && fp_equal(immed_product, Value_t(1))))
                    needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "cMul: Will add new immed " << immed_product << "\n";
                #endif
                    for(size_t a=tree.GetParamCount(); a-->0; )
                        if(tree.GetParam(a).IsImmed())
                        {
                        #ifdef DEBUG_SUBSTITUTIONS
                            std::cout << " - For that, deleting immed " << tree.GetParam(a).GetImmed();
                            std::cout << "\n";
                        #endif
                            tree.DelParam(a);
                        }
                    if(!fp_equal(immed_product, Value_t(1)))
                        tree.AddParam( CodeTreeImmed<Value_t> (immed_product) );
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithOne;
                    case 1: goto ReplaceTreeWithParam0; // Replace self with the single operand
                    default:
                        if(ConstantFolding_MulGrouping(tree)) goto redo;
                        if(ConstantFolding_MulLogicItems(tree)) goto redo;
                }
                break;
            }
            case cAdd:
            {
                AdoptChildrenWithSameOpcode(tree);
                Value_t immed_sum = 0.0;
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    if(!tree.GetParam(a).IsImmed()) continue;
                    // ^ Only check constant values
                    Value_t immed = tree.GetParam(a).GetImmed();
                    immed_sum += immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1 || (n_immeds == 1 && immed_sum == Value_t(0)))
                    needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "cAdd: Will add new immed " << immed_sum << "\n";
                    std::cout << "In: "; DumpTree(tree);
                    std::cout << "\n";
                #endif
                    for(size_t a=tree.GetParamCount(); a-->0; )
                        if(tree.GetParam(a).IsImmed())
                        {
                        #ifdef DEBUG_SUBSTITUTIONS
                            std::cout << " - For that, deleting immed " << tree.GetParam(a).GetImmed();
                            std::cout << "\n";
                        #endif
                            tree.DelParam(a);
                        }
                    if(!(immed_sum == Value_t(0.0)))
                        tree.AddParam( CodeTreeImmed<Value_t> (immed_sum) );
                }
                switch(tree.GetParamCount())
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: goto ReplaceTreeWithParam0; // Replace self with the single operand
                    default:
                        if(ConstantFolding_AddGrouping(tree)) goto redo;
                        if(ConstantFolding_AddLogicItems(tree)) goto redo;
                }
                break;
            }
            case cMin:
            {
                AdoptChildrenWithSameOpcode(tree);
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the larger one.
                 *
                 * Algorithm: 1. figure out the smallest maximum of all operands.
                 *            2. eliminate all operands where their minimum is
                 *               larger than the selected maximum.
                 */
                size_t preserve=0;
                range<Value_t> smallest_maximum;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    while(a+1 < tree.GetParamCount() && tree.GetParam(a).IsIdenticalTo(tree.GetParam(a+1)))
                        tree.DelParam(a+1);
                    range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                    if(p.max.known && (!smallest_maximum.max.known || (p.max.val) < smallest_maximum.max.val))
                    {
                        smallest_maximum.max.val = p.max.val;
                        smallest_maximum.max.known = true;
                        preserve=a;
                }   }
                if(smallest_maximum.max.known)
                    for(size_t a=tree.GetParamCount(); a-- > 0; )
                    {
                        range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                        if(p.min.known && a != preserve && p.min.val >= smallest_maximum.max.val)
                            tree.DelParam(a);
                    }
                //fprintf(stderr, "Remains: %u\n", (unsigned)tree.GetParamCount());
                if(tree.GetParamCount() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }
            case cMax:
            {
                AdoptChildrenWithSameOpcode(tree);
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the smaller one.
                 *
                 * Algorithm: 1. figure out the biggest minimum of all operands.
                 *            2. eliminate all operands where their maximum is
                 *               smaller than the selected minimum.
                 */
                size_t preserve=0;
                range<Value_t> biggest_minimum;
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                {
                    while(a+1 < tree.GetParamCount() && tree.GetParam(a).IsIdenticalTo(tree.GetParam(a+1)))
                        tree.DelParam(a+1);
                    range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                    if(p.min.known && (!biggest_minimum.min.known || p.min.val > biggest_minimum.min.val))
                    {
                        biggest_minimum.min.val = p.min.val;
                        biggest_minimum.min.known = true;
                        preserve=a;
                }   }
                if(biggest_minimum.min.known)
                {
                    //fprintf(stderr, "Removing all where max < %g\n", biggest_minimum.min.val);
                    for(size_t a=tree.GetParamCount(); a-- > 0; )
                    {
                        range<Value_t> p = CalculateResultBoundaries( tree.GetParam(a) );
                        if(p.max.known && a != preserve && (p.max.val) < biggest_minimum.min.val)
                        {
                            //fprintf(stderr, "Removing %g\n", p.max.val);
                            tree.DelParam(a);
                        }
                    }
                }
                //fprintf(stderr, "Remains: %u\n", (unsigned)tree.GetParamCount());
                if(tree.GetParamCount() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }

            case cEqual:
            case cNEqual:
            case cLess:
            case cGreater:
            case cLessOrEq:
            case cGreaterOrEq:
                if(ConstantFolding_Comparison(tree)) goto redo;
                break;

            case cAbs:
            {
                /* If we know the operand is always positive, cAbs is redundant.
                 * If we know the operand is always negative, use actual negation.
                 */
                range<Value_t> p0 = CalculateResultBoundaries( tree.GetParam(0) );
                if(p0.min.known && p0.min.val >= Value_t(0.0))
                    goto ReplaceTreeWithParam0;
                if(p0.max.known && p0.max.val <= fp_const_negativezero<Value_t>())
                {
                    /* abs(negative) = negative*-1 */
                    tree.SetOpcode(cMul);
                    tree.AddParam( CodeTreeImmed(Value_t(1)) );
                    /* The caller of ConstantFolding() will do Sort() and Rehash() next.
                     * Thus, no need to do it here. */
                    /* We were changed into a cMul group. Do cMul folding. */
                    goto NowWeAreMulGroup;
                }
                /* If the operand is a cMul group, find elements
                 * that are always positive and always negative,
                 * and move them out, e.g. abs(p*n*x*y) = p*(-n)*abs(x*y)
                 */
                if(tree.GetParam(0).GetOpcode() == cMul)
                {
                    const CodeTree<Value_t>& p = tree.GetParam(0);
                    std::vector<CodeTree<Value_t> > pos_set;
                    std::vector<CodeTree<Value_t> > neg_set;
                    for(size_t a=0; a<p.GetParamCount(); ++a)
                    {
                        p0 = CalculateResultBoundaries( p.GetParam(a) );
                        if(p0.min.known && p0.min.val >= Value_t(0.0))
                            { pos_set.push_back(p.GetParam(a)); }
                        if(p0.max.known && p0.max.val <= fp_const_negativezero<Value_t>())
                            { neg_set.push_back(p.GetParam(a)); }
                    }
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "Abs: mul group has " << pos_set.size()
                              << " pos, " << neg_set.size() << "neg\n";
                #endif
                    if(!pos_set.empty() || !neg_set.empty())
                    {
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-Before: ";
                        DumpTree(tree);
                        std::cout << "\n" << std::flush;
                        DumpHashes(tree, std::cout);
                #endif
                        CodeTree<Value_t> pclone;
                        pclone.SetOpcode(cMul);
                        for(size_t a=0; a<p.GetParamCount(); ++a)
                        {
                            p0 = CalculateResultBoundaries( p.GetParam(a) );
                            if((p0.min.known && p0.min.val >= Value_t(0.0))
                            || (p0.max.known && p0.max.val <= fp_const_negativezero<Value_t>()))
                                {/*pclone.DelParam(a);*/}
                            else
                                pclone.AddParam( p.GetParam(a) );
                            /* Here, p*n*x*y -> x*y.
                             * p is saved in pos_set[]
                             * n is saved in neg_set[]
                             */
                        }
                        pclone.Rehash();
                        CodeTree<Value_t> abs_mul;
                        abs_mul.SetOpcode(cAbs);
                        abs_mul.AddParamMove(pclone);
                        abs_mul.Rehash();
                        CodeTree<Value_t> mulgroup;
                        mulgroup.SetOpcode(cMul);
                        mulgroup.AddParamMove(abs_mul); // cAbs[whatever remains in p]
                        mulgroup.AddParamsMove(pos_set);
                        /* Now:
                         * mulgroup  = p * Abs(x*y)
                         */
                        if(!neg_set.empty())
                        {
                            if(neg_set.size() % 2)
                                mulgroup.AddParam( CodeTreeImmed(Value_t(-1)) );
                            mulgroup.AddParamsMove(neg_set);
                            /* Now:
                             * mulgroup = p * n * -1 * Abs(x*y)
                             */
                        }
                        tree.Become(mulgroup);
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-After: ";
                        DumpTree(tree, std::cout);
                        std::cout << "\n" << std::flush;
                        DumpHashes(tree, std::cout);
                #endif
                        /* We were changed into a cMul group. Do cMul folding. */
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            }

            #define HANDLE_UNARY_CONST_FUNC(funcname) \
                if(tree.GetParam(0).IsImmed()) \
                    { tree.ReplaceWithImmed( funcname(tree.GetParam(0).GetImmed()) ); \
                      goto do_return; }

            case cLog:
                HANDLE_UNARY_CONST_FUNC(fp_log);
                if(tree.GetParam(0).GetOpcode() == cPow)
                {
                    CodeTree<Value_t> pow = tree.GetParam(0);
                    if(GetPositivityInfo(pow.GetParam(0)) == IsAlways)  // log(posi ^ y) = y*log(posi)
                    {
                        pow.CopyOnWrite();
                        pow.SetOpcode(cLog);
                        tree.SetOpcode(cMul);
                        tree.AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        tree.SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                    if(GetEvennessInfo(pow.GetParam(1)) == IsAlways) // log(x ^ even) = even*log(abs(x))
                    {
                        pow.CopyOnWrite();
                        CodeTree<Value_t> abs;
                        abs.SetOpcode(cAbs);
                        abs.AddParamMove(pow.GetParam(0));
                        abs.Rehash();
                        pow.SetOpcode(cLog);
                        tree.SetOpcode(cMul);
                        pow.SetParamMove(0, abs);
                        tree.AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        tree.SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                }
                else if(tree.GetParam(0).GetOpcode() == cAbs)
                {
                    // log(abs(x^y)) = y*log(abs(x))
                    CodeTree<Value_t> pow = tree.GetParam(0).GetParam(0);
                    if(pow.GetOpcode() == cPow)
                    {
                        pow.CopyOnWrite();
                        CodeTree<Value_t> abs;
                        abs.SetOpcode(cAbs);
                        abs.AddParamMove(pow.GetParam(0));
                        abs.Rehash();
                        pow.SetOpcode(cLog);
                        tree.SetOpcode(cMul);
                        pow.SetParamMove(0, abs);
                        tree.AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        tree.SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            case cAcosh: HANDLE_UNARY_CONST_FUNC(fp_acosh); break;
            case cAsinh: HANDLE_UNARY_CONST_FUNC(fp_asinh); break;
            case cAtanh: HANDLE_UNARY_CONST_FUNC(fp_atanh); break;
            case cAcos: HANDLE_UNARY_CONST_FUNC(fp_acos); break;
            case cAsin: HANDLE_UNARY_CONST_FUNC(fp_asin); break;
            case cAtan: HANDLE_UNARY_CONST_FUNC(fp_atan); break;
            case cCosh: HANDLE_UNARY_CONST_FUNC(fp_cosh); break;
            case cSinh: HANDLE_UNARY_CONST_FUNC(fp_sinh); break;
            case cTanh: HANDLE_UNARY_CONST_FUNC(fp_tanh); break;
            case cSin: HANDLE_UNARY_CONST_FUNC(fp_sin); break;
            case cCos: HANDLE_UNARY_CONST_FUNC(fp_cos); break;
            case cTan: HANDLE_UNARY_CONST_FUNC(fp_tan); break;
            case cCeil:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_ceil); break;
            case cTrunc:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_trunc); break;
            case cFloor:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_floor); break;
            case cInt:
                if(GetIntegerInfo(tree.GetParam(0)) == IsAlways) goto ReplaceTreeWithParam0;
                HANDLE_UNARY_CONST_FUNC(fp_int); break;
            case cCbrt: HANDLE_UNARY_CONST_FUNC(fp_cbrt); break; // converted into cPow x 0.33333
            case cSqrt: HANDLE_UNARY_CONST_FUNC(fp_sqrt); break; // converted into cPow x 0.5
            case cExp: HANDLE_UNARY_CONST_FUNC(fp_exp); break; // convered into cPow CONSTANT_E x
            case cLog2: HANDLE_UNARY_CONST_FUNC(fp_log2); break;
            case cLog10: HANDLE_UNARY_CONST_FUNC(fp_log10); break;

            case cLog2by:
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_log2(tree.GetParam(0).GetImmed()) * tree.GetParam(1).GetImmed() );
                      goto do_return; }
                break;

            case cArg: HANDLE_UNARY_CONST_FUNC(fp_arg); break;
            case cConj: HANDLE_UNARY_CONST_FUNC(fp_conj); break;
            case cImag: HANDLE_UNARY_CONST_FUNC(fp_imag); break;
            case cReal: HANDLE_UNARY_CONST_FUNC(fp_real); break;
            case cPolar:
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_polar(tree.GetParam(0).GetImmed(), tree.GetParam(1).GetImmed() ) );
                      goto do_return; }
                break;

            case cMod: /* Can more be done than this? */
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_mod(tree.GetParam(0).GetImmed(), tree.GetParam(1).GetImmed()) );
                      goto do_return; }
                break;

            case cAtan2:
            {
                /* Range based optimizations for (y,x):
                 * If y is +0 and x <= -0, +pi is returned
                 * If y is -0 and x <= -0, -pi is returned (assumed never happening)
                 * If y is +0 and x >= +0, +0 is returned
                 * If y is -0 and x >= +0, -0 is returned  (assumed never happening)
                 * If x is +-0 and y < 0, -pi/2 is returned
                 * If x is +-0 and y > 0, +pi/2 is returned
                 * Otherwise, perform constant folding when available
                 * If we know x <> 0, convert into atan(y / x)
                 *   TODO: Figure out whether the above step is wise
                 *         It allows e.g. atan2(6*x, 3*y) -> atan(2*x/y)
                 *         when we know y != 0
                 */
                range<Value_t> p0 = CalculateResultBoundaries( tree.GetParam(0) );
                range<Value_t> p1 = CalculateResultBoundaries( tree.GetParam(1) );
                if(tree.GetParam(0).IsImmed()
                && fp_equal(tree.GetParam(0).GetImmed(), Value_t(0)))   // y == 0
                {
                    if(p1.max.known && (p1.max.val) < Value_t(0))    // y == 0 && x < 0
                        { tree.ReplaceWithImmed( fp_const_pi<Value_t>() ); goto do_return; }
                    if(p1.min.known && p1.min.val >= Value_t(0.0))  // y == 0 && x >= 0.0
                        { tree.ReplaceWithImmed( Value_t(0) ); goto do_return; }
                }
                if(tree.GetParam(1).IsImmed()
                && fp_equal(tree.GetParam(1).GetImmed(), Value_t(0)))   // x == 0
                {
                    if(p0.max.known && (p0.max.val) < Value_t(0))   // y < 0 && x == 0
                        { tree.ReplaceWithImmed( -fp_const_pihalf<Value_t>() ); goto do_return; }
                    if(p0.min.known && p0.min.val > Value_t(0))     // y > 0 && x == 0
                        { tree.ReplaceWithImmed(  fp_const_pihalf<Value_t>() ); goto do_return; }
                }
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( fp_atan2(tree.GetParam(0).GetImmed(),
                                             tree.GetParam(1).GetImmed()) );
                      goto do_return; }
                if((p1.min.known && p1.min.val > Value_t(0))            // p1 != 0.0
                || (p1.max.known && (p1.max.val) < fp_const_negativezero<Value_t>())) // become atan(p0 / p1)
                {
                    CodeTree<Value_t> pow_tree;
                    pow_tree.SetOpcode(cPow);
                    pow_tree.AddParamMove(tree.GetParam(1));
                    pow_tree.AddParam(CodeTreeImmed(Value_t(-1)));
                    pow_tree.Rehash();
                    CodeTree<Value_t> div_tree;
                    div_tree.SetOpcode(cMul);
                    div_tree.AddParamMove(tree.GetParam(0));
                    div_tree.AddParamMove(pow_tree);
                    div_tree.Rehash();
                    tree.SetOpcode(cAtan);
                    tree.SetParamMove(0, div_tree);
                    tree.DelParam(1);
                }
                break;
            }

            case cPow:
            {
                if(ConstantFolding_PowOperations(tree)) goto redo;
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.cc and thus
             * they will never occur in the calling context for the
             * most of the parsing context. They may however occur
             * at the late phase, so we deal with them.
             */
            case cDiv: // converted into cPow y -1
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed()
                && tree.GetParam(1).GetImmed() != Value_t(0.0))
                    { tree.ReplaceWithImmed( tree.GetParam(0).GetImmed() / tree.GetParam(1).GetImmed() );
                      goto do_return; }
                break;
            case cInv: // converted into cPow y -1
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(0).GetImmed() != Value_t(0.0))
                    { tree.ReplaceWithImmed( Value_t(1) / tree.GetParam(0).GetImmed() );
                      goto do_return; }
                // Note: Could use (mulgroup)^immed optimization from cPow
                break;
            case cSub: // converted into cMul y -1
                if(tree.GetParam(0).IsImmed()
                && tree.GetParam(1).IsImmed())
                    { tree.ReplaceWithImmed( tree.GetParam(0).GetImmed() - tree.GetParam(1).GetImmed() );
                      goto do_return; }
                break;
            case cNeg: // converted into cMul x -1
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( -tree.GetParam(0).GetImmed() );
                      goto do_return; }
                break;
            case cRad: // converted into cMul x CONSTANT_RD
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( RadiansToDegrees( tree.GetParam(0).GetImmed() ) );
                      goto do_return; }
                break;
            case cDeg: // converted into cMul x CONSTANT_DR
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( DegreesToRadians( tree.GetParam(0).GetImmed() ) );
                      goto do_return; }
                break;
            case cSqr: // converted into cMul x x
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( tree.GetParam(0).GetImmed() * tree.GetParam(0).GetImmed() );
                      goto do_return; }
                break;
            case cExp2: // converted into cPow 2.0 x
                HANDLE_UNARY_CONST_FUNC(fp_exp2); break;
            case cRSqrt: // converted into cPow x -0.5
                if(tree.GetParam(0).IsImmed())
                    { tree.ReplaceWithImmed( Value_t(1) / fp_sqrt(tree.GetParam(0).GetImmed()) );
                      goto do_return; }
                break;
            case cCot: // converted into cMul (cPow (cTan x) -1)
                if(tree.GetParam(0).IsImmed())
                    { Value_t tmp = fp_tan(tree.GetParam(0).GetImmed());
                      if(fp_nequal(tmp, Value_t(0)))
                      { tree.ReplaceWithImmed( Value_t(1) / tmp );
                        goto do_return; } }
                break;
            case cSec: // converted into cMul (cPow (cCos x) -1)
                if(tree.GetParam(0).IsImmed())
                    { Value_t tmp = fp_cos(tree.GetParam(0).GetImmed());
                      if(fp_nequal(tmp, Value_t(0)))
                      { tree.ReplaceWithImmed( Value_t(1) / tmp );
                        goto do_return; } }
                break;
            case cCsc: // converted into cMul (cPow (cSin x) -1)
                if(tree.GetParam(0).IsImmed())
                    { Value_t tmp = fp_sin(tree.GetParam(0).GetImmed());
                      if(fp_nequal(tmp, Value_t(0)))
                      { tree.ReplaceWithImmed( Value_t(1) / tmp );
                        goto do_return; } }
                break;
            case cHypot: // converted into cSqrt(cAdd(cMul(x x), cMul(y y)))
                if(tree.GetParam(0).IsImmed() && tree.GetParam(1).IsImmed())
                {
                    tree.ReplaceWithImmed( fp_hypot(tree.GetParam(0).GetImmed(),
                                           tree.GetParam(1).GetImmed()) );
                    goto do_return;
                }
                break;

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
                break; /* Should never occur */

            /* Opcodes that we can't do anything about */
            case cPCall:
            case cFCall:
                break;
        }
    do_return:;
#ifdef DEBUG_SUBSTITUTIONS
        std::cout << "[" << (&stackptr) << "]Done ConstantFolding, result: ";
        DumpTree(tree);
        std::cout << "\n";
        DumpHashes(tree);
#endif
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template void ConstantFolding(CodeTree<type>& );
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/valuerange.cc"
// line removed for fpoptimizer.cc: #include "valuerange.hh"

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

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template struct range<type>; \
    template bool IsLogicalTrueValue(const range<type> &, bool); \
    template bool IsLogicalFalseValue(const range<type> &, bool);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/rangeestimation.cc"
// line removed for fpoptimizer.cc: #include "rangeestimation.hh"
// line removed for fpoptimizer.cc: #include "consts.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_CodeTree;

//#define DEBUG_SUBSTITUTIONS_extra_verbose

// Intel 11.1 has problems with std::isnan().  In libMesh we have some
// workarounds, but don't need to include all that infrastructure here.
// So, this little hack will do.
namespace
{
  template <typename T>
  inline
  int isnan_workaround(T t)
  {
    return (t != t);
  }
}

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
    range<Value_t> CalculateResultBoundaries_do(const CodeTree<Value_t>& tree)
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
                if(!res2.min.known)
                  res1.min.known = false;
                else if(res1.min.known && (res2.min.val) < res1.min.val)
                  res1.min.val = res2.min.val;
                else if (isnan_workaround(res2.min.val))
                  res1.min.val = res2.min.val;
                if(!res2.max.known)
                  res1.max.known = false;
                else if(res1.max.known && (res2.max.val) > res1.max.val)
                  res1.max.val = res2.max.val;
                else if (isnan_workaround(res2.max.val))
                  res1.max.val = res2.max.val;
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

                            // we've already determined the result to be positive, but these boundaries would result in a min = inf
                            if(p0.min.val == Value_t(0) || p1.min.val < Value_t(0))
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

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template range<type> CalculateResultBoundaries(const CodeTree<type> &); \
    template bool IsLogicalValue(const CodeTree<type> &); \
    template TriTruthValue GetIntegerInfo(const CodeTree<type> &);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/transformations.cc"
// line removed for fpoptimizer.cc: #include "codetree.hh"

#ifdef FP_SUPPORT_OPTIMIZER

// line removed for fpoptimizer.cc: #include "bytecodesynth.hh"
// line removed for fpoptimizer.cc: #include "rangeestimation.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh" // For DEBUG_SUBSTITUTIONS

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

//#define DEBUG_POWI

#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif

#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif

namespace FPoptimizer_ByteCode
{
    extern const unsigned char powi_table[256];
}
namespace
{
    using namespace FPoptimizer_CodeTree;

    template<typename Value_t>
    bool IsOptimizableUsingPowi(long immed, long penalty = 0)
    {
        FPoptimizer_ByteCode::ByteCodeSynth<Value_t> synth;
        synth.PushVar(VarBegin);
        // Ignore the size generated by subtree
        size_t bytecodesize_backup = synth.GetByteCodeSize();
        FPoptimizer_ByteCode::AssembleSequence(immed,
            FPoptimizer_ByteCode::SequenceOpcodes<Value_t>::MulSequence, synth);

        size_t bytecode_grow_amount = synth.GetByteCodeSize() - bytecodesize_backup;

        return bytecode_grow_amount < size_t(MAX_POWI_BYTECODE_LENGTH - penalty);
    }

    template<typename Value_t>
    void ChangeIntoRootChain(
        CodeTree<Value_t>& tree,
        bool inverted,
        long sqrt_count,
        long cbrt_count)
    {
        while(cbrt_count > 0)
        {
            CodeTree<Value_t> tmp;
            tmp.SetOpcode(cCbrt);
            tmp.AddParamMove(tree);
            tmp.Rehash();
            tree.swap(tmp);
            --cbrt_count;
        }
        while(sqrt_count > 0)
        {
            CodeTree<Value_t> tmp;
            tmp.SetOpcode(cSqrt);
            if(inverted)
            {
                tmp.SetOpcode(cRSqrt);
                inverted = false;
            }
            tmp.AddParamMove(tree);
            tmp.Rehash();
            tree.swap(tmp);
            --sqrt_count;
        }
        if(inverted)
        {
            CodeTree<Value_t> tmp;
            tmp.SetOpcode(cInv);
            tmp.AddParamMove(tree);
            tree.swap(tmp);
        }
    }

    template<typename Value_t>
    struct RootPowerTable
    {
        static const Value_t RootPowers[(1+4)*(1+3)];
    };
    template<typename Value_t>
    const Value_t RootPowerTable<Value_t>::RootPowers[(1+4)*(1+3)] =
    {
        // (sqrt^n(x)) // Workaround for gcc-4.6 bug by RHS
        Value_t(1),
        Value_t(1.L/2.L),
        Value_t(1.L/ (2.L*2.L)),
        Value_t(1.L/ (2.L*2.L*2.L)),
        Value_t(1.L/ (2.L*2.L*2.L*2.L)),
        // cbrt^1(sqrt^n(x))
        Value_t(1.L/ (3.L)),
        Value_t(1.L/ (3.L*2.L)),
        Value_t(1.L/ (3.L*2.L*2.L)),
        Value_t(1.L/ (3.L*2.L*2.L*2.L)),
        Value_t(1.L/ (3.L*2.L*2.L*2.L*2.L)),
        // cbrt^2(sqrt^n(x))
        Value_t(1.L/ (3.L*3.L)),
        Value_t(1.L/ (3.L*3.L*2.L)),
        Value_t(1.L/ (3.L*3.L*2.L*2.L)),
        Value_t(1.L/ (3.L*3.L*2.L*2.L*2.L)),
        Value_t(1.L/ (3.L*3.L*2.L*2.L*2.L*2.L)),
        // cbrt^3(sqrt^n(x))
        Value_t(1.L/ (3.L*3.L*3.L)),
        Value_t(1.L/ (3.L*3.L*3.L*2.L)),
        Value_t(1.L/ (3.L*3.L*3.L*2.L*2.L)),
        Value_t(1.L/ (3.L*3.L*3.L*2.L*2.L*2.L)),
        Value_t(1.L/ (3.L*3.L*3.L*2.L*2.L*2.L*2.L))
    };

    struct PowiResolver
    {
        /* Any exponentiation can be turned into one of these:
         *
         *   x^y  -> sqrt(x)^(y*2)         = x Sqrt       y*2  Pow
         *   x^y  -> cbrt(x)^(y*3)         = x Cbrt       y*3  Pow
         *   x^y  -> rsqrt(x)^(y*-2)       = x RSqrt     y*-2  Pow
         *   x^y  -> x^(y-1/2) * sqrt(x)   = x Sqrt   x y-0.5  Pow Mul
         *   x^y  -> x^(y-1/3) * cbrt(x)   = x Cbrt   x y-0.33 Pow Mul
         *   x^y  -> x^(y+1/2) * rsqrt(x)  = x Sqrt   x y+0.5  Pow Mul
         *   x^y  -> inv(x)^(-y)           = x Inv      -y     Pow
         *
         * These rules can be applied recursively.
         * The goal is to find the optimal chain of operations
         * that results in the least number of sqrt,cbrt operations;
         * an integer value of y, and that the integer is as close
         * to zero as possible.
         */
        static const unsigned MaxSep = 4;
        static const int      MaxOp  = 5;

        typedef int factor_t;
        typedef long cost_t;
        typedef long int_exponent_t;

        struct PowiResult
        {
            PowiResult() :
                n_int_sqrt(0),
                n_int_cbrt(0),
                sep_list(),
                resulting_exponent(0) { }

            int n_int_sqrt; // totals
            int n_int_cbrt; // totals
            int sep_list[MaxSep]; // action list. Each element is (n_sqrt + MaxOp * n_cbrt).
            int_exponent_t resulting_exponent;
        };

        template<typename Value_t>
        static PowiResult CreatePowiResult(Value_t exponent)
        {
            PowiResult result;

            factor_t best_factor = FindIntegerFactor(exponent);
            if(best_factor == 0)
            {
        #ifdef DEBUG_POWI
                printf("no factor found for %Lg\n", (long double)exponent);
        #endif
                return result; // Unoptimizable
            }

            result.resulting_exponent = MultiplyAndMakeLong(exponent, best_factor);
            cost_t best_cost =
                EvaluateFactorCost(best_factor, 0, 0, 0)
              + CalculatePowiFactorCost( result.resulting_exponent );
            int s_count = 0;
            int c_count = 0;
            int mul_count = 0;

        #ifdef DEBUG_POWI
            printf("orig = %Lg\n", (long double) exponent);
            printf("plain factor = %d, cost %ld\n", (int) best_factor, (long) best_cost);
        #endif

            for(unsigned n_s=0; n_s<MaxSep; ++n_s)
            {
                int best_selected_sep = 0;
                cost_t best_sep_cost  = best_cost;
                factor_t best_sep_factor = best_factor;
                for(int s=1; s<MaxOp*4; ++s)
                {
#ifdef CBRT_IS_SLOW
                    if(s >= MaxOp) break;
                    // When cbrt is implemented through exp and log,
                    // there is no advantage over exp(log()), so don't support it.
#endif
                    int n_sqrt = s%MaxOp;
                    int n_cbrt = s/MaxOp;
                    if(n_sqrt + n_cbrt > 4) continue;

                    Value_t changed_exponent = exponent;
                    changed_exponent -= RootPowerTable<Value_t>::RootPowers[s];

                    factor_t factor = FindIntegerFactor(changed_exponent);
                    if(factor != 0)
                    {
                        int_exponent_t int_exponent = MultiplyAndMakeLong(changed_exponent, factor);
                        cost_t cost =
                            EvaluateFactorCost(factor, s_count + n_sqrt, c_count + n_cbrt, mul_count + 1)
                          + CalculatePowiFactorCost(int_exponent);

        #ifdef DEBUG_POWI
                        printf("Candidate sep %u (%d*sqrt %d*cbrt)factor = %d, cost %ld (for %Lg to %ld)\n",
                            s, n_sqrt, n_cbrt, factor,
                            (long) cost,
                            (long double) changed_exponent,
                            (long) int_exponent);
        #endif
                        if(cost < best_sep_cost)
                        {
                            best_selected_sep = s;
                            best_sep_factor   = factor;
                            best_sep_cost     = cost;
                        }
                    }
                }
                if(!best_selected_sep) break;

        #ifdef DEBUG_POWI
                printf("CHOSEN sep %u (%d*sqrt %d*cbrt)factor = %d, cost %ld, exponent %Lg->%Lg\n",
                       best_selected_sep,
                       best_selected_sep % MaxOp,
                       best_selected_sep / MaxOp,
                       best_sep_factor, best_sep_cost,
                       (long double)(exponent),
                       (long double)(exponent-RootPowerTable<Value_t>::RootPowers[best_selected_sep]));
        #endif
                result.sep_list[n_s] = best_selected_sep;
                exponent -= RootPowerTable<Value_t>::RootPowers[best_selected_sep];
                s_count += best_selected_sep % MaxOp;
                c_count += best_selected_sep / MaxOp;
                best_cost   = best_sep_cost;
                best_factor = best_sep_factor;
                mul_count += 1;
            }

            result.resulting_exponent = MultiplyAndMakeLong(exponent, best_factor);
        #ifdef DEBUG_POWI
            printf("resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n",
                result.resulting_exponent,
                (long double) exponent,
                (long double) best_factor);
        #endif
            while(best_factor % 2 == 0)
            {
                ++result.n_int_sqrt;
                best_factor /= 2;
            }
            while(best_factor % 3 == 0)
            {
                ++result.n_int_cbrt;
                best_factor /= 3;
            }
            return result;
        }

    private:
        static cost_t CalculatePowiFactorCost(int_exponent_t int_exponent)
        {
            static std::map<int_exponent_t, cost_t> cache;
            if(int_exponent < 0)
            {
                cost_t cost = 22; // division cost
                return cost + CalculatePowiFactorCost(-int_exponent);
            }
            std::map<int_exponent_t,cost_t>::iterator i = cache.lower_bound(int_exponent);
            if(i != cache.end() && i->first == int_exponent)
                return i->second;
            std::pair<int_exponent_t, cost_t> result(int_exponent, 0.0);
            cost_t& cost = result.second;

            while(int_exponent > 1)
            {
                int factor = 0;
                if(int_exponent < 256)
                {
                    factor = FPoptimizer_ByteCode::powi_table[int_exponent];
                    if(factor & 128) factor &= 127; else factor = 0;
                    if(factor & 64) factor = -(factor&63) - 1;
                }
                if(factor)
                {
                    cost += CalculatePowiFactorCost(factor);
                    int_exponent /= factor;
                    continue;
                }
                if(!(int_exponent & 1))
                {
                    int_exponent /= 2;
                    cost += 6; // sqr
                }
                else
                {
                    cost += 7; // dup+mul
                    int_exponent -= 1;
                }
            }

            cache.insert(i, result);
            return cost;
        }

        template<typename Value_t>
        static int_exponent_t MultiplyAndMakeLong(const Value_t& value, factor_t factor)
        {
            return makeLongInteger( value * Value_t(factor) );
        }

        // Find the integer that "value" must be multiplied
        // with to produce an integer...
        // Consisting of factors 2 and 3 only.
        template<typename Value_t>
        static bool MakesInteger(const Value_t& value, factor_t factor)
        {
            /* Does value, multiplied by factor, result in an integer? */
            Value_t v = value * Value_t(factor);
            return isLongInteger(v);
            /*
            Value_t diff = fp_abs(v - fp_int(v));
            //printf("factor %d: v=%.20f, diff=%.20f\n", factor,v, diff);
            return diff < Value_t(1e-9l);
            */
        }

        template<typename Value_t>
        static factor_t FindIntegerFactor(const Value_t& value)
        {
            factor_t factor = (2*2*2*2);
#ifdef CBRT_IS_SLOW
            // When cbrt is implemented through exp and log,
            // there is no advantage over exp(log()), so don't support it.
#else
            factor *= (3*3*3);
#endif
            factor_t result = 0;
            if(MakesInteger(value, factor))
            {
                result = factor;
                while((factor % 2) == 0 && MakesInteger(value, factor/2))
                    result = factor /= 2;
                while((factor % 3) == 0 && MakesInteger(value, factor/3))
                    result = factor /= 3;
            }
#ifdef CBRT_IS_SLOW
            if(result == 0)
            {
                /* Note: Even if we allow one cbrt,
                 *        cbrt(cbrt(x)) still gets turned into
                 *        exp(log(x)*0.111111)
                 *        which gives an error when x < 0...
                 *        should we use a special system here?
                 *        i.e. exp(log(-5)*y)
                 *      =      -exp(log(5)*y)
                 *        except when y is an even integer,
                 *      when  = exp(log(5)*y)
                 * We use a custom fp_pow() function
                 * in order to handle these situations.
                 */
                if(MakesInteger(value, 3)) return 3; // single cbrt opcode
            }
#endif
            return result;
        }

        static int EvaluateFactorCost(int factor, int s, int c, int nmuls)
        {
            const int sqrt_cost = 6;
#ifdef CBRT_IS_SLOW
            const int cbrt_cost = 25;
#else
            const int cbrt_cost = 8;
#endif
            int result = s * sqrt_cost + c * cbrt_cost;
            while(factor % 2 == 0) { factor /= 2; result += sqrt_cost; }
            while(factor % 3 == 0) { factor /= 3; result += cbrt_cost; }
            result += nmuls;
            return result;
        }
    };
}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    bool CodeTree<Value_t>::RecreateInversionsAndNegations(bool prefer_base2)
    {
        bool changed = false;

        for(size_t a=0; a<GetParamCount(); ++a)
            if(GetParam(a).RecreateInversionsAndNegations(prefer_base2))
                changed = true;

        if(changed)
        {
        exit_changed:
            Mark_Incompletely_Hashed();
            return true;
        }

        switch(GetOpcode()) // Recreate inversions and negations
        {
            case cMul:
            {
                std::vector<CodeTree<Value_t> > div_params;
                CodeTree<Value_t> found_log2, found_log2by;

                if(true)
                {
                    /* This lengthy bit of code
                     * changes log2(x)^3 * 5
                     * to      log2by(x, 5^(1/3)) ^ 3
                     * which is better for runtime
                     * than    log2by(x,1)^3 * 5
                     */
                    bool found_log2_on_exponent = false;
                    Value_t log2_exponent = 0;
                    for(size_t a = GetParamCount(); a-- > 0; )
                    {
                        const CodeTree<Value_t>& powgroup = GetParam(a);
                        if(powgroup.GetOpcode() == cPow
                        && powgroup.GetParam(0).GetOpcode() == cLog2
                        && powgroup.GetParam(1).IsImmed())
                        {
                            // Found log2 on exponent
                            found_log2_on_exponent = true;
                            log2_exponent = powgroup.GetParam(1).GetImmed();
                            break;
                        }
                    }
                    if(found_log2_on_exponent)
                    {
                        Value_t immeds = 1.0;
                        for(size_t a = GetParamCount(); a-- > 0; )
                        {
                            const CodeTree<Value_t>& powgroup = GetParam(a);
                            if(powgroup.IsImmed())
                            {
                                immeds *= powgroup.GetImmed();
                                DelParam(a);
                            }
                        }
                        for(size_t a = GetParamCount(); a-- > 0; )
                        {
                            CodeTree<Value_t>& powgroup = GetParam(a);
                            if(powgroup.GetOpcode() == cPow
                            && powgroup.GetParam(0).GetOpcode() == cLog2
                            && powgroup.GetParam(1).IsImmed())
                            {
                                CodeTree<Value_t>& log2 = powgroup.GetParam(0);
                                log2.CopyOnWrite();
                                log2.SetOpcode(cLog2by);
                                log2.AddParam( CodeTreeImmed(
                                    fp_pow(immeds, Value_t(1) / log2_exponent) ) );
                                log2.Rehash();
                                break;
                            }
                        }
                    }
                }

                for(size_t a = GetParamCount(); a-- > 0; )
                {
                    const CodeTree<Value_t>& powgroup = GetParam(a);

                    if(powgroup.GetOpcode() == cPow
                    && powgroup.GetParam(1).IsImmed())
                    {
                        const CodeTree<Value_t>& exp_param = powgroup.GetParam(1);
                        Value_t exponent = exp_param.GetImmed();
                        if(fp_equal(exponent, Value_t(-1)))
                        {
                            CopyOnWrite();
                            div_params.push_back(GetParam(a).GetParam(0));
                            DelParam(a); // delete the pow group
                        }
                        else if(exponent < Value_t(0) && isInteger(exponent))
                        {
                            CodeTree<Value_t> edited_powgroup;
                            edited_powgroup.SetOpcode(cPow);
                            edited_powgroup.AddParam(powgroup.GetParam(0));
                            edited_powgroup.AddParam(CodeTreeImmed( -exponent ));
                            edited_powgroup.Rehash();
                            div_params.push_back(edited_powgroup);
                            CopyOnWrite();
                            DelParam(a); // delete the pow group
                        }
                    }
                    else if(powgroup.GetOpcode() == cLog2 && !found_log2.IsDefined())
                    {
                        found_log2 = powgroup.GetParam(0);
                        CopyOnWrite();
                        DelParam(a);
                    }
                    else if(powgroup.GetOpcode() == cLog2by && !found_log2by.IsDefined())
                    {
                        found_log2by = powgroup;
                        CopyOnWrite();
                        DelParam(a);
                    }
                }
                if(!div_params.empty())
                {
                    changed = true;

                    CodeTree<Value_t> divgroup;
                    divgroup.SetOpcode(cMul);
                    divgroup.SetParamsMove(div_params);
                    divgroup.Rehash(); // will reduce to div_params[0] if only one item
                    CodeTree<Value_t> mulgroup;
                    mulgroup.SetOpcode(cMul);
                    mulgroup.SetParamsMove(GetParams());
                    mulgroup.Rehash(); // will reduce to 1.0 if none remained in this cMul
                    if(mulgroup.IsImmed() && fp_equal(mulgroup.GetImmed(), Value_t(1)))
                    {
                        SetOpcode(cInv);
                        AddParamMove(divgroup);
                    }
                    /*else if(mulgroup.IsImmed() && fp_equal(mulgroup.GetImmed(), Value_t(-1)))
                    {
                        CodeTree<Value_t> invgroup;
                        invgroup.SetOpcode(cInv);
                        invgroup.AddParamMove(divgroup);
                        invgroup.Rehash();
                        SetOpcode(cNeg);
                        AddParamMove(invgroup);
                    }*/
                    else
                    {
                        if(mulgroup.GetDepth() >= divgroup.GetDepth())
                        {
                            SetOpcode(cDiv);
                            AddParamMove(mulgroup);
                            AddParamMove(divgroup);
                        }
                        else
                        {
                            SetOpcode(cRDiv);
                            AddParamMove(divgroup);
                            AddParamMove(mulgroup);
                        }
                    }
                }
                if(found_log2.IsDefined())
                {
                    CodeTree<Value_t> mulgroup;
                    mulgroup.SetOpcode(GetOpcode());
                    mulgroup.SetParamsMove(GetParams());
                    mulgroup.Rehash();
                    while(mulgroup.RecreateInversionsAndNegations(prefer_base2))
                        mulgroup.FixIncompleteHashes();
                    SetOpcode(cLog2by);
                    AddParamMove(found_log2);
                    AddParamMove(mulgroup);
                    changed = true;
                }
                if(found_log2by.IsDefined())
                {
                    CodeTree<Value_t> mulgroup;
                    mulgroup.SetOpcode(cMul);
                    mulgroup.AddParamMove(found_log2by.GetParam(1));
                    mulgroup.AddParamsMove(GetParams());
                    mulgroup.Rehash();
                    while(mulgroup.RecreateInversionsAndNegations(prefer_base2))
                        mulgroup.FixIncompleteHashes();
                    DelParams();
                    SetOpcode(cLog2by);
                    AddParamMove(found_log2by.GetParam(0));
                    AddParamMove(mulgroup);
                    changed = true;
                }
                break;
            }
            case cAdd:
            {
                std::vector<CodeTree<Value_t> > sub_params;

                for(size_t a = GetParamCount(); a-- > 0; )
                    if(GetParam(a).GetOpcode() == cMul)
                    {
                        bool is_signed = false; // if the mul group has a -1 constant...

                    Recheck_RefCount_Mul:;
                        CodeTree<Value_t>& mulgroup = GetParam(a);
                        bool needs_cow = GetRefCount() > 1;

                        for(size_t b=mulgroup.GetParamCount(); b-- > 0; )
                        {
                            if(mulgroup.GetParam(b).IsImmed())
                            {
                                Value_t factor = mulgroup.GetParam(b).GetImmed();
                                if(fp_equal(factor, Value_t(-1)))
                                {
                                    if(needs_cow) { CopyOnWrite(); goto Recheck_RefCount_Mul; }
                                    mulgroup.CopyOnWrite();
                                    mulgroup.DelParam(b);
                                    is_signed = !is_signed;
                                }
                                else if(fp_equal(factor, Value_t(-2)))
                                {
                                    if(needs_cow) { CopyOnWrite(); goto Recheck_RefCount_Mul; }
                                    mulgroup.CopyOnWrite();
                                    mulgroup.DelParam(b);
                                    mulgroup.AddParam( CodeTreeImmed( Value_t(2) ) );
                                    is_signed = !is_signed;
                                }
                            }
                        }
                        if(is_signed)
                        {
                            mulgroup.Rehash();
                            sub_params.push_back(mulgroup);
                            DelParam(a);
                        }
                    }
                    else if(GetParam(a).GetOpcode() == cDiv && !IsIntType<Value_t>::result)
                    {
                        bool is_signed = false;

                    Recheck_RefCount_Div:;
                        CodeTree<Value_t>& divgroup = GetParam(a);
                        bool needs_cow = GetRefCount() > 1;

                        if(divgroup.GetParam(0).IsImmed())
                        {
                            if(fp_equal(divgroup.GetParam(0).GetImmed(), Value_t(-1)))
                            {
                                if(needs_cow) { CopyOnWrite(); goto Recheck_RefCount_Div; }
                                divgroup.CopyOnWrite();
                                divgroup.DelParam(0);
                                divgroup.SetOpcode(cInv);
                                is_signed = !is_signed;
                            }
                        }
                        if(is_signed)
                        {
                            if(needs_cow) { CopyOnWrite(); goto Recheck_RefCount_Div; }
                            divgroup.Rehash();
                            sub_params.push_back(divgroup);
                            DelParam(a);
                        }
                    }
                    else if(GetParam(a).GetOpcode() == cRDiv && !IsIntType<Value_t>::result)
                    {
                        bool is_signed = false;

                    Recheck_RefCount_RDiv:;
                        CodeTree<Value_t>& divgroup = GetParam(a);
                        bool needs_cow = GetRefCount() > 1;

                        if(divgroup.GetParam(1).IsImmed())
                        {
                            if(fp_equal(divgroup.GetParam(1).GetImmed(), Value_t(-1)))
                            {
                                if(needs_cow) { CopyOnWrite(); goto Recheck_RefCount_RDiv; }
                                divgroup.CopyOnWrite();
                                divgroup.DelParam(1);
                                divgroup.SetOpcode(cInv);
                                is_signed = !is_signed;
                            }
                        }
                        if(is_signed)
                        {
                            if(needs_cow) { CopyOnWrite(); goto Recheck_RefCount_RDiv; }
                            divgroup.Rehash();
                            sub_params.push_back(divgroup);
                            DelParam(a);
                        }
                    }
                if(!sub_params.empty())
                {
                  #ifdef DEBUG_SUBSTITUTIONS
                    printf("Will make a Sub conversion in:\n"); fflush(stdout);
                    DumpTreeWithIndent(*this);
                  #endif
                    CodeTree<Value_t> subgroup;
                    subgroup.SetOpcode(cAdd);
                    subgroup.SetParamsMove(sub_params);
                    subgroup.Rehash(); // will reduce to sub_params[0] if only one item
                    CodeTree<Value_t> addgroup;
                    addgroup.SetOpcode(cAdd);
                    addgroup.SetParamsMove(GetParams());
                    addgroup.Rehash(); // will reduce to 0.0 if none remained in this cAdd
                    if(addgroup.IsImmed() && fp_equal(addgroup.GetImmed(), Value_t(0)))
                    {
                        SetOpcode(cNeg);
                        AddParamMove(subgroup);
                    }
                    else
                    {
                        if(addgroup.GetDepth() == 1)
                        {
                            /* 5 - (x+y+z) is best expressed as rsub(x+y+z, 5);
                             * this has lowest stack usage.
                             * This is identified by addgroup having just one member.
                             */
                            SetOpcode(cRSub);
                            AddParamMove(subgroup);
                            AddParamMove(addgroup);
                        }
                        else if(subgroup.GetOpcode() == cAdd)
                        {
                            /* a+b-(x+y+z) is expressed as a+b-x-y-z.
                             * Making a long chain of cSubs is okay, because the
                             * cost of cSub is the same as the cost of cAdd.
                             * Thus we get the lowest stack usage.
                             * This approach cannot be used for cDiv.
                             */
                            SetOpcode(cSub);
                            AddParamMove(addgroup);
                            AddParamMove(subgroup.GetParam(0));
                            for(size_t a=1; a<subgroup.GetParamCount(); ++a)
                            {
                                CodeTree<Value_t> innersub;
                                innersub.SetOpcode(cSub);
                                innersub.SetParamsMove(GetParams());
                                innersub.Rehash(false);
                                //DelParams();
                                AddParamMove(innersub);
                                AddParamMove(subgroup.GetParam(a));
                            }
                        }
                        else
                        {
                            SetOpcode(cSub);
                            AddParamMove(addgroup);
                            AddParamMove(subgroup);
                        }
                    }
                  #ifdef DEBUG_SUBSTITUTIONS
                    printf("After Sub conversion:\n"); fflush(stdout);
                    DumpTreeWithIndent(*this);
                  #endif
                }
                break;
            }
            case cPow:
            {
                const CodeTree<Value_t>& p0 = GetParam(0);
                const CodeTree<Value_t>& p1 = GetParam(1);
                if(p1.IsImmed())
                {
                    if(p1.GetImmed() != Value_t(0)
                    && !isInteger(p1.GetImmed()))
                    {
                        PowiResolver::PowiResult
                            r = PowiResolver::CreatePowiResult(fp_abs(p1.GetImmed()));

                        if(r.resulting_exponent != 0)
                        {
                            bool signed_chain = false;

                            if(p1.GetImmed() < Value_t(0)
                            && r.sep_list[0] == 0
                            && r.n_int_sqrt > 0)
                            {
                                // If one of the internal sqrts can be changed into rsqrt
                                signed_chain = true;
                            }

                        #ifdef DEBUG_POWI
                            printf("Will resolve powi %Lg as powi(chain(%d,%d),%ld)",
                                (long double) fp_abs(p1.GetImmed()),
                                r.n_int_sqrt,
                                r.n_int_cbrt,
                                r.resulting_exponent);
                            for(unsigned n=0; n<PowiResolver::MaxSep; ++n)
                            {
                                if(r.sep_list[n] == 0) break;
                                int n_sqrt = r.sep_list[n] % PowiResolver::MaxOp;
                                int n_cbrt = r.sep_list[n] / PowiResolver::MaxOp;
                                printf("*chain(%d,%d)", n_sqrt,n_cbrt);
                            }
                            printf("\n");
                        #endif

                            CodeTree<Value_t> source_tree = GetParam(0);

                            CodeTree<Value_t> pow_item = source_tree;
                            pow_item.CopyOnWrite();
                            ChangeIntoRootChain(pow_item,
                                signed_chain,
                                r.n_int_sqrt,
                                r.n_int_cbrt);
                            pow_item.Rehash();

                            CodeTree<Value_t> pow;
                            if(r.resulting_exponent != 1)
                            {
                                pow.SetOpcode(cPow);
                                pow.AddParamMove(pow_item);
                                pow.AddParam(CodeTreeImmed( Value_t(r.resulting_exponent) ));
                            }
                            else
                                pow.swap(pow_item);

                            CodeTree<Value_t> mul;
                            mul.SetOpcode(cMul);
                            mul.AddParamMove(pow);

                            for(unsigned n=0; n<PowiResolver::MaxSep; ++n)
                            {
                                if(r.sep_list[n] == 0) break;
                                int n_sqrt = r.sep_list[n] % PowiResolver::MaxOp;
                                int n_cbrt = r.sep_list[n] / PowiResolver::MaxOp;

                                CodeTree<Value_t> mul_item = source_tree;
                                mul_item.CopyOnWrite();
                                ChangeIntoRootChain(mul_item, false, n_sqrt, n_cbrt);
                                mul_item.Rehash();
                                mul.AddParamMove(mul_item);
                            }

                            if(p1.GetImmed() < Value_t(0) && !signed_chain)
                            {
                                mul.Rehash();
                                SetOpcode(cInv);
                                SetParamMove(0, mul);
                                DelParam(1);
                            }
                            else
                            {
                                SetOpcode(cMul);
                                SetParamsMove(mul.GetParams());
                            }
                        #ifdef DEBUG_POWI
                            DumpTreeWithIndent(*this);
                        #endif
                            changed = true;
                            break;
                        }
                    }
                }
                if(GetOpcode() == cPow
                && (!p1.IsImmed()
                 || !isLongInteger(p1.GetImmed())
                 || !IsOptimizableUsingPowi<Value_t>( makeLongInteger(p1.GetImmed()) )))
                {
                    if(p0.IsImmed() && p0.GetImmed() > Value_t(0.0))
                    {
                        // Convert into cExp or Exp2.
                        //    x^y = exp(log(x) * y) =
                        //    Can only be done when x is positive, though.
                        if(prefer_base2)
                        {
                            Value_t mulvalue = fp_log2( p0.GetImmed() );
                            if(fp_equal(mulvalue, Value_t(1)))
                            {
                                // exp2(1)^x becomes exp2(x)
                                DelParam(0);
                            }
                            else
                            {
                                // exp2(4)^x becomes exp2(4*x)
                                CodeTree<Value_t> exponent;
                                exponent.SetOpcode(cMul);
                                exponent.AddParam( CodeTreeImmed<Value_t>( mulvalue ) );
                                exponent.AddParam(p1);
                                exponent.Rehash();
                                SetParamMove(0, exponent);
                                DelParam(1);
                            }
                            SetOpcode(cExp2);
                            changed = true;
                        }
                        else
                        {
                            Value_t mulvalue = fp_log( p0.GetImmed() );
                            if(fp_equal(mulvalue, Value_t(1)))
                            {
                                // exp(1)^x becomes exp(x)
                                DelParam(0);
                            }
                            else
                            {
                                // exp(4)^x becomes exp(4*x)
                                CodeTree<Value_t> exponent;
                                exponent.SetOpcode(cMul);
                                exponent.AddParam( CodeTreeImmed<Value_t>( mulvalue ) );
                                exponent.AddParam(p1);
                                exponent.Rehash();
                                SetParamMove(0, exponent);
                                DelParam(1);
                            }
                            SetOpcode(cExp);
                            changed = true;
                        }
                    }
                    else if(GetPositivityInfo(p0) == IsAlways)
                    {
                        if(prefer_base2)
                        {
                            CodeTree<Value_t> log;
                            log.SetOpcode(cLog2);
                            log.AddParam(p0);
                            log.Rehash();
                            CodeTree<Value_t> exponent;
                            exponent.SetOpcode(cMul);
                            exponent.AddParam(p1);
                            exponent.AddParamMove(log);
                            exponent.Rehash();
                            SetOpcode(cExp2);
                            SetParamMove(0, exponent);
                            DelParam(1);
                            changed = true;
                        }
                        else
                        {
                            CodeTree<Value_t> log;
                            log.SetOpcode(cLog);
                            log.AddParam(p0);
                            log.Rehash();
                            CodeTree<Value_t> exponent;
                            exponent.SetOpcode(cMul);
                            exponent.AddParam(p1);
                            exponent.AddParamMove(log);
                            exponent.Rehash();
                            SetOpcode(cExp);
                            SetParamMove(0, exponent);
                            DelParam(1);
                            changed = true;
                        }
                    }
                }
                break;
            }
            case cDiv:
            {
                // Change 1/x into inv(x)
                // Needed in integer mode, no other code does it.
                if(GetParam(0).IsImmed()
                && fp_equal(GetParam(0).GetImmed(), Value_t(1)))
                {
                    SetOpcode(cInv);
                    DelParam(0);
                }
                break;
            }

            default: break;
        }

        if(changed)
            goto exit_changed;

        return changed;
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template \
    bool CodeTree<type>::RecreateInversionsAndNegations(bool prefer_base2);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/cse.cc"
// line removed for fpoptimizer.cc: #include "bytecodesynth.hh"
// line removed for fpoptimizer.cc: #include "codetree.hh"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

//#define DEBUG_SUBSTITUTIONS_CSE

namespace
{
    using namespace FPoptimizer_CodeTree;

    class TreeCountItem
    {
        size_t n_occurrences;
        size_t n_as_cos_param;
        size_t n_as_sin_param;
        size_t n_as_tan_param;
        size_t n_as_cosh_param;
        size_t n_as_sinh_param;
        size_t n_as_tanh_param;
    public:
        TreeCountItem() :
            n_occurrences(0),
            n_as_cos_param(0),
            n_as_sin_param(0),
            n_as_tan_param(0),
            n_as_cosh_param(0),
            n_as_sinh_param(0),
            n_as_tanh_param(0) { }

        void AddFrom(OPCODE op)
        {
            n_occurrences += 1;
            if(op == cCos) ++n_as_cos_param;
            if(op == cSin) ++n_as_sin_param;
            if(op == cSec) ++n_as_cos_param;
            if(op == cCsc) ++n_as_sin_param;
            if(op == cTan) ++n_as_tan_param;
            if(op == cCot) ++n_as_tan_param;
            if(op == cSinh) ++n_as_sinh_param;
            if(op == cCosh) ++n_as_cosh_param;
            if(op == cTanh) ++n_as_tanh_param;
        }

        size_t GetCSEscore() const
        {
            //size_t n_sincos = std::min(n_as_cos_param, n_as_sin_param);
            size_t result = n_occurrences;// - n_sincos;
            return result;
        }

        /* Calculate whether a sincos() would be useful.
         * Return values: 0 = not useful
         *                1,2 = yes
         * 1 = the tree is always a sin/cos parameter,
         *     so once a sincos() is synthesized, the
         *     tree itself does not need to be synthesized
         */
        int NeedsSinCos() const
        {
            bool always_sincostan =
                (n_occurrences == (n_as_cos_param + n_as_sin_param + n_as_tan_param));
            if((n_as_tan_param && (n_as_sin_param || n_as_cos_param))
            || (n_as_sin_param && n_as_cos_param))
            {
                if(always_sincostan)
                    return 1;
                return 2;
            }
            return 0;
        }

        /* Calculate whether a sinhcosh() would be useful.
         */
        int NeedsSinhCosh() const
        {
            bool always_sincostan =
                (n_occurrences == (n_as_cosh_param + n_as_sinh_param + n_as_tanh_param));
            if((n_as_tanh_param && (n_as_sinh_param || n_as_cosh_param))
            || (n_as_sinh_param && n_as_cosh_param))
            {
                if(always_sincostan)
                    return 1;
                return 2;
            }
            return 0;
        }

        size_t MinimumDepth() const
        {
            size_t n_sincos   = std::min(n_as_cos_param, n_as_sin_param);
            size_t n_sinhcosh = std::min(n_as_cosh_param, n_as_sinh_param);
            if(n_sincos == 0 && n_sinhcosh == 0)
                return 2;
            return 1;
        }
    };

    template<typename Value_t>
    class TreeCountType:
        public std::multimap<fphash_t,  std::pair<TreeCountItem, CodeTree<Value_t> > >
    {
    };

    template<typename Value_t>
    void FindTreeCounts(
        TreeCountType<Value_t>& TreeCounts,
        const CodeTree<Value_t>& tree,
        OPCODE parent_opcode,
        bool skip_root = false)
    {
        typename TreeCountType<Value_t>::iterator
            i = TreeCounts.lower_bound(tree.GetHash());
        if(!skip_root)
        {
            bool found = false;
            for(; i != TreeCounts.end() && i->first == tree.GetHash(); ++i)
            {
                if(tree.IsIdenticalTo( i->second.second ) )
                {
                    i->second.first.AddFrom(parent_opcode);
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                TreeCountItem count;
                count.AddFrom(parent_opcode);
                TreeCounts.insert(i, std::make_pair(tree.GetHash(),
                    std::make_pair(count, tree)));
            }
        }
        for(size_t a=0; a<tree.GetParamCount(); ++a)
            FindTreeCounts(TreeCounts, tree.GetParam(a),
                           tree.GetOpcode());
    }

    struct BalanceResultType
    {
        bool BalanceGood;
        bool FoundChild;
    };

    template<typename Value_t>
    BalanceResultType IfBalanceGood(const CodeTree<Value_t>& root,
                                    const CodeTree<Value_t>& child)
    {
        if(root.IsIdenticalTo(child))
        {
            BalanceResultType result = {true,true};
            return result;
        }

        BalanceResultType result = {true,false};

        if(root.GetOpcode() == cIf
        || root.GetOpcode() == cAbsIf)
        {
            BalanceResultType cond    = IfBalanceGood(root.GetParam(0), child);
            BalanceResultType branch1 = IfBalanceGood(root.GetParam(1), child);
            BalanceResultType branch2 = IfBalanceGood(root.GetParam(2), child);

            if(cond.FoundChild || branch1.FoundChild || branch2.FoundChild)
                { result.FoundChild = true; }

            // balance is good if:
            //      branch1.found = branch2.found OR (cond.found AND cond.goodbalance)
            // AND  cond.goodbalance OR (branch1.found AND branch2.found)
            // AND  branch1.goodbalance OR (cond.found AND cond.goodbalance)
            // AND  branch2.goodbalance OR (cond.found AND cond.goodbalance)

            result.BalanceGood =
                (   (branch1.FoundChild == branch2.FoundChild)
                 || (cond.FoundChild && cond.BalanceGood) )
             && (cond.BalanceGood || (branch1.FoundChild && branch2.FoundChild))
             && (branch1.BalanceGood || (cond.FoundChild && cond.BalanceGood))
             && (branch2.BalanceGood || (cond.FoundChild && cond.BalanceGood));
        }
        else
        {
            bool has_bad_balance        = false;
            bool has_good_balance_found = false;

            // Balance is bad if one of the children has bad balance
            // Unless one of the children has good balance & found

            for(size_t b=root.GetParamCount(), a=0; a<b; ++a)
            {
                BalanceResultType tmp = IfBalanceGood(root.GetParam(a), child);
                if(tmp.FoundChild)
                    result.FoundChild = true;

                if(tmp.BalanceGood == false)
                    has_bad_balance = true;
                else if(tmp.FoundChild)
                    has_good_balance_found = true;

                // if the expression is
                //   if(x, sin(x), 0) + sin(x)
                // then sin(x) is a good subexpression
                // even though it occurs in unbalance.
            }
            if(has_bad_balance && !has_good_balance_found)
                result.BalanceGood = false;
        }
        return result;
    }

    template<typename Value_t>
    bool ContainsOtherCandidates(
        const CodeTree<Value_t>& within,
        const CodeTree<Value_t>& tree,
        const FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth,
        const TreeCountType<Value_t>& TreeCounts)
    {
        for(size_t b=tree.GetParamCount(), a=0; a<b; ++a)
        {
            const CodeTree<Value_t>& leaf = tree.GetParam(a);

            typename TreeCountType<Value_t>::iterator synth_it;
            for(typename TreeCountType<Value_t>::const_iterator
                i = TreeCounts.begin();
                i != TreeCounts.end();
                ++i)
            {
                if(i->first != leaf.GetHash())
                    continue;

                const TreeCountItem& occ  = i->second.first;
                size_t          score     = occ.GetCSEscore();
                const CodeTree<Value_t>& candidate = i->second.second;

                // It must not yet have been synthesized
                if(synth.Find(candidate))
                    continue;

                // And it must not be a simple expression
                // Because cImmed, VarBegin are faster than cFetch
                if(leaf.GetDepth() < occ.MinimumDepth())
                    continue;

                // It must always occur at least twice
                if(score < 2)
                    continue;

                // And it must either appear on both sides
                // of a cIf, or neither
                if(IfBalanceGood(within, leaf).BalanceGood == false)
                    continue;

                return true;
            }
            if(ContainsOtherCandidates(within, leaf, synth, TreeCounts))
                return true;
        }
        return false;
    }

    template<typename Value_t>
    bool IsDescendantOf(const CodeTree<Value_t>& parent, const CodeTree<Value_t>& expr)
    {
        for(size_t a=0; a<parent.GetParamCount(); ++a)
            if(parent.GetParam(a).IsIdenticalTo(expr))
                return true;

        for(size_t a=0; a<parent.GetParamCount(); ++a)
            if(IsDescendantOf(parent.GetParam(a), expr))
                return true;

        return false;
    }

    template<typename Value_t>
    bool GoodMomentForCSE(const CodeTree<Value_t>& parent, const CodeTree<Value_t>& expr)
    {
        if(parent.GetOpcode() == cIf)
            return true;

        // Good if it's one of our direct children
        // Bad if it is a descendant of only one of our children

        for(size_t a=0; a<parent.GetParamCount(); ++a)
            if(parent.GetParam(a).IsIdenticalTo(expr))
                return true;

        size_t leaf_count = 0;
        for(size_t a=0; a<parent.GetParamCount(); ++a)
            if(IsDescendantOf(parent.GetParam(a), expr))
                ++leaf_count;

        return leaf_count != 1;
    }

}

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    size_t CodeTree<Value_t>::SynthCommonSubExpressions(
        FPoptimizer_ByteCode::ByteCodeSynth<Value_t>& synth) const
    {
        if(GetParamCount() == 0) return 0; // No subexpressions to synthesize.

        size_t stacktop_before = synth.GetStackTop();

        /* Find common subtrees */
        TreeCountType<Value_t> TreeCounts;
        FindTreeCounts(TreeCounts, *this, GetOpcode(), true);

        /* Synthesize some of the most common ones */
        for(;;)
        {
            size_t best_score = 0;
    #ifdef DEBUG_SUBSTITUTIONS_CSE
            std::cout << "Finding a CSE candidate, root is:" << std::endl;
            DumpHashes(*this);
    #endif
            typename TreeCountType<Value_t>::iterator cs_it ( TreeCounts.end() );
            for(typename TreeCountType<Value_t>::iterator
                j = TreeCounts.begin();
                j != TreeCounts.end(); )
            {
                typename TreeCountType<Value_t>::iterator i( j++ );

                const TreeCountItem& occ  = i->second.first;
                size_t          score     = occ.GetCSEscore();
                const CodeTree<Value_t>& tree = i->second.second;

    #ifdef DEBUG_SUBSTITUTIONS_CSE
                std::cout << "Score " << score << ":\n" << std::flush;
                DumpTreeWithIndent(tree);
    #endif

                // It must not yet have been synthesized
                if(synth.Find(tree))
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // And it must not be a simple expression
                // Because cImmed, VarBegin are faster than cFetch
                if(tree.GetDepth() < occ.MinimumDepth())
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // It must always occur at least twice
                if(score < 2)
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // And it must either appear on both sides
                // of a cIf, or neither
                if(IfBalanceGood(*this, tree).BalanceGood == false)
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // It must not contain other candidates
                if(ContainsOtherCandidates(*this, tree, synth, TreeCounts))
                {
                    // Don't erase it; it may be a proper candidate later
                    continue;
                }

                if(!GoodMomentForCSE(*this, tree))
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // Is a candidate.
                score *= tree.GetDepth();
                if(score > best_score)
                    { best_score = score; cs_it = i; }
            }

            if(best_score <= 0)
            {
    #ifdef DEBUG_SUBSTITUTIONS_CSE
                std::cout << "No more CSE candidates.\n" << std::flush;
    #endif
                break; // Didn't find anything.
            }

            //const TreeCountItem& occ    = cs_it->second.first;
            const CodeTree<Value_t>& tree = cs_it->second.second;
    #ifdef DEBUG_SUBSTITUTIONS_CSE
            std::cout << "Found Common Subexpression:"; DumpTree<Value_t>(tree); std::cout << std::endl;
    #endif
          #if 0
            int needs_sincos   = occ.NeedsSinCos();
            int needs_sinhcosh = occ.NeedsSinhCosh();
            CodeTree<Value_t> sintree, costree, sinhtree, coshtree;
            if(needs_sincos)
            {
                sintree.AddParam(tree);
                sintree.SetOpcode(cSin);
                sintree.Rehash();
                costree.AddParam(tree);
                costree.SetOpcode(cCos);
                costree.Rehash();
                if(synth.Find(sintree) || synth.Find(costree))
                {
                    if(needs_sincos == 2)
                    {
                        // sin, cos already found, and we don't
                        // actually need _this_ tree by itself
                        TreeCounts.erase(cs_it);
                        continue;
                    }
                    needs_sincos = 0;
                }
            }
            if(needs_sinhcosh)
            {
                sinhtree.AddParam(tree);
                sinhtree.SetOpcode(cSinh);
                sinhtree.Rehash();
                coshtree.AddParam(tree);
                coshtree.SetOpcode(cCosh);
                coshtree.Rehash();
                if(synth.Find(sinhtree) || synth.Find(coshtree))
                {
                    if(needs_sinhcosh == 2)
                    {
                        // sinh, cosh already found, and we don't
                        // actually need _this_ tree by itself
                        TreeCounts.erase(cs_it);
                        continue;
                    }
                    needs_sinhcosh = 0;
                }
            }
          #endif

            /* Synthesize the selected tree */
            tree.SynthesizeByteCode(synth, false);
            TreeCounts.erase(cs_it);
    #ifdef DEBUG_SUBSTITUTIONS_CSE
            synth.template Dump<0> ();
            std::cout << "Done with Common Subexpression:"; DumpTree<Value_t>(tree); std::cout << std::endl;
    #endif
          #if 0
            if(needs_sincos)
            {
                if(needs_sincos == 2 || needs_sinhcosh)
                {
                    // make a duplicate of the value, since it
                    // is also needed in addition to the sin/cos.
                    synth.FindAndDup(tree);
                }
                synth.AddOperation(cSinCos, 1, 2);

                synth.StackTopIs(sintree, 1);
                synth.StackTopIs(costree, 0);
            }
            if(needs_sinhcosh)
            {
                if(needs_sincos) synth.FindAndDup(tree);
                if(needs_sinhcosh == 2)
                {
                    // make a duplicate of the value, since it
                    // is also needed in addition to the sin/cos.
                    synth.FindAndDup(tree);
                }
                synth.AddOperation(cSinhCosh, 1, 2);

                synth.StackTopIs(sinhtree, 1);
                synth.StackTopIs(coshtree, 0);
            }
          #endif
        }

        return synth.GetStackTop() - stacktop_before;
    }
}

/* 
// line removed for fpoptimizer.cc: #include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template \
    size_t CodeTree<type>::SynthCommonSubExpressions( \
        FPoptimizer_ByteCode::ByteCodeSynth<type>& synth) const;
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
 */

#endif

#line 1 "fpoptimizer/optimize_main.cc"
// line removed for fpoptimizer.cc: #include "fpconfig.hh"
// line removed for fpoptimizer.cc: #include "fparser.hh"
// line removed for fpoptimizer.cc: #include "extrasrc/fptypes.hh"

// line removed for fpoptimizer.cc: #include "codetree.hh"
// line removed for fpoptimizer.cc: #include "optimize.hh"

#ifndef FP_DUMMY_OPTIMIZER

template<typename Value_t>
void FunctionParserBase<Value_t>::Optimize()
{
    using namespace FPoptimizer_CodeTree;

    CopyOnWrite();

    //PrintByteCode(std::cout);
    /*std::fprintf(stderr,
        "O:refCount:%u mVarCount:%u mfuncPtrs:%u mFuncParsers:%u mByteCode:%u mImmed:%u\n",
        mData->mReferenceCounter,
        mData->mVariablesAmount,
        (unsigned)mData->mFuncPtrs.size(),
        (unsigned)mData->mFuncParsers.size(),
        (unsigned)mData->mByteCode.size(),
        (unsigned)mData->mImmed.size()
    );*/

    CodeTree<Value_t> tree;
    tree.GenerateFrom(*mData);

    FPoptimizer_Optimize::ApplyGrammars(tree);

    std::vector<unsigned> byteCode;
    std::vector<Value_t> immed;
    size_t stacktop_max = 0;
    tree.SynthesizeByteCode(byteCode, immed, stacktop_max);

    /*std::cout << std::flush;
    std::cerr << std::flush;
    fprintf(stderr, "Estimated stacktop %u\n", (unsigned)stacktop_max);
    fflush(stderr);*/

    if(mData->mStackSize != stacktop_max)
    {
        mData->mStackSize = unsigned(stacktop_max); // Note: Ignoring GCC warning here.
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
        mData->mStack.resize(stacktop_max);
#endif
    }

    mData->mByteCode.swap(byteCode);
    mData->mImmed.swap(immed);

    //PrintByteCode(std::cout);
}

#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) \
    template<> void FunctionParserBase< type >::Optimize() {}

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(MpfrFloat)
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(GmpInt)
#endif

#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(std::complex<double>)
#endif

#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(std::complex<float>)
#endif

#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(std::complex<long double>)
#endif

#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) \
    template void FunctionParserBase<type>::Optimize();

#ifndef FP_DISABLE_DOUBLE_TYPE
FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(double)
#endif

#ifdef FP_SUPPORT_FLOAT_TYPE
FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(float)
#endif

#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(long double)
#endif

#ifdef FP_SUPPORT_LONG_INT_TYPE
FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(long)
#endif

#endif //FP_DUMMY_OPTIMIZER


#endif

#include "instantiate_for_ad.hh"
