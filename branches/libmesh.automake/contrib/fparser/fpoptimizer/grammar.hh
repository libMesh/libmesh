#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH

#include <iostream>

#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

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
        /* BEGIN_EXPLICIT_INSTANTATIONS */
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
        /* END_EXPLICIT_INSTANTATIONS */
    }

    template<typename Value_t>
    void DumpParam(const ParamSpec& p, std::ostream& o = std::cout);

    template<typename Value_t>
    void DumpParams(unsigned paramlist, unsigned count, std::ostream& o = std::cout);
}

#endif
