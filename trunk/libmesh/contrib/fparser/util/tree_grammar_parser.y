%{
#define YYDEBUG 1
#define YYERROR_VERBOSE 1
#include <string.h> // for error reporting

#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

#include "../fpoptimizer/grammar.hh"
#include "../fpoptimizer/consts.hh"

#include "../fpoptimizer/grammar.cc"
/* ^Note: including .cc file here in order to be able
 *  to instantiate DumpParam and DumpParams for complex types.
 */

#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <complex>
#include <map>
#include <set>
#include <algorithm>
#include <assert.h>

#include "../lib/crc32.hh"

#ifdef __GNUC__
# define likely(x)       __builtin_expect(!!(x), 1)
# define unlikely(x)     __builtin_expect(!!(x), 0)
#else
# define likely(x)   (x)
# define unlikely(x) (x)
#endif

static const unsigned PARAM_INDEX_BITS = 10;

/*********/
using namespace FPoptimizer_Grammar;

class GrammarDumper;

static void yyerror(const char* msg);
static int yylex(union YYSTYPE* lval);

namespace
{
    /* This function generated with make_identifier_parser.cc */
    unsigned readOpcode(const char* input)
    {
        using namespace FUNCTIONPARSERTYPES;
#include "extrasrc/fp_identifier_parser.inc"
        return 0;
    }
}

namespace
{
    struct mycomplex
    {
        double real, imag;
    };
    mycomplex operator -(const mycomplex& v)
        { mycomplex res = {-v.real, -v.imag }; return res; }

    typedef std::complex<double> stdcomplex;
}

namespace GrammarData
{
    class ParamSpec;

    class MatchedParams
    {
    public:
        ParamMatchingType Type;
        std::vector<ParamSpec*> Params;
        unsigned RestHolderIndex;

    public:
        MatchedParams()                    : Type(PositionalParams), Params(), RestHolderIndex(0) { }
        MatchedParams(ParamMatchingType t) : Type(t),                Params(), RestHolderIndex(0) { }
        MatchedParams(ParamSpec* p)        : Type(PositionalParams), Params(), RestHolderIndex(0) { Params.push_back(p); }

        MatchedParams* SetType(ParamMatchingType t) { Type=t; return this; }
        MatchedParams* AddParam(ParamSpec* p) { Params.push_back(p); return this; }

        const std::vector<ParamSpec*>& GetParams() const { return Params; }

        void RecursivelySetDefaultParamMatchingType();
        bool EnsureNoRepeatedNamedHolders(std::set<unsigned>& used) const;
        bool EnsureNoRepeatedNamedHolders() const;
        bool EnsureNoVariableCoverageParams_InPositionalParamLists();

        unsigned CalcRequiredParamsCount() const;

        unsigned BuildDepMask();
        void BuildFinalDepMask();
    };

    class FunctionType
    {
    public:
        FUNCTIONPARSERTYPES::OPCODE Opcode;
        MatchedParams Params;
    public:
        FunctionType(FUNCTIONPARSERTYPES::OPCODE o, const MatchedParams& p)
            : Opcode(o), Params(p) { }

        void RecursivelySetDefaultParamMatchingType()
        {
            using namespace FUNCTIONPARSERTYPES;
            Params.RecursivelySetDefaultParamMatchingType();
            if((Opcode == cAdd || Opcode == cMul
            || Opcode == cAnd || Opcode == cOr
            || Opcode == cAbsAnd || Opcode == cAbsOr)
            && Params.Type == PositionalParams)
                Params.Type = SelectedParams;
        }

        bool EnsureNoRepeatedNamedHolders() const
            { return Params.EnsureNoRepeatedNamedHolders(); }
    };

    class ParamSpec
    {
    public:
        unsigned DepMask;

        SpecialOpcode Opcode;      // specifies the type of the function
        union
        {
            mycomplex ConstantValue;// for NumConstant
            unsigned Index;                 // for ParamHolder
            FunctionType* Func;             // for SubFunction
        };
        unsigned ImmedConstraint;
        bool     IsConst;                   // when SubFunction

    public:
        struct ParamHolderTag{};

        ParamSpec(FunctionType* f)
            : DepMask(),
              Opcode(SubFunction),
              Func(f),
              ImmedConstraint(0),
              IsConst(false)
        {
        }

        ParamSpec(mycomplex d, unsigned constraints)
            : DepMask(),
              Opcode(NumConstant),
              ConstantValue(d),
              ImmedConstraint(constraints),
              IsConst(true)
        {
        }

        ParamSpec(FUNCTIONPARSERTYPES::OPCODE o, const std::vector<ParamSpec*>& p)
            : DepMask(),
              Opcode(SubFunction),
              Func(new FunctionType(o, MatchedParams(PositionalParams))),
              ImmedConstraint(0),
              IsConst(true)
        {
            if(o == FUNCTIONPARSERTYPES::cNeg && p[0]->Opcode == NumConstant)
            {
                delete Func;
                Opcode        = NumConstant;
                ConstantValue  = -p[0]->ConstantValue;
                ImmedConstraint = p[0]->ImmedConstraint;
            }
            else
            {
                Func->Params.Params = p;
                /*
                if(o == cAdd && p[1]->Opcode == SubFunction
                             && p[1]->Func->Opcode == cNeg
                             && p.size() == 2)
                {
                    Func->Opcode = cSub;
                    Func->Params.Params[1] = p[1]->Func->Params.Params[0];
                } -- not done because ConstantFolding() cannot handle cSub
                */
            }
        }

        ParamSpec(unsigned i, ParamHolderTag)
            : DepMask(),
              Opcode(ParamHolder), Index(i),
              ImmedConstraint(0),
              IsConst(true)
        {
        }

/*
        // Order:
        //  NumConstant { ConstantValue }
        //  ParamHolder { Index }
        //  SubFunction { Opcode, IsConst }
        bool operator< (const ParamSpec& b) const
        {
            if(Opcode == NumConstant)
                return (b.Opcode == NumConstant)
                        ? ConstantValue < b.ConstantValue
                        : true;
            if(Opcode == ParamHolder)
                return (b.Opcode == ParamHolder)
                        ? Index < b.Index
                        : (b.Opcode == SubFunction)
                            ? true
                            : false;
            if(Opcode == SubFunction)
                return (b.Opcode == SubFunction)
                    ? (Func->Opcode != b.Func->Opcode
                         ? Func->Opcode < b.Func->Opcode
                         : IsConst < b.IsConst
                      )
                    : false;
            return false;
        }
        bool operator!= (const ParamSpec& b) const { return !operator==(b); }
        bool operator== (const ParamSpec& b) const
        {
            switch(Opcode)
            {
                case NumConstant:
                    return b.Opcode == Opcode && fp_equal(ConstantValue, b.ConstantValue);
                case ParamHolder:
                    return b.Opcode == Opcode && ImmedConstraint == b.ImmedConstraint
                        && b.DepMask == DepMask && Index == b.Index;
                case SubFunction:
                    if(b.Opcode != SubFunction) return false;
                    if(Func->Opcode != b.Func->Opcode) return false;
                    if(ImmedConstraint != b.ImmedConstraint) return false;
                    if(DepMask != b.DepMask) return false;
                    if(IsConst != b.IsConst) return false;
                    if(Func->Params.Type != b.Func->Params.Type
                    || Func->Params.RestHolderIndex != b.Func->Params.RestHolderIndex
                    || Func->Params.Params.size() != b.Func->Params.Params.size())
                        return false;
                    for(size_t a=0; a<Func->Params.Params.size(); ++a)
                        if(*Func->Params.Params[a] != *b.Func->Params.Params[a])
                            return false;
            }
            return true;
        }
*/
        ParamSpec* SetConstraint(unsigned mask)
            { ImmedConstraint |= mask; return this; }

        unsigned BuildDepMask();

        void RecursivelySetDefaultParamMatchingType()
        {
            if(Opcode == SubFunction)
                Func->RecursivelySetDefaultParamMatchingType();
        }
        bool VerifyIsConstant()
        {
            switch(Opcode)
            {
                case NumConstant: return true;
                case ParamHolder: return
                    (ImmedConstraint & ConstnessMask) == Constness_Const;
                case SubFunction:
                    if(!IsConst) return false; // subfunctions are not constant
            }
            // For const-subfunctions, all params must be const.
            for(size_t a=0; a<Func->Params.Params.size(); ++a)
                if(!Func->Params.Params[a]->VerifyIsConstant()) return false;
            return true;
        }

        bool EnsureNoRepeatedNamedHolders() const
        {
            if(Opcode != SubFunction) return true;
            MatchedParams tmp;
            tmp.Params = Func->Params.Params;
            return tmp.EnsureNoRepeatedNamedHolders();
        }

    private:
        ParamSpec(const ParamSpec&);
        ParamSpec& operator= (const ParamSpec&);
    };

    class Rule
    {
    public:
        friend class GrammarDumper;
        RuleType Type;

        FunctionType  Input;
        MatchedParams Replacement; // length should be 1 if ProduceNewTree is used
        unsigned SituationFlags;
    public:
        Rule(RuleType t, const FunctionType& f, const MatchedParams& r)
            : Type(t), Input(f), Replacement(r), SituationFlags(0)
        { }

        Rule(RuleType t, const FunctionType& f, ParamSpec* p)
            : Type(t), Input(f), Replacement(), SituationFlags(0)
        { Replacement.AddParam(p); }

        void BuildFinalDepMask()
        {
            Input.Params.BuildFinalDepMask();
            //Replacement.BuildFinalDepMask(); -- not needed, though not wrong either.
        }
        void SetSituationFlags(unsigned flags)
        {
            SituationFlags = flags;
        }
    };

    class Grammar
    {
    public:
        std::vector<Rule> rules;
    public:
        Grammar(): rules() { }

        void AddRule(const Rule& r) { rules.push_back(r); }
        void BuildFinalDepMask()
        {
            for(size_t a=0; a<rules.size(); ++a)
                rules[a].BuildFinalDepMask();
        }
    };

    ////////////////////

    void MatchedParams::RecursivelySetDefaultParamMatchingType()
    {
        Type = PositionalParams;
        if(RestHolderIndex != 0)
            Type = AnyParams;

        for(size_t a=0; a<Params.size(); ++a)
            Params[a]->RecursivelySetDefaultParamMatchingType();
    }

    bool MatchedParams::EnsureNoRepeatedNamedHolders(std::set<unsigned>& used) const
    {
        for(size_t a=0; a<Params.size(); ++a)
        {
            if(Params[a]->Opcode == ParamHolder)
            {
                unsigned index = Params[a]->Index;
                std::set<unsigned>::iterator i = used.lower_bound(index);
                if(i != used.end() && *i == index)
                    return false;
                used.insert(i, index);
            }
            if(Params[a]->Opcode == SubFunction)
                if(!Params[a]->Func->Params.EnsureNoRepeatedNamedHolders(used))
                    return false;
        }
        return true;
    }

    bool MatchedParams::EnsureNoRepeatedNamedHolders() const
    {
        std::set<unsigned> used;
        return EnsureNoRepeatedNamedHolders(used);
    }

    bool MatchedParams::EnsureNoVariableCoverageParams_InPositionalParamLists()
    {
        if(Type != PositionalParams
        && Type != SelectedParams) return true;

        if(RestHolderIndex != 0) return false;

        for(size_t a=0; a<Params.size(); ++a)
        {
            if(Params[a]->Opcode == SubFunction)
                if(!Params[a]->Func->Params.EnsureNoVariableCoverageParams_InPositionalParamLists())
                    return false;
        }
        return true;
    }
    unsigned MatchedParams::CalcRequiredParamsCount() const
    {
        return (unsigned)Params.size();
    }

    unsigned MatchedParams::BuildDepMask()
    {
        unsigned result = 0;
        for(size_t a=0; a<Params.size(); ++a)
            result |= Params[a]->BuildDepMask();
        return result;
    }

    void MatchedParams::BuildFinalDepMask()
    {
        unsigned all_bits = BuildDepMask();

        // For each bit that is set in all_bits, unset
        // all of them that are only set in one of the parameters.
        for(unsigned bit=1; all_bits >= bit; bit <<= 1)
            if(all_bits & bit)
            {
                unsigned count_found = 0;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    unsigned param_bitmask = Params[a]->DepMask;
                    if(param_bitmask & bit) ++count_found;
                }
                if(count_found <= 1)
                {
                    for(size_t a=0; a<Params.size(); ++a)
                        Params[a]->DepMask &= ~bit;
                }
            }

        // Recurse
        for(size_t a=0; a<Params.size(); ++a)
            if(Params[a]->Opcode == SubFunction)
                Params[a]->Func->Params.BuildFinalDepMask();
    }
}

namespace FPoptimizer_Grammar
{
    template<typename Value_t> // Used only by tree_grammar_parser.y
    bool ParamSpec_Compare(const void* aa, const void* bb, SpecialOpcode type)
    {
        switch(type)
        {
            case ParamHolder:
            {
                ParamSpec_ParamHolder& a = *(ParamSpec_ParamHolder*) aa;
                ParamSpec_ParamHolder& b = *(ParamSpec_ParamHolder*) bb;
                return a.constraints == b.constraints
                    && a.index       == b.index
                    && a.depcode     == b.depcode;
            }
            case NumConstant:
            {
                ParamSpec_NumConstant<Value_t>& a = *(ParamSpec_NumConstant<Value_t>*) aa;
                ParamSpec_NumConstant<Value_t>& b = *(ParamSpec_NumConstant<Value_t>*) bb;
                return a.constvalue == b.constvalue
                    && a.modulo == b.modulo;
            }
            case SubFunction:
            {
                ParamSpec_SubFunction& a = *(ParamSpec_SubFunction*) aa;
                ParamSpec_SubFunction& b = *(ParamSpec_SubFunction*) bb;
                return a.constraints    == b.constraints
                    && a.data.subfunc_opcode   == b.data.subfunc_opcode
                    && a.data.match_type       == b.data.match_type
                    && a.data.param_count      == b.data.param_count
                    && a.data.param_list       == b.data.param_list
                    && a.data.restholder_index == b.data.restholder_index
                    && a.depcode               == b.depcode;
            }
        }
        return true;
    }
}

GrammarData::Grammar grammar;
std::vector<ParamSpec> plist;
std::vector<Rule>      rlist;

struct RuleComparer
{
    bool operator() (const Rule& a, const Rule& b) const
    {
        if(a.match_tree.subfunc_opcode != b.match_tree.subfunc_opcode)
            return a.match_tree.subfunc_opcode < b.match_tree.subfunc_opcode;

        // Other rules to break ties
        if(a.situation_flags != b.situation_flags)
            return a.situation_flags < b.situation_flags;

        if(a.ruletype != b.ruletype)
            return a.ruletype < b.ruletype;

        if(a.match_tree.match_type != b.match_tree.match_type)
            return a.match_tree.match_type < b.match_tree.match_type;

        if(a.match_tree.param_count != b.match_tree.param_count)
            return a.match_tree.param_count < b.match_tree.param_count;

        if(a.repl_param_count != b.repl_param_count)
            return a.repl_param_count < b.repl_param_count;

        if(a.match_tree.param_list != b.match_tree.param_list)
            return a.match_tree.param_list < b.match_tree.param_list;

        if(a.repl_param_list != b.repl_param_list)
            return a.repl_param_list < b.repl_param_list;

        return false;
    }

    bool operator() (unsigned a, unsigned b) const
    {
        return this->operator() ( rlist[a], rlist[b] );
    }
};

class GrammarDumper
{
private:
    std::string GenName(const char* prefix)
    {
        static unsigned counter = 0;
        std::ostringstream tmp;
        tmp << prefix << ++counter;
        return tmp.str();
    }
private:
    std::map<std::string, size_t> n_index;

    std::vector<std::string>        nlist;
    std::map<std::string, Grammar>  glist;
public:
    GrammarDumper():
        n_index(),
        nlist(),glist()
    {
        plist.reserve(16384);
        nlist.reserve(16);
        rlist.reserve(16384);
    }

    unsigned ConvertNamedHolderNameIntoIndex(const std::string& n)
    {
        std::map<std::string, size_t>::const_iterator i = n_index.find(n);
        if(i != n_index.end()) return i->second;
        nlist.push_back(n);
        return n_index[n] = (unsigned)(nlist.size()-1);
    }
    size_t GetNumNamedHolderNames() const { return nlist.size(); }

    void DumpParamList(const std::vector<GrammarData::ParamSpec*>& Params,
                       unsigned&       param_count,
                       unsigned&       param_list)
    {
        param_count = (unsigned)Params.size();
        param_list  = 0;
        for(unsigned a=0; a<param_count; ++a)
        {
            ParamSpec p = CreateParam(*Params[a]);

            unsigned paramno = (unsigned)plist.size();

            for(size_t b = 0; b < plist.size(); ++b)
                if(plist[b].first == p.first
                && ParamSpec_Compare<stdcomplex>(plist[b].second, p.second, p.first))
                {
                    paramno = (unsigned)b;
                    break;
                }

            if(paramno == plist.size()) plist.push_back(p);

            param_list |= paramno << (a * PARAM_INDEX_BITS);
        }
    }

    ParamSpec CreateParam(const GrammarData::ParamSpec& p)
    {
        unsigned    pcount;
        unsigned    plist;
        switch(p.Opcode)
        {
            case SubFunction:
            {
                ParamSpec_SubFunction* result = new ParamSpec_SubFunction;
                result->constraints    = p.ImmedConstraint;
                result->data.subfunc_opcode = p.Func->Opcode;
                result->data.match_type     = p.Func->Params.Type;
                DumpParamList(p.Func->Params.Params, pcount, plist);
                result->data.param_count = pcount;
                result->data.param_list  = plist;
                result->depcode        = p.DepMask;
                result->data.restholder_index = p.Func->Params.RestHolderIndex;
                if(p.IsConst)
                {
                    result->data.match_type = GroupFunction;
                    result->constraints |= Constness_Const;
                }
                return std::make_pair(SubFunction, (void*)result);
            }
            case NumConstant:
            {
                typedef stdcomplex v;
                ParamSpec_NumConstant<v>* result = new ParamSpec_NumConstant<v>;
                result->constvalue     = v(p.ConstantValue.real, p.ConstantValue.imag);
                result->modulo         = p.ImmedConstraint;
                return std::make_pair(NumConstant, (void*)result);
            }
            case ParamHolder:
            {
                ParamSpec_ParamHolder* result = new ParamSpec_ParamHolder;
                result->constraints    = p.ImmedConstraint;
                result->index          = p.Index;
                result->depcode        = p.DepMask;
                return std::make_pair(ParamHolder, (void*)result);
            }
        }
        std::cout << "???\n";
        return std::make_pair(SubFunction, (void*) 0);
    }

    Rule CreateRule(const GrammarData::Rule& r)
    {
        unsigned min_params = r.Input.Params.CalcRequiredParamsCount();

        Rule ritem;
        memset(&ritem, 0, sizeof(ritem));
        //ritem.n_minimum_params          = min_params;
        ritem.ruletype                  = r.Type;
        ritem.situation_flags           = r.SituationFlags;
        ritem.match_tree.subfunc_opcode = r.Input.Opcode;
        ritem.match_tree.match_type     = r.Input.Params.Type;
        ritem.match_tree.restholder_index = r.Input.Params.RestHolderIndex;
        unsigned         pcount;
        unsigned         plist;
        DumpParamList(r.Input.Params.Params, pcount, plist);
        ritem.match_tree.param_count = pcount;
        ritem.match_tree.param_list  = plist;

        DumpParamList(r.Replacement.Params,  pcount, plist);
        ritem.repl_param_count = pcount;
        ritem.repl_param_list  = plist;
        return ritem;
    }

    void RegisterGrammar(const std::vector<GrammarData::Grammar>& gset)
    {
        using namespace FUNCTIONPARSERTYPES;
        std::vector<Rule> this_rules;

        for(size_t a=0; a<gset.size(); ++a)
        {
            const GrammarData::Grammar& g = gset[a];

            for(size_t a=0; a<g.rules.size(); ++a)
            {
                if(g.rules[a].Input.Opcode == cNop) continue;
                this_rules.push_back( CreateRule(g.rules[a]) );
            }
        }

        std::sort(this_rules.begin(), this_rules.end(),
                  RuleComparer());

        for(size_t a=0; a<this_rules.size(); ++a)
        {
            const Rule& r = this_rules[a];

            // Add to global rule list, unless it's already there
            bool dup=false;
            for(size_t c=0; c<rlist.size(); ++c)
                if(memcmp(&r, &rlist[c], sizeof(r)) == 0)
                {
                    // Already in global rule list...
                    dup = true;
                    break;
                }
            if(!dup)
                rlist.push_back(r);
        }
    }

    void DumpGrammar(const std::string& grammarname,
                     const std::vector<GrammarData::Grammar>& gset)
    {
        using namespace FUNCTIONPARSERTYPES;
        std::vector<unsigned> rule_list;

        std::vector<Rule> this_rules;

        for(size_t a=0; a<gset.size(); ++a)
        {
            const GrammarData::Grammar& g = gset[a];

            for(size_t a=0; a<g.rules.size(); ++a)
            {
                if(g.rules[a].Input.Opcode == cNop) continue;
                this_rules.push_back( CreateRule(g.rules[a]) );
            }
        }

        std::sort(this_rules.begin(), this_rules.end(),
                  RuleComparer());

        for(size_t a=0; a<this_rules.size(); ++a)
        {
            const Rule& r = this_rules[a];

            // Add to global rule list, unless it's already there
            bool dup=false;
            for(size_t c=0; c<rlist.size(); ++c)
                if(memcmp(&r, &rlist[c], sizeof(r)) == 0)
                {
                    // Already in global rule list...
                    // Add to grammar's rule list unless it's already there
                    dup = false;
                    for(size_t b=0; b<rule_list.size(); ++b)
                        if(c == rule_list[b])
                        {
                            dup = true;
                            break;
                        }
                    if(!dup)
                    {
                        // Global duplicate, but not yet in grammar.
                        rule_list.push_back(c);
                    }
                    dup = true;
                    break;
                }
            if(!dup)
            {
                // Not in global rule list. Add there and in grammar.
                rule_list.push_back( (unsigned) rlist.size() );
                rlist.push_back(r);
            }
        }

        Grammar& gitem = glist[grammarname];

        gitem.rule_count = (unsigned) rule_list.size();

        std::sort(rule_list.begin(), rule_list.end(),
                  RuleComparer());

        for(size_t a=0; a<rule_list.size(); ++a)
            gitem.rule_list[a] = rule_list[a];
    }

    static std::string ConstraintsToString(unsigned constraints)
    {
        std::ostringstream result;
        const char* sep = "";
        static const char s[] = " | ";
        switch( ImmedConstraint_Value( constraints & ValueMask ) )
        {
            case ValueMask: case Value_AnyNum: break;
            case Value_EvenInt: result << sep << "Value_EvenInt"; sep=s; break;
            case Value_OddInt: result << sep << "Value_OddInt"; sep=s; break;
            case Value_IsInteger: result << sep << "Value_IsInteger"; sep=s; break;
            case Value_NonInteger: result << sep << "Value_NonInteger"; sep=s; break;
            case Value_Logical: result << sep << "Value_Logical"; sep=s; break;
        }
        switch( ImmedConstraint_Sign( constraints & SignMask ) )
        {
            /*case SignMask:*/ case Sign_AnySign: break;
            case Sign_Positive: result << sep << "Sign_Positive"; sep=s; break;
            case Sign_Negative: result << sep << "Sign_Negative"; sep=s; break;
            case Sign_NoIdea:   result << sep << "Sign_NoIdea"; sep=s; break;
        }
        switch( ImmedConstraint_Oneness( constraints & OnenessMask ) )
        {
            case OnenessMask: case Oneness_Any: break;
            case Oneness_One: result << sep << "Oneness_One"; sep=s; break;
            case Oneness_NotOne: result << sep << "Oneness_NotOne"; sep=s; break;
        }
        switch( ImmedConstraint_Constness( constraints & ConstnessMask ) )
        {
            case ConstnessMask: case Oneness_Any: break;
            case Constness_Const: result << sep << "Constness_Const"; sep=s; break;
            case Constness_NotConst: result << sep << "Constness_NotConst"; sep=s; break;
        }
        if(!*sep) result << "0";
        return result.str();
    }
    static std::string ModuloToString(unsigned constraints)
    {
        std::ostringstream result;
        const char* sep = "";
        static const char s[] = " | ";
        switch( Modulo_Mode(constraints) )
        {
            case Modulo_None: break;
            case Modulo_Radians: result << sep << "Modulo_Radians"; sep=s; break;
        }
        if(!*sep) result << "0";
        return result.str();
    }

    static std::string ConstValueToString(const stdcomplex& value)
    {
        using namespace FUNCTIONPARSERTYPES;
        std::ostringstream result;
        result.precision(50);
        double dvalue = value.real();
        if(value.imag() != 0.0) goto NotAnyKnownConstant;
        #define Value_t double
        #define if_const(n) \
            if(fp_equal(dvalue, n)) result << #n; \
            else if(fp_equal(dvalue, -n)) result << "-" #n;
        if_const(fp_const_e<Value_t>())
        else if_const(fp_const_einv<Value_t>())
        else if_const(fp_const_twoe<Value_t>())
        else if_const(fp_const_twoeinv<Value_t>())
        else if_const(fp_const_pi<Value_t>())
        else if_const(fp_const_pihalf<Value_t>())
        else if_const(fp_const_twopi<Value_t>())
        else if_const(fp_const_log2<Value_t>())
        else if_const(fp_const_log2inv<Value_t>())
        else if_const(fp_const_log10<Value_t>())
        else if_const(fp_const_log10inv<Value_t>())
        else if_const(fp_const_rad_to_deg<Value_t>())
        else if_const(fp_const_deg_to_rad<Value_t>())
        #undef if_const
        #undef Value_t
        else
        {
        NotAnyKnownConstant:
            result << "Value_t(" << value.real() << ")";
            if(value.imag() != 0.0)
                result << " + fp_make_imag(Value_t(" << value.imag() << "))";
        }
        return result.str();
    }

    struct ParamCollection
    {
        std::vector<ParamSpec_ParamHolder>           plist_p;
        std::vector<ParamSpec_NumConstant<stdcomplex> >  plist_n;
        std::vector<ParamSpec_SubFunction>           plist_s;

        void Populate(const ParamSpec& param)
        {
            #define set(when, type, list, code) \
                case when: \
                  { for(size_t a=0; a<list.size(); ++a) \
                        if(ParamSpec_Compare<stdcomplex>(param.second, (const void*) &list[a], when)) \
                            return; \
                    list.push_back( *(type*) param.second ); \
                    code; \
                    break; }
            switch(param.first)
            {
                set(ParamHolder, ParamSpec_ParamHolder,         plist_p, {} );
                set(NumConstant, ParamSpec_NumConstant<stdcomplex>, plist_n, {} );
                set(SubFunction, ParamSpec_SubFunction,         plist_s,
                     ParamSpec_SubFunction* p = (ParamSpec_SubFunction*)param.second;
                     for(size_t a=0; a<p->data.param_count; ++a)
                         Populate( ParamSpec_Extract<stdcomplex>( p->data.param_list, a) );
                    );
            }
            #undef set
        }

        struct p_compare { int kind(
            const ParamSpec_ParamHolder& a,
            const ParamSpec_ParamHolder& b) const
        {
            if((a.index^2) != (b.index^2)) return (a.index^2) < (b.index^2) ? -1 : 1;
            // xor-2 is here to tweak the sorting order such that
            // the most used parameters (x,y) are first in the list,
            // resulting in smaller numbers for the parameter indexes,
            // and thus a smaller source code size for grammar data.
            return 0;
        } };
        struct n_compare { int kind(
            const ParamSpec_NumConstant<stdcomplex>& a,
            const ParamSpec_NumConstant<stdcomplex>& b) const
        {
            if(a.modulo != b.modulo) return a.modulo < b.modulo ? -1 : 1;
            double av = std::norm(a.constvalue), bv = std::norm(b.constvalue);
            if(a.constvalue.real() < 0) av = -av;
            if(b.constvalue.real() < 0) bv = -bv;
            if(av != bv) return av < bv ? -1 : 1;
            return 0;
        } };
        struct s_compare { int kind(
            const ParamSpec_SubFunction& a,
            const ParamSpec_SubFunction& b) const
        {
            unsigned a_opcode = a.data.subfunc_opcode;
            unsigned b_opcode = b.data.subfunc_opcode;

            if(a_opcode == FUNCTIONPARSERTYPES::cAdd) a_opcode = 2;
            else if(a_opcode == FUNCTIONPARSERTYPES::cMul) a_opcode = 3;
            else if(a_opcode == FUNCTIONPARSERTYPES::cPow) a_opcode = 4;
            else if(a_opcode == FUNCTIONPARSERTYPES::cNeg) a_opcode = 0;
            else if(a_opcode == FUNCTIONPARSERTYPES::cInv) a_opcode = 1;
            else a_opcode += 5;
            if(b_opcode == FUNCTIONPARSERTYPES::cAdd) b_opcode = 2;
            else if(b_opcode == FUNCTIONPARSERTYPES::cMul) b_opcode = 3;
            else if(b_opcode == FUNCTIONPARSERTYPES::cPow) b_opcode = 4;
            else if(b_opcode == FUNCTIONPARSERTYPES::cNeg) b_opcode = 0;
            else if(b_opcode == FUNCTIONPARSERTYPES::cInv) b_opcode = 1;
            else b_opcode += 5;

            if(a_opcode != b_opcode)
                return a_opcode < b_opcode ? -1 : 1;
            if(a.constraints != b.constraints)
                return a.constraints < b.constraints ? -1 : 1;
            if(a.data.match_type != b.data.match_type)
                return a.data.match_type < b.data.match_type ? -1 : 1;

            size_t min_param_count = std::min(a.data.param_count, b.data.param_count);

            for(size_t c=0; c< min_param_count; ++c)
            {
                ParamSpec aa = ParamSpec_Extract<stdcomplex>(a.data.param_list, (unsigned)c);
                ParamSpec bb = ParamSpec_Extract<stdcomplex>(b.data.param_list, (unsigned)c);
                if(aa.first != bb.first)
                    return aa.first < bb.first;
                switch(aa.first)
                {
                    case ParamHolder: {
                        int k = p_compare().kind
                            (*(const ParamSpec_ParamHolder*)aa.second,
                             *(const ParamSpec_ParamHolder*)bb.second);
                        if(k) return k;
                        break;
                   }case NumConstant: {
                        int k = n_compare().kind
                            (*(const ParamSpec_NumConstant<stdcomplex>*)aa.second,
                             *(const ParamSpec_NumConstant<stdcomplex>*)bb.second);
                        if(k) return k;
                        break;
                   }case SubFunction:{
                        int k = s_compare().kind
                            (*(const ParamSpec_SubFunction*)aa.second,
                             *(const ParamSpec_SubFunction*)bb.second);
                        if(k) return k;
                        break;
                }  }
            }
            if(a.data.param_count != b.data.param_count)
                return a.data.param_count < b.data.param_count ? -1 : 1;
            return 0;
        } };
        template<typename T>
        struct kind_compare
        {
            template<typename K>
            bool operator() (const K& a, const K& b) const
            {
                return T().kind(a,b) < 0;
            }
        };

        void Sort()
        {
            std::stable_sort(plist_p.begin(), plist_p.end(), kind_compare<p_compare>());
            std::stable_sort(plist_n.begin(), plist_n.end(), kind_compare<n_compare>());
            std::stable_sort(plist_s.begin(), plist_s.end(), kind_compare<s_compare>());
        }

        unsigned ParamPtrToParamIndex(unsigned paramlist, unsigned index) const
        {
            const ParamSpec& p = ParamSpec_Extract<stdcomplex> (paramlist, index);
            if(p.second)
            {
                #define set(when, list, c) \
                    case when: \
                        for(size_t a=0; a<list.size(); ++a) \
                            if(ParamSpec_Compare<stdcomplex> (p.second, (const void*)&list[a], when)) \
                                return (a + c##offset); \
                        break;
                unsigned Poffset = 0;
                unsigned Noffset = plist_p.size();
                unsigned Soffset = plist_n.size() + Noffset;
                switch(p.first)
                {
                    set(ParamHolder, plist_p, P);
                    set(NumConstant, plist_n, N);
                    set(SubFunction, plist_s, S);
                }
                #undef set
            }
            return (1 << 10)-1;
        }

        std::string ParamListToString(unsigned paramlist, unsigned paramcount) const
        {
            std::ostringstream result, comment;
            unsigned value = 0;
            for(unsigned p=0; p<paramcount; ++p)
            {
                unsigned index = ParamPtrToParamIndex(paramlist, p);
                if(p) comment << ',';
                comment << index;
                value += index << (p*PARAM_INDEX_BITS);
            }
            std::string commentstr = comment.str();
            commentstr.resize(3*3+2, ' ');
            result << "/*" << commentstr << "*/" << value;

            std::string res = result.str();
            if(res.size() < 25) res.resize(25, ' ');
            /* 999*x+999*x+999 = 15 characters */
            /* (*999,999,999*)1048551399 = 25 characters */
            return res;
        }
        std::string ParamHolderToString(const ParamSpec_ParamHolder& i) const
        {
            std::ostringstream result;
            result << "{" << i.index
                   << ", " << ConstraintsToString(i.constraints)
                   << ", 0x" << i.depcode
                   << "}";
            return result.str();
        }

        std::string NumConstantToString(const ParamSpec_NumConstant<stdcomplex>& i) const
        {
            std::ostringstream result;
            result << "{" << ConstValueToString(i.constvalue)
                   << ", " << ModuloToString(i.modulo)
                   << "}";
            return result.str();
        }

        std::string SubFunctionDataToString(const ParamSpec_SubFunctionData& i) const
        {
            std::ostringstream result;
            result << "{"  << i.param_count
                   <<  "," << ParamListToString(i.param_list, i.param_count)
                   << ", " << FP_GetOpcodeName(i.subfunc_opcode, true)
                   << ","  << (i.match_type == PositionalParams ? "PositionalParams"
                            :  i.match_type == SelectedParams   ? "SelectedParams  "
                            :  i.match_type == AnyParams        ? "AnyParams       "
                            :/*i.match_type == GroupFunction  ?*/ "GroupFunction   "
                            )
                   << "," << i.restholder_index
                   << "}";
            return result.str();
        }

        std::string SubFunctionToString(const ParamSpec_SubFunction& i) const
        {
            std::ostringstream result;
            result << "{" << SubFunctionDataToString(i.data)
                   << ", " << ConstraintsToString(i.constraints)
                   << ", 0x" << i.depcode
                   << "}";
            return result.str();
        }
    };

    ParamCollection collection;

    void Flush()
    {
        for(size_t a=0; a<rlist.size(); ++a)
        {
            for(unsigned b=0; b < rlist[a].match_tree.param_count; ++b)
                collection.Populate( ParamSpec_Extract<stdcomplex>(rlist[a].match_tree.param_list, b) );
            for(unsigned b=0; b < rlist[a].repl_param_count; ++b)
                collection.Populate( ParamSpec_Extract<stdcomplex>(rlist[a].repl_param_list, b) );
        }
        collection.Sort();

        std::cout << "/* BEGIN_EXPLICIT_INSTANTATIONS */\n";
        for(std::map<std::string, Grammar>::const_iterator
             i = glist.begin(); i != glist.end(); ++i)
            std::cout << "#define grammar_" << i->first << " grammar_" << i->first << "_tweak\n";
        std::cout <<
            "#include \"../fpoptimizer/grammar.hh\"\n";
        for(std::map<std::string, Grammar>::const_iterator
             i = glist.begin(); i != glist.end(); ++i)
            std::cout << "#undef grammar_" << i->first << "\n";
        std::cout << "/* END_EXPLICIT_INSTANTATIONS */\n";

        std::cout <<
            "\n"
            "using namespace FPoptimizer_Grammar;\n"
            "using namespace FUNCTIONPARSERTYPES;\n"
            "\n"
            "namespace\n"
            "{\n";

        {

        #define set(type, listprefix, list, c) \
            std::cout << \
            "    const ParamSpec_" #type " " listprefix #list "[" << collection.list.size() << "] =\n" \
            "    {\n"; \
            for(size_t a=0; a<collection.list.size(); ++a) \
            { \
                std::cout << "    /* " << offset++ << "\t*/ " \
                          << collection.type##ToString(collection.list[a]) \
                          << ", /* "; \
                FPoptimizer_Grammar::DumpParam<stdcomplex>( ParamSpec(type, (const void*) &collection.list[a]), std::cout); \
                std::cout << " */\n"; \
            } \
            std::cout << \
            "    };\n" \
            "\n";

        unsigned offset = 0;
        set(ParamHolder, "", plist_p, P) // Must be first one
        std::cout <<
            "    template<typename Value_t>\n"
            "    struct plist_n_container\n"
            "    {\n"
            "        static const ParamSpec_NumConstant<Value_t> plist_n[" << collection.plist_n.size() << "];\n"
            "    };\n"
            "    template<typename Value_t>\n";
        set(NumConstant, "<Value_t> plist_n_container<Value_t>::", plist_n, N)
        set(SubFunction, "", plist_s, S)

        std::cout <<
            "}\n";
        }

        #undef set

        std::cout <<
            "namespace FPoptimizer_Grammar\n"
            "{\n";
        std::cout <<
            "    const Rule grammar_rules[" << rlist.size() << "] =\n"
            "    {\n";
        for(size_t a=0; a<rlist.size(); ++a)
        {
            std::cout <<
            "        /* " << a << ":\t";
            ParamSpec_SubFunction tmp = {rlist[a].match_tree,0,0};
            if(rlist[a].situation_flags & OnlyForComplex)
                std::cout << "@C ";
            if(rlist[a].situation_flags & NotForComplex)
                std::cout << "@R ";
            if(rlist[a].situation_flags & LogicalContextOnly)
                std::cout << "@L ";
            if(rlist[a].situation_flags & NotForIntegers)
                std::cout << "@F ";
            if(rlist[a].situation_flags & OnlyForIntegers)
                std::cout << "@I ";
            FPoptimizer_Grammar::DumpParam<stdcomplex>
                ( ParamSpec(SubFunction, (const void*) &tmp) );
            switch(rlist[a].ruletype)
            {
                case ProduceNewTree:
                    std::cout <<
                    "\n"
                    "         *\t->\t";
                    FPoptimizer_Grammar::DumpParam<stdcomplex>(
                        ParamSpec_Extract<stdcomplex>(rlist[a].repl_param_list, 0) );
                    break;
                case ReplaceParams: default:
                    std::cout <<
                    "\n"
                    "         *\t:\t";
                    FPoptimizer_Grammar::DumpParams<stdcomplex>
                        ( rlist[a].repl_param_list, rlist[a].repl_param_count);
                    break;
            }
            std::cout <<
            "\n"
            "         */\t\t "
                        "{"
                        << (rlist[a].ruletype == ProduceNewTree  ? "ProduceNewTree"
                         :/*rlist[a].ruletype == ReplaceParams ?*/ "ReplaceParams "
                           )
                        << ", " << rlist[a].situation_flags
                        << ", " << rlist[a].repl_param_count
                        <<  "," << collection.ParamListToString(rlist[a].repl_param_list, rlist[a].repl_param_count)
                        << ", " << collection.SubFunctionDataToString(rlist[a].match_tree)
                        << "},\n";
        }
        std::cout <<
            "    };\n"
            <<
            "\n";
        for(std::map<std::string, Grammar>::const_iterator
             i = glist.begin(); i != glist.end(); ++i)
        {
            std::cout << "    struct grammar_" << i->first << "_type\n"
                         "    {\n"
                         "        unsigned c;\n"
                         "        unsigned short l[" << i->second.rule_count << "];\n"
                         "    };\n"
                         "    extern \"C\"\n"
                         "    {\n"
                         "        grammar_" << i->first << "_type grammar_" << i->first << " =\n"
                         "        {\n"
                         "            " << i->second.rule_count << ",\n"
                         "            { ";
            for(size_t p=0; p<i->second.rule_count; ++p)
            {
                std::cout << (unsigned) i->second.rule_list[p];
                if(p+1 == i->second.rule_count) std::cout << "\n";
                else
                {
                    std::cout << ',';
                    if(p%10 == 9)
                        std::cout << "\n              ";
                }
            }
            std::cout << "    }   };  }\n";
        }
        std::cout <<
            "}\n";
    }
private:
};

static GrammarDumper dumper;

%}

%pure_parser

%union {
    /* Note: Because bison's token type is an union or a simple type,
     *       anything that has constructors and destructors must be
     *       carried behind pointers here.
     */
    GrammarData::Rule*          r;
    GrammarData::FunctionType*  f;
    GrammarData::MatchedParams* p;
    GrammarData::ParamSpec*     a;

    mycomplex         num;
    unsigned                     index;
    FUNCTIONPARSERTYPES::OPCODE  opcode;
}

/* See documentation about syntax and token meanings in fpoptimizer.dat */
%token <num>       NUMERIC_CONSTANT
%token <index>     NAMEDHOLDER_TOKEN
%token <index>     RESTHOLDER_TOKEN
%token <index>     IMMEDHOLDER_TOKEN
%token <opcode>    BUILTIN_FUNC_NAME
%token <opcode>    OPCODE_TOKEN
%token <opcode>    UNARY_TRANSFORMATION
%token <index>     PARAM_CONSTRAINT
%token <index>     CONST_CONSTRAINT
%token NEWLINE

%token SUBST_OP_COLON /* ':' */
%token SUBST_OP_ARROW /* '->'  */

%type <r> substitution
%type <f> function function_match
%type <p> paramlist
%type <a> param
%type <index> param_constraints const_constraints rule_constraints

%%
    grammar:
      grammar substitution
      {
        grammar.AddRule(*$2);
        delete $2;
      }
    | grammar rule_constraints substitution
      {
        $3->SetSituationFlags($2);
        grammar.AddRule(*$3);
        delete $3;
      }
    | grammar NEWLINE
    | /* empty */
    ;

    rule_constraints:
      rule_constraints PARAM_CONSTRAINT
      {
        // Translate constraint flags to Rule Constraints.
        // They were lexed as param constraints, which
        // is why they have a different value.
        if($2 == Value_Logical)         // @L
          $$ = $1 | LogicalContextOnly;
        else if($2 == Value_NonInteger) // @F
          $$ = $1 | NotForIntegers;
        else if($2 == Value_IsInteger) // @I
          $$ = $1 | OnlyForIntegers;
        else if($2 == Constness_Const)
          $$ = $1 | OnlyForComplex;
        else
        {
          char msg[] = "Only @L, @F, @I, @C and @R rule constraints are allowed for now";
          yyerror(msg); YYERROR;
        }
      }
    | rule_constraints CONST_CONSTRAINT
      {
        if($2 == Modulo_Radians)
          $$ = $1 | NotForComplex;
        else
        {
          char msg[] = "Only @L, @F, @I, @C and @R rule constraints are allowed for now";
          yyerror(msg); YYERROR;
        }
      }
    |  /* empty */
      {
        $$ = 0;
      }
    ;

    substitution:
      function_match SUBST_OP_ARROW param NEWLINE
      /* Entire function is changed into the particular param */
      {
        $3->RecursivelySetDefaultParamMatchingType();

        $$ = new GrammarData::Rule(ProduceNewTree, *$1, $3);
        delete $1;
      }

    | function_match SUBST_OP_ARROW function NEWLINE
      /* Entire function changes, the param_notinv_list is rewritten */
      /* NOTE: "p x -> o y"  is a shortcut for "p x -> (o y)"  */
      {
        GrammarData::ParamSpec* p = new GrammarData::ParamSpec($3);
        p->RecursivelySetDefaultParamMatchingType();
        /*if(!$3->Params.EnsureNoRepeatedNamedHolders())
        {
            char msg[] = "The replacement function may not specify the same variable twice";
            yyerror(msg); YYERROR;
        }*/

        $$ = new GrammarData::Rule(ProduceNewTree, *$1, p);

        //std::cout << GrammarDumper().Dump(*new GrammarData::ParamSpec($3)) << "\n";
        delete $1;
      }

    | function_match SUBST_OP_COLON  paramlist NEWLINE
      /* The params provided are replaced with the new param_maybeinv_list */
      {
        /*if($1->Params.RestHolderIndex != 0)
        {
            char msg[] = "Restholder is not valid in the outermost function when ReplaceParams is used";
            yyerror(msg); YYERROR;
        }*/
        $3->RecursivelySetDefaultParamMatchingType();
        /*if(!$3->EnsureNoRepeatedNamedHolders())
        {
            char msg[] = "The replacement function may not specify the same variable twice";
            yyerror(msg); YYERROR;
        }*/

        $$ = new GrammarData::Rule(ReplaceParams, *$1, *$3);
        delete $1;
        delete $3;
      }
    ;

    function_match:
       function
       {
           if(!$1->Params.EnsureNoVariableCoverageParams_InPositionalParamLists())
           {
               char msg[] = "Restholders such as <1>, must not occur in bracketed param lists on the matching side";
               yyerror(msg); YYERROR;
           }
           $$ = $1;
       }
    ;

    function:
       OPCODE_TOKEN '[' paramlist ']'
       /* Match a function with opcode=opcode,
        * and the exact parameter list as specified
        */
       {
         $$ = new GrammarData::FunctionType($1, *$3);
         delete $3;
       }
    |  OPCODE_TOKEN '{' paramlist '}'
       /* Match a function with opcode=opcode,
        * and the exact parameter list in any order
        */
       {
         $$ = new GrammarData::FunctionType($1, *$3->SetType(SelectedParams));
         delete $3;
       }
    |  OPCODE_TOKEN paramlist
       /* Match a function with opcode=opcode and the given way of matching params */
       /* There may be more parameters, don't care about them */
       {
         $$ = new GrammarData::FunctionType($1, *$2->SetType(AnyParams));
         delete $2;
       }
    ;

    paramlist: /* left-recursive list of 0-n params with no delimiter */
        paramlist param    /* param */
        {
          $$ = $1->AddParam($2);
        }
      | paramlist RESTHOLDER_TOKEN /* a placeholder for all remaining params */
        {
          if($1->RestHolderIndex != 0)
          {
              char msg[] = "Illegal attempt to specify two restholders for the same param list";
              yyerror(msg); YYERROR;
          }
          $1->RestHolderIndex = $2;
          $$ = $1;
        }
      | /* empty */
        {
          $$ = new GrammarData::MatchedParams;
        }
    ;

    param:
       NUMERIC_CONSTANT const_constraints /* particular immed */
       {
         $$ = new GrammarData::ParamSpec($1, $2);
       }
    |  IMMEDHOLDER_TOKEN param_constraints  /* a placeholder for some immed */
       {
         $$ = new GrammarData::ParamSpec($1, GrammarData::ParamSpec::ParamHolderTag());
         $$->SetConstraint($2 | Constness_Const);
       }
    |  BUILTIN_FUNC_NAME '(' paramlist ')'  /* literal logarithm/sin/etc. of the provided immed-type params -- also sum/product/minimum/maximum */
       {
         /* Verify that $3 consists of constants */
         $$ = new GrammarData::ParamSpec($1, $3->GetParams() );
         if(!$$->VerifyIsConstant())
         {
             char msg[] = "Not constant";
             yyerror(msg); YYERROR;
         }
         delete $3;
       }
    |  NAMEDHOLDER_TOKEN param_constraints /* any expression, indicated by "x", "a" etc. */
       {
         $$ = new GrammarData::ParamSpec($1 + 2, GrammarData::ParamSpec::ParamHolderTag());
         $$->SetConstraint($2);
       }
    |  '(' function ')' param_constraints    /* a subtree */
       {
         $$ = new GrammarData::ParamSpec($2);
         $$->SetConstraint($4);
       }
    |  UNARY_TRANSFORMATION param   /* the negated/inverted literal value of the param */
       {
         /* Verify that $2 is constant */
         if(!$2->VerifyIsConstant())
         {
             char msg[] = "Not constant";
             yyerror(msg); YYERROR;
         }
         std::vector<GrammarData::ParamSpec*> tmp;
         tmp.push_back($2);
         $$ = new GrammarData::ParamSpec($1, tmp);
       }
    ;

    param_constraints: /* List of possible constraints to the given param, eg. odd,int,etc */
       param_constraints PARAM_CONSTRAINT
       {
         $$ = $1 | $2;
       }
    |  /* empty */
       {
         $$ = 0;
       }
    ;

    const_constraints: /* List of possible constraints to the given param */
       const_constraints CONST_CONSTRAINT
       {
         $$ = $1 | $2;
       }
    |  /* empty */
       {
         $$ = 0;
       }
    ;
%%

#ifndef FP_SUPPORT_OPTIMIZER
enum { cVar,cFetch,cPopNMov };
#endif

static void yyerror(const char* msg)
{
    std::cerr << msg << std::endl;
    for(;;)
    {
        int c = std::fgetc(stdin);
        if(c == EOF) break;
        std::fputc(c, stderr);
    }
    exit(1);
}

static int yylex(YYSTYPE* lval)
{
    int c = std::fgetc(stdin);
    switch(c)
    {
        case EOF: break;
        case '#':
            while(c != EOF && c != '\n') c = std::fgetc(stdin);
            return NEWLINE;
        case '\n':
        {
            c = std::fgetc(stdin);
            std::ungetc(c, stdin);
            if(c == '['
            || c == '$')
                return EOF;
            return NEWLINE;
        }
        case '+':
        {
            c = std::fgetc(stdin);
            std::ungetc(c, stdin);
            if(c == '(') { lval->opcode = FUNCTIONPARSERTYPES::cAdd; return BUILTIN_FUNC_NAME; }
            return '+';
        }
        case '*':
        {
            c = std::fgetc(stdin);
            std::ungetc(c, stdin);
            if(c == '(') { lval->opcode = FUNCTIONPARSERTYPES::cMul; return BUILTIN_FUNC_NAME; }
            return '*';
        }
        case '-':
        {
            int c2 = std::fgetc(stdin);
            if(c2 == '>') return SUBST_OP_ARROW;
            std::ungetc(c2, stdin);
            if(c2 >= '0' && c2 <= '9')
            {
                goto GotNumeric;
            }
            lval->opcode = FUNCTIONPARSERTYPES::cNeg;
            return UNARY_TRANSFORMATION;
        }
        case '/':
            lval->opcode = FUNCTIONPARSERTYPES::cInv;
            return UNARY_TRANSFORMATION;

        case '=':
        {
            int c2 = std::fgetc(stdin);
            std::ungetc(c2, stdin);
            return '=';
        }
        case '[': case '{':
        case ']': case '}':
        case '(':
        case ')':
            return c;
        case ' ':
        case '\t':
        case '\v':
        case '\r':
            return yylex(lval); // Counts as tail recursion, I hope
        case ':':
            return SUBST_OP_COLON;
        case '%': { lval->index = 0; return IMMEDHOLDER_TOKEN; }
        case '&': { lval->index = 1; return IMMEDHOLDER_TOKEN; }

        case '@':
        {
            int c2 = std::fgetc(stdin);
            switch(c2)
            {
                case 'E': { lval->index = Value_EvenInt; return PARAM_CONSTRAINT; }
                case 'O': { lval->index = Value_OddInt; return PARAM_CONSTRAINT; }
                case 'I': { lval->index = Value_IsInteger; return PARAM_CONSTRAINT; }
                case 'F': { lval->index = Value_NonInteger; return PARAM_CONSTRAINT; }
                case 'L': { lval->index = Value_Logical; return PARAM_CONSTRAINT; }
                case 'P': { lval->index = Sign_Positive; return PARAM_CONSTRAINT; }
                case 'N': { lval->index = Sign_Negative; return PARAM_CONSTRAINT; }
                case 'Q': { lval->index = Sign_NoIdea; return PARAM_CONSTRAINT; }
                case '1': { lval->index = Oneness_One; return PARAM_CONSTRAINT; }
                case 'M': { lval->index = Oneness_NotOne; return PARAM_CONSTRAINT; }
                case 'C': { lval->index = Constness_Const; return PARAM_CONSTRAINT; }
                case 'V': { lval->index = Constness_NotConst; return PARAM_CONSTRAINT; }
                case 'R': { lval->index = Modulo_Radians; return CONST_CONSTRAINT; }
            }
            std::ungetc(c2, stdin);
            return '@';
        }
        case '<':
        {
            lval->index  = 0;
            for(;;)
            {
                c = std::fgetc(stdin);
                if(c < '0' || c > '9') { std::ungetc(c, stdin); break; }
                lval->index = lval->index * 10 + (c-'0');
            }
            c = std::fgetc(stdin);
            if(c != '>') std::ungetc(c, stdin);
            return RESTHOLDER_TOKEN;
        }
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
        {
        GotNumeric:;
            std::string NumBuf;
            NumBuf += (char)c;
            bool had_comma = false;
            for(;;)
            {
                c = std::fgetc(stdin);
                if(c >= '0' && c <= '9')  { NumBuf += (char)c; continue; }
                if(c == '.' && !had_comma){ had_comma = true; NumBuf += (char)c; continue; }
                std::ungetc(c, stdin);
                break;
            }
            lval->num.real = std::strtod(NumBuf.c_str(), 0);
            lval->num.imag = 0.0;
            if(c == 'i')
            {
                std::fgetc(stdin);
                lval->num.imag = lval->num.real;
                lval->num.real = 0.0;
            }
            else if(c == '-' || c == '+')
            {
                NumBuf.clear();
                NumBuf += (char)c;
                bool had_comma = false;
                for(;;)
                {
                    c = std::fgetc(stdin);
                    if(c >= '0' && c <= '9')  { NumBuf += (char)c; continue; }
                    if(c == '.' && !had_comma){ had_comma = true; NumBuf += (char)c; continue; }
                    std::ungetc(c, stdin);
                    break;
                }
                if(c == 'i')
                {
                    lval->num.imag = std::strtod(NumBuf.c_str(), 0);
                    std::fgetc(stdin);
                }
            }
            return NUMERIC_CONSTANT;
        }
        case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
        case 'G': case 'H': case 'I': case 'J': case 'K': case 'L':
        case 'M': case 'N': case 'O': case 'P': case 'Q': case 'R':
        case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
        case 'Y': case 'Z': case '_':
        case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
        case 'g': case 'h': case 'i': case 'j': case 'k': case 'l':
        case 'm': case 'n': case 'o': case 'p': case 'q': case 'r':
        case 's': case 't': case 'u': case 'v': case 'w': case 'x':
        case 'y': case 'z':
        {
            std::string IdBuf;
            IdBuf += (char)c;
            for(;;)
            {
                c = std::fgetc(stdin);
                if((c >= '0' && c <= '9')
                || c == '_'
                || (c >= 'a' && c <= 'z')
                || (c >= 'A' && c <= 'Z')) { IdBuf += (char)c; continue; }
                std::ungetc(c, stdin);
                break;
            }
            lval->num.real = 0;
            lval->num.imag = 0;

            /* This code figures out if this is a named constant,
               an opcode, or a parse-time function name,
               or just an identifier
             */

            /* Detect named constants */
            if(IdBuf == "i")          { lval->num.imag = 1.0; return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_E") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_e<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_EI") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_einv<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_2E") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_twoe<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_2EI") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_twoeinv<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_RD") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_rad_to_deg<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_DR") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_deg_to_rad<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_PI") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_pi<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_PIHALF") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_pihalf<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_TWOPI") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_twopi<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_L2I") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_log2inv<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_L10I") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_log10inv<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_L2") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_log2<double>(); return NUMERIC_CONSTANT; }
            if(IdBuf == "CONSTANT_L10") { lval->num.real = FUNCTIONPARSERTYPES::fp_const_log10<double>(); return NUMERIC_CONSTANT; }

            /* Detect opcodes */
            if(IdBuf == "cAdd") { lval->opcode = FUNCTIONPARSERTYPES::cAdd; return OPCODE_TOKEN; }
            if(IdBuf == "cAnd") { lval->opcode = FUNCTIONPARSERTYPES::cAnd; return OPCODE_TOKEN; }
            if(IdBuf == "cMul") { lval->opcode = FUNCTIONPARSERTYPES::cMul; return OPCODE_TOKEN; }
            if(IdBuf == "cOr")  { lval->opcode = FUNCTIONPARSERTYPES::cOr; return OPCODE_TOKEN; }

            if(IdBuf == "cNeg") { lval->opcode = FUNCTIONPARSERTYPES::cNeg; return OPCODE_TOKEN; }
            if(IdBuf == "cSub") { lval->opcode = FUNCTIONPARSERTYPES::cSub; return OPCODE_TOKEN; }
            if(IdBuf == "cDiv") { lval->opcode = FUNCTIONPARSERTYPES::cDiv; return OPCODE_TOKEN; }
            if(IdBuf == "cMod") { lval->opcode = FUNCTIONPARSERTYPES::cMod; return OPCODE_TOKEN; }
            if(IdBuf == "cEqual") { lval->opcode = FUNCTIONPARSERTYPES::cEqual; return OPCODE_TOKEN; }
            if(IdBuf == "cNEqual") { lval->opcode = FUNCTIONPARSERTYPES::cNEqual; return OPCODE_TOKEN; }
            if(IdBuf == "cLess") { lval->opcode = FUNCTIONPARSERTYPES::cLess; return OPCODE_TOKEN; }
            if(IdBuf == "cLessOrEq") { lval->opcode = FUNCTIONPARSERTYPES::cLessOrEq; return OPCODE_TOKEN; }
            if(IdBuf == "cGreater") { lval->opcode = FUNCTIONPARSERTYPES::cGreater; return OPCODE_TOKEN; }
            if(IdBuf == "cGreaterOrEq") { lval->opcode = FUNCTIONPARSERTYPES::cGreaterOrEq; return OPCODE_TOKEN; }
            if(IdBuf == "cNot") { lval->opcode = FUNCTIONPARSERTYPES::cNot; return OPCODE_TOKEN; }
            if(IdBuf == "cNotNot") { lval->opcode = FUNCTIONPARSERTYPES::cNotNot; return OPCODE_TOKEN; }
            if(IdBuf == "cAbsNot") { lval->opcode = FUNCTIONPARSERTYPES::cAbsNot; return OPCODE_TOKEN; }
            if(IdBuf == "cAbsNotNot") { lval->opcode = FUNCTIONPARSERTYPES::cAbsNotNot; return OPCODE_TOKEN; }
            if(IdBuf == "cAbsAnd") { lval->opcode = FUNCTIONPARSERTYPES::cAbsAnd; return OPCODE_TOKEN; }
            if(IdBuf == "cAbsOr") { lval->opcode = FUNCTIONPARSERTYPES::cAbsOr; return OPCODE_TOKEN; }
            if(IdBuf == "cAbsIf") { lval->opcode = FUNCTIONPARSERTYPES::cAbsIf; return OPCODE_TOKEN; }
            if(IdBuf == "cDeg")  { lval->opcode = FUNCTIONPARSERTYPES::cDeg; return OPCODE_TOKEN; }
            if(IdBuf == "cRad")  { lval->opcode = FUNCTIONPARSERTYPES::cRad; return OPCODE_TOKEN; }
            if(IdBuf == "cInv")  { lval->opcode = FUNCTIONPARSERTYPES::cInv; return OPCODE_TOKEN; }
            if(IdBuf == "cSqr")  { lval->opcode = FUNCTIONPARSERTYPES::cSqr; return OPCODE_TOKEN; }
            if(IdBuf == "cRDiv") { lval->opcode = FUNCTIONPARSERTYPES::cRDiv; return OPCODE_TOKEN; }
            if(IdBuf == "cRSub") { lval->opcode = FUNCTIONPARSERTYPES::cRSub; return OPCODE_TOKEN; }
            if(IdBuf == "cRSqrt") { lval->opcode = FUNCTIONPARSERTYPES::cRSqrt; return OPCODE_TOKEN; }
#ifdef FP_SUPPORT_OPTIMIZER
            if(IdBuf == "cLog2by") { lval->opcode = FUNCTIONPARSERTYPES::cLog2by; return OPCODE_TOKEN; }
#else
            if(IdBuf == "cLog2by") { lval->opcode = FUNCTIONPARSERTYPES::cNop; return OPCODE_TOKEN; }
#endif

            /* Detect other function opcodes */
            if(IdBuf[0] == 'c' && std::isupper(IdBuf[1]))
            {
                // This has a chance of being an opcode token
                std::string opcodetoken = IdBuf.substr(1);
                opcodetoken[0] = (char) std::tolower( opcodetoken[0] );

                unsigned nameLength = readOpcode(opcodetoken.c_str());
                if(nameLength & 0x80000000U)
                {
                    lval->opcode = FUNCTIONPARSERTYPES::OPCODE(
                        (nameLength >> 16) & 0x7FFF );
                    return OPCODE_TOKEN;
                }
                std::cerr <<
                    "Warning: Unrecognized opcode '" << IdBuf << "' interpreted as cNop\n";
                lval->opcode = FUNCTIONPARSERTYPES::cNop;
                return OPCODE_TOKEN;
            }

            // If it is typed entirely in capitals, it has a chance of being
            // a group token
            if(true)
            {
                std::string grouptoken = IdBuf;
                for(size_t a=0; a<grouptoken.size(); ++a)
                {
                    if(std::islower(grouptoken[a])) goto NotAGroupToken;
                    grouptoken[a] = (char) std::tolower(grouptoken[a]);
                }
                if(1) // scope
                {
                    unsigned nameLength = readOpcode(grouptoken.c_str());
                    if(nameLength & 0x80000000U)
                    {
                        lval->opcode = FUNCTIONPARSERTYPES::OPCODE(
                            (nameLength >> 16) & 0x7FFF );
                        return BUILTIN_FUNC_NAME;
                    }
                    if(IdBuf == "MOD")
                    {
                        lval->opcode = FUNCTIONPARSERTYPES::cMod;
                        return BUILTIN_FUNC_NAME;
                    }
                    if(IdBuf == "DIV")
                    {
                        lval->opcode = FUNCTIONPARSERTYPES::cDiv;
                        return BUILTIN_FUNC_NAME;
                    }
                    if(IdBuf == "SUB")
                    {
                        lval->opcode = FUNCTIONPARSERTYPES::cSub;
                        return BUILTIN_FUNC_NAME;
                    }

                    std::cerr << "Warning: Unrecognized constant function '" << IdBuf
                              << "' interpreted as cNop\n";
                    lval->opcode = FUNCTIONPARSERTYPES::cNop;
                    return BUILTIN_FUNC_NAME;
                }
            NotAGroupToken:;
            }
            // Anything else is an identifier
            lval->index = dumper.ConvertNamedHolderNameIntoIndex(IdBuf);
            // std::cerr << "'" << IdBuf << "'' interpreted as PARAM\n";

            return NAMEDHOLDER_TOKEN;
        }
        default:
        {
            std::cerr << "Ignoring unidentifier character '" << char(c) << "'\n";
            return yylex(lval); // tail recursion
        }
    }
    return EOF;
}

unsigned GrammarData::ParamSpec::BuildDepMask()
{
    DepMask = 0;
    switch(Opcode)
    {
        case ParamHolder:
            DepMask |= 1 << Index;
            break;
        case SubFunction:
            DepMask = Func->Params.BuildDepMask();
            break;
        default: break;
    }
    return DepMask;
}

namespace FPoptimizer_Grammar
{
    template<typename Value_t>
    ParamSpec ParamSpec_Extract(unsigned paramlist, unsigned index)
    {
        unsigned plist_index = (paramlist >> (index*PARAM_INDEX_BITS))
                               % (1 << PARAM_INDEX_BITS);
        return plist[plist_index];
    }
    template ParamSpec ParamSpec_Extract<stdcomplex>(unsigned paramlist, unsigned index);
}

int main()
{
    std::map<std::string, GrammarData::Grammar> sections;

    std::string sectionname;

    for(;;)
    {
        grammar = GrammarData::Grammar();

        yyparse();

        grammar.BuildFinalDepMask();
        sections[sectionname] = grammar;

        int c = std::fgetc(stdin);
        if(c != '[')
        {
            std::ungetc(c, stdin);
            break;
        }

        sectionname.clear();
        for(;;)
        {
            c = std::fgetc(stdin);
            if(c == ']' || c == EOF) break;
            sectionname += (char)c;
        }
        std::cerr << "Parsing [" << sectionname << "]\n";
    }

    std::map<std::string, std::vector<std::string> > grammar_components;
    sectionname = "";
    for(;;)
    {
        int c = std::fgetc(stdin);
        if(c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
        if(c == '#')
            { do { c = std::fgetc(stdin); } while(!(c == '\n' || c == EOF));
              continue; }
        if(c == '$')
        {
            sectionname = "";
            for(;;)
            {
                c = std::fgetc(stdin);
                if(c == EOF) break;
                if(c == ' ' || c == '\t' || c == '\r' || c == '\n') break;
                if(c == ':') break;
                sectionname += char(c);
            }
            std::cerr << "Parsing $" << sectionname << "\n";
            continue;
        }
        if((c >= 'A' && c <= 'Z') || c == '_' || (c >= '0' && c <= '9'))
        {
            std::string componentname;
            for(;;)
            {
                if(c == EOF) break;
                if(c == ' ' || c == '\t' || c == '\r' || c == '\n') break;
                componentname += char(c);
                c = std::fgetc(stdin);
            }
            std::cerr << "- Has [" << componentname << "]\n";
            grammar_components[sectionname].push_back(componentname);
            //dumper.AddRulesFrom(sections[componentname]);
        }
        else break;
    }

    std::cout <<
        "/* This file is automatically generated. Do not edit... */\n"
        "#include \"../fpoptimizer/consts.hh\"\n"
        "#include \"fpconfig.hh\"\n"
        "#include \"extrasrc/fptypes.hh\"\n"
        "#include <algorithm>\n"
        "\n";

    std::vector<GrammarData::Grammar> components;
    for(std::map<std::string, std::vector<std::string> >::const_iterator
        i = grammar_components.begin();
        i != grammar_components.end();
        ++i)
    {
        for(size_t a=0; a<i->second.size(); ++a)
            components.push_back(sections[ i->second[a] ]);
    }
    dumper.RegisterGrammar(components);

    for(std::map<std::string, std::vector<std::string> >::const_iterator
        i = grammar_components.begin();
        i != grammar_components.end();
        ++i)
    {
        components.clear();
        for(size_t a=0; a<i->second.size(); ++a)
            components.push_back(sections[ i->second[a] ]);
        dumper.DumpGrammar(i->first, components);
    }
    dumper.Flush();

    unsigned mask = (1 << PARAM_INDEX_BITS)-1;
    const unsigned p_begin = 0;
    const unsigned n_begin = p_begin + dumper.collection.plist_p.size();
    const unsigned s_begin = n_begin + dumper.collection.plist_n.size();
    std::cout <<
        "namespace FPoptimizer_Grammar\n"
        "{\n"
        "    template<typename Value_t>\n"
        "    ParamSpec ParamSpec_Extract(unsigned paramlist, unsigned index)\n"
        "    {\n"
        "        index = (paramlist >> (index * " << PARAM_INDEX_BITS << ")) & " << mask << " /* % (1 << " << PARAM_INDEX_BITS << ") */;\n"
        "        if(index >= " << s_begin << ")\n"
        "            return ParamSpec(SubFunction,(const void*)&plist_s[index-" << s_begin << "]);\n"
        "        if(index >= " << n_begin << ")\n"
        "            return ParamSpec(NumConstant,(const void*)&plist_n_container<Value_t>::plist_n[index-" << n_begin << "]);\n"
        "        return ParamSpec(ParamHolder,(const void*)&plist_p[index"/*"-" << p_begin << */"]);\n"
        "    }\n"
        "}\n"
        "/* BEGIN_EXPLICIT_INSTANTATION */\n"
        "#include \"instantiate.hh\"\n"
        "namespace FPoptimizer_Grammar\n"
        "{\n"
        "#define FP_INSTANTIATE(type) \\\n"
        "    template ParamSpec ParamSpec_Extract<type>(unsigned paramlist, unsigned index);\n"
        "    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)\n"
        "#undef FP_INSTANTIATE\n"
        "}\n"
        "/* END_EXPLICIT_INSTANTATION */\n";

    return 0;
}
