
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 1 "util/tree_grammar_parser.y"

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



/* Line 189 of yacc.c  */
#line 1287 "util/tree_grammar_parser.cc"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUMERIC_CONSTANT = 258,
     NAMEDHOLDER_TOKEN = 259,
     RESTHOLDER_TOKEN = 260,
     IMMEDHOLDER_TOKEN = 261,
     BUILTIN_FUNC_NAME = 262,
     OPCODE_TOKEN = 263,
     UNARY_TRANSFORMATION = 264,
     PARAM_CONSTRAINT = 265,
     CONST_CONSTRAINT = 266,
     NEWLINE = 267,
     SUBST_OP_COLON = 268,
     SUBST_OP_ARROW = 269
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 1216 "util/tree_grammar_parser.y"

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



/* Line 214 of yacc.c  */
#line 1354 "util/tree_grammar_parser.cc"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 1366 "util/tree_grammar_parser.cc"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   86

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  21
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  10
/* YYNRULES -- Number of rules.  */
#define YYNRULES  28
/* YYNRULES -- Number of states.  */
#define YYNSTATES  47

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   269

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      19,    20,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    15,     2,    16,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    17,     2,    18,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     6,    10,    13,    14,    17,    20,    21,
      26,    31,    36,    38,    43,    48,    51,    54,    57,    58,
      61,    64,    69,    72,    77,    80,    83,    84,    87
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      22,     0,    -1,    22,    24,    -1,    22,    23,    24,    -1,
      22,    12,    -1,    -1,    23,    10,    -1,    23,    11,    -1,
      -1,    25,    14,    28,    12,    -1,    25,    14,    26,    12,
      -1,    25,    13,    27,    12,    -1,    26,    -1,     8,    15,
      27,    16,    -1,     8,    17,    27,    18,    -1,     8,    27,
      -1,    27,    28,    -1,    27,     5,    -1,    -1,     3,    30,
      -1,     6,    29,    -1,     7,    19,    27,    20,    -1,     4,
      29,    -1,    19,    26,    20,    29,    -1,     9,    28,    -1,
      29,    10,    -1,    -1,    30,    11,    -1,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,  1254,  1254,  1259,  1265,  1266,  1270,  1289,  1300,  1306,
    1315,  1333,  1355,  1367,  1375,  1383,  1393,  1397,  1408,  1414,
    1418,  1423,  1434,  1439,  1444,  1459,  1464,  1470,  1475
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUMERIC_CONSTANT", "NAMEDHOLDER_TOKEN",
  "RESTHOLDER_TOKEN", "IMMEDHOLDER_TOKEN", "BUILTIN_FUNC_NAME",
  "OPCODE_TOKEN", "UNARY_TRANSFORMATION", "PARAM_CONSTRAINT",
  "CONST_CONSTRAINT", "NEWLINE", "SUBST_OP_COLON", "SUBST_OP_ARROW", "'['",
  "']'", "'{'", "'}'", "'('", "')'", "$accept", "grammar",
  "rule_constraints", "substitution", "function_match", "function",
  "paramlist", "param", "param_constraints", "const_constraints", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,    91,    93,   123,   125,    40,
      41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    21,    22,    22,    22,    22,    23,    23,    23,    24,
      24,    24,    25,    26,    26,    26,    27,    27,    27,    28,
      28,    28,    28,    28,    28,    29,    29,    30,    30
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     3,     2,     0,     2,     2,     0,     4,
       4,     4,     1,     4,     4,     2,     2,     2,     0,     2,
       2,     4,     2,     4,     2,     2,     0,     2,     0
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       5,     8,     1,    18,     4,     0,     2,     0,    12,    18,
      18,    15,     6,     7,     3,    18,     0,     0,     0,    28,
      26,    17,    26,     0,     0,     0,    16,     0,     0,     0,
      13,    14,    19,    22,    20,    18,    24,     0,    11,    10,
       9,    27,    25,     0,    26,    21,    23
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,     5,     6,     7,     8,    11,    26,    33,    32
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -21
static const yytype_int8 yypact[] =
{
     -21,    74,   -21,   -12,   -21,    23,   -21,    63,   -21,   -21,
     -21,    52,   -21,   -21,   -21,   -21,    59,    16,    33,   -21,
     -21,   -21,   -21,   -15,    66,    -1,   -21,    41,    -2,    17,
     -21,   -21,    19,    31,    31,   -21,   -21,    29,   -21,   -21,
     -21,   -21,   -21,     8,   -21,   -21,    31
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -21,   -21,   -21,    38,   -21,    -7,    -9,    -8,   -20,   -21
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      17,    18,    34,     9,    35,    10,    27,     3,    29,    28,
      39,    19,    20,    21,    22,    23,    36,    24,    37,    19,
      20,    21,    22,    23,    46,    24,    43,    25,    45,    40,
      41,     3,    30,    12,    13,    25,    19,    20,    21,    22,
      23,    42,    24,    14,    19,    20,    21,    22,    23,    44,
      24,    31,    25,    38,     0,    19,    20,    21,    22,    23,
      25,    24,    19,    20,     0,    22,    23,     3,    24,    19,
      20,    25,    22,    23,     2,    24,    15,    16,    25,     0,
       0,     0,     3,     0,     0,    25,     4
};

static const yytype_int8 yycheck[] =
{
       9,    10,    22,    15,    19,    17,    15,     8,    16,    16,
      12,     3,     4,     5,     6,     7,    24,     9,    25,     3,
       4,     5,     6,     7,    44,     9,    35,    19,    20,    12,
      11,     8,    16,    10,    11,    19,     3,     4,     5,     6,
       7,    10,     9,     5,     3,     4,     5,     6,     7,    20,
       9,    18,    19,    12,    -1,     3,     4,     5,     6,     7,
      19,     9,     3,     4,    -1,     6,     7,     8,     9,     3,
       4,    19,     6,     7,     0,     9,    13,    14,    19,    -1,
      -1,    -1,     8,    -1,    -1,    19,    12
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    22,     0,     8,    12,    23,    24,    25,    26,    15,
      17,    27,    10,    11,    24,    13,    14,    27,    27,     3,
       4,     5,     6,     7,     9,    19,    28,    27,    26,    28,
      16,    18,    30,    29,    29,    19,    28,    26,    12,    12,
      12,    11,    10,    27,    20,    20,    29
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */





/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 1255 "util/tree_grammar_parser.y"
    {
        grammar.AddRule(*(yyvsp[(2) - (2)].r));
        delete (yyvsp[(2) - (2)].r);
      ;}
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 1260 "util/tree_grammar_parser.y"
    {
        (yyvsp[(3) - (3)].r)->SetSituationFlags((yyvsp[(2) - (3)].index));
        grammar.AddRule(*(yyvsp[(3) - (3)].r));
        delete (yyvsp[(3) - (3)].r);
      ;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 1271 "util/tree_grammar_parser.y"
    {
        // Translate constraint flags to Rule Constraints.
        // They were lexed as param constraints, which
        // is why they have a different value.
        if((yyvsp[(2) - (2)].index) == Value_Logical)         // @L
          (yyval.index) = (yyvsp[(1) - (2)].index) | LogicalContextOnly;
        else if((yyvsp[(2) - (2)].index) == Value_NonInteger) // @F
          (yyval.index) = (yyvsp[(1) - (2)].index) | NotForIntegers;
        else if((yyvsp[(2) - (2)].index) == Value_IsInteger) // @I
          (yyval.index) = (yyvsp[(1) - (2)].index) | OnlyForIntegers;
        else if((yyvsp[(2) - (2)].index) == Constness_Const)
          (yyval.index) = (yyvsp[(1) - (2)].index) | OnlyForComplex;
        else
        {
          char msg[] = "Only @L, @F, @I, @C and @R rule constraints are allowed for now";
          yyerror(msg); YYERROR;
        }
      ;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 1290 "util/tree_grammar_parser.y"
    {
        if((yyvsp[(2) - (2)].index) == Modulo_Radians)
          (yyval.index) = (yyvsp[(1) - (2)].index) | NotForComplex;
        else
        {
          char msg[] = "Only @L, @F, @I, @C and @R rule constraints are allowed for now";
          yyerror(msg); YYERROR;
        }
      ;}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 1300 "util/tree_grammar_parser.y"
    {
        (yyval.index) = 0;
      ;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 1308 "util/tree_grammar_parser.y"
    {
        (yyvsp[(3) - (4)].a)->RecursivelySetDefaultParamMatchingType();

        (yyval.r) = new GrammarData::Rule(ProduceNewTree, *(yyvsp[(1) - (4)].f), (yyvsp[(3) - (4)].a));
        delete (yyvsp[(1) - (4)].f);
      ;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 1318 "util/tree_grammar_parser.y"
    {
        GrammarData::ParamSpec* p = new GrammarData::ParamSpec((yyvsp[(3) - (4)].f));
        p->RecursivelySetDefaultParamMatchingType();
        /*if(!$3->Params.EnsureNoRepeatedNamedHolders())
        {
            char msg[] = "The replacement function may not specify the same variable twice";
            yyerror(msg); YYERROR;
        }*/

        (yyval.r) = new GrammarData::Rule(ProduceNewTree, *(yyvsp[(1) - (4)].f), p);

        //std::cout << GrammarDumper().Dump(*new GrammarData::ParamSpec($3)) << "\n";
        delete (yyvsp[(1) - (4)].f);
      ;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 1335 "util/tree_grammar_parser.y"
    {
        /*if($1->Params.RestHolderIndex != 0)
        {
            char msg[] = "Restholder is not valid in the outermost function when ReplaceParams is used";
            yyerror(msg); YYERROR;
        }*/
        (yyvsp[(3) - (4)].p)->RecursivelySetDefaultParamMatchingType();
        /*if(!$3->EnsureNoRepeatedNamedHolders())
        {
            char msg[] = "The replacement function may not specify the same variable twice";
            yyerror(msg); YYERROR;
        }*/

        (yyval.r) = new GrammarData::Rule(ReplaceParams, *(yyvsp[(1) - (4)].f), *(yyvsp[(3) - (4)].p));
        delete (yyvsp[(1) - (4)].f);
        delete (yyvsp[(3) - (4)].p);
      ;}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 1356 "util/tree_grammar_parser.y"
    {
           if(!(yyvsp[(1) - (1)].f)->Params.EnsureNoVariableCoverageParams_InPositionalParamLists())
           {
               char msg[] = "Restholders such as <1>, must not occur in bracketed param lists on the matching side";
               yyerror(msg); YYERROR;
           }
           (yyval.f) = (yyvsp[(1) - (1)].f);
       ;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 1371 "util/tree_grammar_parser.y"
    {
         (yyval.f) = new GrammarData::FunctionType((yyvsp[(1) - (4)].opcode), *(yyvsp[(3) - (4)].p));
         delete (yyvsp[(3) - (4)].p);
       ;}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 1379 "util/tree_grammar_parser.y"
    {
         (yyval.f) = new GrammarData::FunctionType((yyvsp[(1) - (4)].opcode), *(yyvsp[(3) - (4)].p)->SetType(SelectedParams));
         delete (yyvsp[(3) - (4)].p);
       ;}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 1386 "util/tree_grammar_parser.y"
    {
         (yyval.f) = new GrammarData::FunctionType((yyvsp[(1) - (2)].opcode), *(yyvsp[(2) - (2)].p)->SetType(AnyParams));
         delete (yyvsp[(2) - (2)].p);
       ;}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 1394 "util/tree_grammar_parser.y"
    {
          (yyval.p) = (yyvsp[(1) - (2)].p)->AddParam((yyvsp[(2) - (2)].a));
        ;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 1398 "util/tree_grammar_parser.y"
    {
          if((yyvsp[(1) - (2)].p)->RestHolderIndex != 0)
          {
              char msg[] = "Illegal attempt to specify two restholders for the same param list";
              yyerror(msg); YYERROR;
          }
          (yyvsp[(1) - (2)].p)->RestHolderIndex = (yyvsp[(2) - (2)].index);
          (yyval.p) = (yyvsp[(1) - (2)].p);
        ;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 1408 "util/tree_grammar_parser.y"
    {
          (yyval.p) = new GrammarData::MatchedParams;
        ;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 1415 "util/tree_grammar_parser.y"
    {
         (yyval.a) = new GrammarData::ParamSpec((yyvsp[(1) - (2)].num), (yyvsp[(2) - (2)].index));
       ;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 1419 "util/tree_grammar_parser.y"
    {
         (yyval.a) = new GrammarData::ParamSpec((yyvsp[(1) - (2)].index), GrammarData::ParamSpec::ParamHolderTag());
         (yyval.a)->SetConstraint((yyvsp[(2) - (2)].index) | Constness_Const);
       ;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 1424 "util/tree_grammar_parser.y"
    {
         /* Verify that $3 consists of constants */
         (yyval.a) = new GrammarData::ParamSpec((yyvsp[(1) - (4)].opcode), (yyvsp[(3) - (4)].p)->GetParams() );
         if(!(yyval.a)->VerifyIsConstant())
         {
             char msg[] = "Not constant";
             yyerror(msg); YYERROR;
         }
         delete (yyvsp[(3) - (4)].p);
       ;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 1435 "util/tree_grammar_parser.y"
    {
         (yyval.a) = new GrammarData::ParamSpec((yyvsp[(1) - (2)].index) + 2, GrammarData::ParamSpec::ParamHolderTag());
         (yyval.a)->SetConstraint((yyvsp[(2) - (2)].index));
       ;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 1440 "util/tree_grammar_parser.y"
    {
         (yyval.a) = new GrammarData::ParamSpec((yyvsp[(2) - (4)].f));
         (yyval.a)->SetConstraint((yyvsp[(4) - (4)].index));
       ;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 1445 "util/tree_grammar_parser.y"
    {
         /* Verify that $2 is constant */
         if(!(yyvsp[(2) - (2)].a)->VerifyIsConstant())
         {
             char msg[] = "Not constant";
             yyerror(msg); YYERROR;
         }
         std::vector<GrammarData::ParamSpec*> tmp;
         tmp.push_back((yyvsp[(2) - (2)].a));
         (yyval.a) = new GrammarData::ParamSpec((yyvsp[(1) - (2)].opcode), tmp);
       ;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 1460 "util/tree_grammar_parser.y"
    {
         (yyval.index) = (yyvsp[(1) - (2)].index) | (yyvsp[(2) - (2)].index);
       ;}
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 1464 "util/tree_grammar_parser.y"
    {
         (yyval.index) = 0;
       ;}
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 1471 "util/tree_grammar_parser.y"
    {
         (yyval.index) = (yyvsp[(1) - (2)].index) | (yyvsp[(2) - (2)].index);
       ;}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 1475 "util/tree_grammar_parser.y"
    {
         (yyval.index) = 0;
       ;}
    break;



/* Line 1455 of yacc.c  */
#line 2903 "util/tree_grammar_parser.cc"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 1479 "util/tree_grammar_parser.y"


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

