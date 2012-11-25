#include "fpconfig.hh"
#include "fparser.hh"
#include "extrasrc/fptypes.hh"

#ifdef FP_SUPPORT_OPTIMIZER

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <cmath>

#include <memory> /* for auto_ptr */

#include "grammar.hh"
#include "optimize.hh"
#include "rangeestimation.hh"
#include "consts.hh"

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
                    if(!&*start_at)
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
        if(&*start_at)
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
            if(&*(*position)[a].start_at) // is there another try?
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
                if(&*start_at)
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
                    if(&*(*position)[a].start_at) // is there another try?
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
                if(&*start_at)
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
                    if(&*(*position)[a].start_at) // is there another try?
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


/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
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
/* END_EXPLICIT_INSTANTATION */

#endif
