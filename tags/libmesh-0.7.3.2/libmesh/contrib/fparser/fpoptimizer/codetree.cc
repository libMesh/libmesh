#include <list>
#include <algorithm>

#include "rangeestimation.hh"
#include "optimize.hh" // for DEBUG_SUBSTITUTIONS
#include "codetree.hh"
#include "consts.hh"

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

#ifdef __GXX_EXPERIMENTAL_CXX0X__
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

#ifdef __GXX_EXPERIMENTAL_CXX0X__
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
    #ifdef __GXX_EXPERIMENTAL_CXX0X__
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

#ifdef __GXX_EXPERIMENTAL_CXX0X__
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

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template class CodeTree<type>; \
    template struct CodeTreeData<type>;
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif

