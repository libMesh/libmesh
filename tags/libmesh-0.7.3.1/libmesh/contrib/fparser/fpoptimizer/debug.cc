#include "codetree.hh"
#include "opcodename.hh"

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

/* BEGIN_EXPLICIT_INSTANTATION */
#include "instantiate.hh"
namespace FPoptimizer_CodeTree
{
#define FP_INSTANTIATE(type) \
    template void DumpHashes(const CodeTree<type>& tree, std::ostream& o); \
    template void DumpTree(const CodeTree<type>& tree, std::ostream& o); \
    template void DumpTreeWithIndent(const CodeTree<type>& tree, std::ostream& o, const std::string& indent);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
}
/* END_EXPLICIT_INSTANTATION */

#endif
