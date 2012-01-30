#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH

#include "codetree.hh"

namespace FPoptimizer_CodeTree
{
    template<typename Value_t>
    void ConstantFolding(CodeTree<Value_t>& tree);
}

#endif
