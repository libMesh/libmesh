// the space is needed to preserve the include in the postprocessor
# include "instantiate.hh"

namespace FPoptimizer_CodeTree
{

#define FP_INSTANTIATE(type) \
    template void CodeTree<type>::SynthesizeByteCode( \
        std::vector<unsigned>& ByteCode, \
        std::vector<type>&   Immed, \
        size_t& stacktop_max);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE

#define FP_INSTANTIATE(type) \
    template \
    void CodeTree<type>::GenerateFrom( \
        const FunctionParserBase<type>::Data& fpdata, \
        bool keep_powi);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE

#define FP_INSTANTIATE(type) \
    template class CodeTree<type>; \
    template struct CodeTreeData<type>;
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE

#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
#define FP_INSTANTIATE(type) \
    template void DumpTreeWithIndent<type>(const CodeTree<type>&, std::ostream&, const std::string&);
    FPOPTIMIZER_EXPLICITLY_INSTANTIATE(FP_INSTANTIATE)
#undef FP_INSTANTIATE
#endif

}
