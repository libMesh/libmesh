#include <string>

class CPPcompressor
{
public:
    /* Compress the given C++ source code. Warning: Function not re-entrant. */
    std::string Compress(const std::string& input);
    std::string Compress(const std::string& input,
        const std::string& macro_prefix_chars);
};
