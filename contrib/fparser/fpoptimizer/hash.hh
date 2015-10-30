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
