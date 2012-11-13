#include <utility>

struct Compare2ndRev
{
    template<typename T>
    inline bool operator() (const T& a, const T& b) const
    {
        return a.second > b.second;
    }
};

struct Compare1st
{
    template<typename T1, typename T2>
    inline bool operator() (const std::pair<T1,T2>& a,
                            const std::pair<T1,T2>& b) const
    {
        return a.first < b.first;
    }

    template<typename T1, typename T2>
    inline bool operator() (const std::pair<T1,T2>& a, T1 b) const
    {
        return a.first < b;
    }

    template<typename T1, typename T2>
    inline bool operator() (T1 a, const std::pair<T1,T2>& b) const
    {
        return a < b.first;
    }
};
