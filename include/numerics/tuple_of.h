#ifndef LIBMESH_TUPLE_OF_H
#define LIBMESH_TUPLE_OF_H

#include <cstddef> // size_t
#include <tuple>

namespace libMesh
{

// Recursive tuple scheme
template <std::size_t Index, typename T>
struct tuple_n
{
  template< typename...Args> using type = typename tuple_n<Index-1, T>::template type<T, Args...>;
};

template <typename T>
struct tuple_n<0, T>
{
  template<typename...Args> using type = std::tuple<Args...>;
};
template <std::size_t Index, typename T>  using tuple_of = typename tuple_n<Index,T>::template type<>;

}

#endif
