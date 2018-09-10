#ifndef LIBMESH_TUPLE_OF_H
#define LIBMESH_TUPLE_OF_H

#include <tuple>

namespace libMesh
{

// Recursive tuple scheme
template <size_t I, typename T>
struct tuple_n
{
  template< typename...Args> using type = typename tuple_n<I-1, T>::template type<T, Args...>;
};

template <typename T>
struct tuple_n<0, T>
{
  template<typename...Args> using type = std::tuple<Args...>;
};
template <size_t I, typename T>  using tuple_of = typename tuple_n<I,T>::template type<>;

}

#endif
