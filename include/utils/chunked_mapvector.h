// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef LIBMESH_CHUNKED_MAPVECTOR_H
#define LIBMESH_CHUNKED_MAPVECTOR_H

// C++ Includes   -----------------------------------
#include <map>

namespace libMesh
{

/**
 * This \p chunked_mapvector templated class is intended to provide
 * the asymptotic performance characteristics of a std::map with an
 * interface more closely resembling that of a std::vector, for use
 * with DistributedMesh.
 *
 * The intermediate array "chunks" give better constants on
 * performance, for the typical sparsity structure we see on meshes in
 * practice, where large swaths of ids are contiguous.
 *
 * \author  Roy H. Stogner
 */

template <typename Val, typename index_t=unsigned int, unsigned int N=16>
class chunked_mapvector : public std::map<index_t, std::array<Val,N>>
{
public:
  typedef std::map<index_t, std::array<Val,N>> maptype;

  typedef unsigned int iter_t; // Only has to hold 0 through N

  template <typename MapIter>
  class veclike_iterator_base
  {
  public:
    veclike_iterator_base(const MapIter & i,
                     const iter_t idx)
      : it(i), array_index(idx) {}

    veclike_iterator_base(const veclike_iterator_base & i) = default;

    template <typename T>
    veclike_iterator_base(const veclike_iterator_base<T> & i)
      : it(i.it), array_index(i.idx) {}

    veclike_iterator_base & operator++() {
      ++array_index;
      if (array_index >= N)
        {
          array_index = 0;
          ++it;
        }
      return *this;
    }

    veclike_iterator_base operator++(int) {
      veclike_iterator i = *this;
      ++(*this);
      return i;
    }

    bool operator==(const veclike_iterator_base & other) const {
      return (it == other.it && array_index == other.array_index);
    }

    bool operator!=(const veclike_iterator_base & other) const {
      return (it != other.it || array_index != other.array_index);
    }

    index_t index() const { return this->it->first * N + this->array_index; }

  private:
    friend class chunked_mapvector;

    MapIter it;

    iter_t array_index;
  };

  class veclike_iterator :
    public veclike_iterator_base<typename maptype::iterator>
  {
  public:
    veclike_iterator(const typename maptype::iterator & i,
                     const iter_t idx)
      : veclike_iterator_base<typename maptype::iterator>(i,idx) {}

    Val & operator*() const { return (this->it->second)[this->array_index]; }

    Val * operator->() const { return &((this->it->second)[this->array_index]); }
  };


  class const_veclike_iterator :
    public veclike_iterator_base<typename maptype::const_iterator>
  {
  public:
    const_veclike_iterator(const typename maptype::const_iterator & i,
                           const iter_t idx)
      : veclike_iterator_base<typename maptype::const_iterator>(i,idx) {}

    const Val & operator*() const { return (this->it->second)[this->array_index]; }

    const Val * operator->() const { return &((this->it->second)[this->array_index]); }
  };

  class const_reverse_veclike_iterator
  {
  public:
    const_reverse_veclike_iterator(const typename maptype::const_reverse_iterator & i,
                                   const iter_t idx)
      : it(i), array_index(idx) {}

    const_reverse_veclike_iterator(const const_reverse_veclike_iterator & i) = default;

    const_reverse_veclike_iterator & operator++() {
      if (array_index == 0)
        {
          array_index = N-1;
          ++it;
        }
      else
        --array_index;
      return *this;
    }

    const_reverse_veclike_iterator operator++(int) {
      veclike_iterator i = *this;
      ++(*this);
      return i;
    }

    const Val & operator*() const { return (this->it->second)[this->array_index]; }

    const Val * operator->() const { return &((this->it->second)[this->array_index]); }

    index_t index() const { return this->it->first * N + this->array_index; }

    bool operator==(const const_reverse_veclike_iterator & other) const {
      return (it == other.it && array_index == other.array_index);
    }

    bool operator!=(const const_reverse_veclike_iterator & other) const {
      return (it != other.it || array_index != other.array_index);
    }

  private:
    typename maptype::const_reverse_iterator it;

    iter_t array_index;
  };

  veclike_iterator find (const index_t & k)
  {
    auto sub_it = maptype::find(k/N);
    if (sub_it == maptype::end())
      return veclike_iterator(sub_it, 0);

    veclike_iterator it {sub_it, iter_t(k%N)};
    if (*it)
      return it;
    return this->end();
  }

  const_veclike_iterator find (const index_t & k) const
  {
    auto sub_it = maptype::find(k/N);
    if (sub_it == maptype::end())
      return const_veclike_iterator(sub_it, 0);

    const_veclike_iterator it {sub_it, iter_t(k%N)};
    if (*it)
      return it;
    return this->end();
  }

  Val & operator[] (const index_t & k)
  {
    return maptype::operator[](k/N)[k%N];
  }

  Val operator[] (const index_t & k) const
  {
    typename maptype::const_iterator it = maptype::find(k/N);
    return it == this->end().it? Val() : (it->second)[k%N];
  }

  void erase(index_t i) {
    typename maptype::iterator it = maptype::find(i/N);
    if (it == maptype::end())
      return;

    (it->second)[i%N] = Val();
    for (auto v : it->second)
      if (v != Val())
        return;

    maptype::erase(it);
  }

  veclike_iterator erase(const veclike_iterator & pos) {
    if (pos.it == maptype::end())
      return pos;
    *pos = Val();

    veclike_iterator newpos = pos;
    do {
      ++newpos;
    } while (newpos.it == pos.it &&
             *newpos == Val());

    for (auto v : pos.it->second)
      if (v != Val())
        return newpos;

    maptype::erase(pos.it);

    return newpos;
  }

  veclike_iterator begin() {
    return veclike_iterator(maptype::begin(), 0);
  }

  const_veclike_iterator begin() const {
    return const_veclike_iterator(maptype::begin(), 0);
  }

  veclike_iterator end() {
    return veclike_iterator(maptype::end(), 0);
  }

  const_veclike_iterator end() const {
    return const_veclike_iterator(maptype::end(), 0);
  }

  const_reverse_veclike_iterator rbegin() const {
    return const_reverse_veclike_iterator(maptype::rbegin(), N-1);
  }

  const_reverse_veclike_iterator rend() const {
    return const_reverse_veclike_iterator(maptype::rend(), N-1);
  }
};

} // namespace libMesh

#endif // LIBMESH_CHUNKED_MAPVECTOR_H
