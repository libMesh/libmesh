// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_VECTORMAP_H
#define LIBMESH_VECTORMAP_H

// C++ Includes   -----------------------------------
#include <vector>
#include <algorithm>

// libMesh includes
#include "libmesh/libmesh_common.h"


namespace libMesh
{

  /**
   * This \p vectormap templated class is intended to provide the
   * performance characteristics of a sorted std::vector with an
   * interface more closely resembling that of a std::map, for use in
   * particular when memory is tight.
   *
   * This class is limited in its applicability.  The typical use case is:
   *
   \verbatim
   vectormap<KeyType,ValType> vmap;
   for ( ; ;)
   vmap.insert (std::make_pair(key,val));

   val1 = vmap[key1];
   val2 = vmap[key2];
   \endverbatim
   *
   * Note in particular the two-phase usage.  It is not advised to do
   * intermediate insert/lookups, as each time an insertion is done the
   * sort is invalidated, and must be reperformed. Additionally, during
   * the insersion phase, there is no accounting for duplicate keys.
   * Each insertion is accepted regardless of key value, and then
   * duplicate keys are removed later.
   *
   * \author  Benjamin S. Kirk
   */
  template <typename Key, typename Tp>
  class vectormap : public std::vector<std::pair<Key, Tp> >
  {

  public:

    typedef Key                     key_type;
    typedef Tp                      mapped_type;
    typedef std::pair<Key, Tp>      value_type;
    typedef std::vector<value_type> vector_type;

  private:

    /**
     * Strict weak ordering, based solely on first element in a pair.
     */
    struct FirstOrder
    {
      bool operator()(const value_type &lhs,
		      const value_type &rhs) const
      { return lhs.first < rhs.first; }
    };

    /**
     * Equality comparison, based solely on first element in a pair.
     */
    struct FirstCompare
    {
      bool operator()(const value_type &lhs,
		      const value_type &rhs) const
      { return lhs.first == rhs.first; }
    };

  public:

    /**
     * Default constructor.  Initializes sorted member to false.
     */
    vectormap() :
      _sorted(false)
    {}

    /**
     * Copy constructor.
     */
    vectormap(const vectormap<Key,Tp> &other) :
      std::vector<std::pair<Key, Tp> > (other),
      _sorted(other._sorted)
    {}

    /**
     * Inserts \p x into the vectormap.
     */
    void insert (const value_type &x)
    {
      _sorted = false;
      this->push_back(x);
    }

    /**
     * Sort & unique the vectormap, preparing for use.
     */
    void sort()
    {
      FirstOrder   order;
      FirstCompare comp;

      std::sort (this->begin(), this->end(), order);

      this->erase (std::unique (this->begin(), this->end(), comp), this->end());

      _sorted = true;
    }

    /**
     * @returns the value corresponding to \p key
     */
    const Tp & operator[](const key_type &key) const
    {
      if (!_sorted)
	const_cast<vectormap<Key, Tp>*>(this)->sort();

      libmesh_assert (_sorted);

      value_type to_find;
      to_find.first = key;

      FirstOrder order;

      typename vectormap<Key,Tp>::const_iterator
	lower_bound = std::lower_bound (this->begin(), this->end(), to_find, order);

      libmesh_assert (lower_bound != this->end());
      libmesh_assert_equal_to (lower_bound->first, key);

      return lower_bound->second;
    }

    /**
     * *returns the number of occurances of \p key.  For a map-like object, this should
     * be 1 or 0.
     */
    unsigned int count (const key_type &key) const
    {
      if (!_sorted)
	const_cast<vectormap<Key, Tp>*>(this)->sort();

      libmesh_assert (_sorted);

      value_type to_find;
      to_find.first = key;

      FirstOrder order;

      std::pair<typename vectormap<Key,Tp>::const_iterator,
		typename vectormap<Key,Tp>::const_iterator>
	bounds = std::equal_range (this->begin(), this->end(), to_find, order);

      return std::distance (bounds.first, bounds.second);
    }


  private:

    bool _sorted;
  };

} // namespace libMesh

#endif // LIBMESH_VECTORMAP_H
