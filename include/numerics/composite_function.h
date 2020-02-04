// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_COMPOSITE_FUNCTION_H
#define LIBMESH_COMPOSITE_FUNCTION_H

// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/function_base.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh.h"
#include "libmesh/point.h"

// C++ includes
#include <algorithm>
#include <utility>
#include <vector>

namespace libMesh
{

/**
 * \brief A function that returns a vector whose components are defined by
 * multiple functions.
 *
 * A function which is defined by composing the result of different functions
 * into a single vector.  All overridden virtual functions are documented in
 * function_base.h.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief Function which is a function of another function.
 */
template <typename Output=Number>
class CompositeFunction : public FunctionBase<Output>
{
public:
  typedef typename GetScalarType<Output>::type RealType;
  typedef PointTempl<RealType> Point;

  explicit
  CompositeFunction () = default;

  /**
   * This class can be default move constructed and assigned.
   */
  CompositeFunction (CompositeFunction &&) = default;
  CompositeFunction & operator= (CompositeFunction &&) = default;

  /**
   * This class contains unique_ptr members so it can't be default
   * copied or assigned.
   */
  CompositeFunction (const CompositeFunction &) = delete;
  CompositeFunction & operator= (const CompositeFunction &) = delete;

  /**
   * The subfunctions vector is automatically cleaned up.
   */
  virtual ~CompositeFunction () = default;

  /**
   * Attach a new subfunction, along with a map from the indices of
   * the attached subfunction to the indices of the composed function.
   *
   * The composed function will return a vector whose value at index
   * \p index_map[i] is the value of the attached function at index i,
   * i.e.,
   * (*this)(x, t)(index_map[i]) will return f(x, t)(i).
   */
  void attach_subfunction (const FunctionBase<Output> & f,
                           const std::vector<unsigned int> & index_map)
  {
    const unsigned int subfunction_index =
      cast_int<unsigned int>(subfunctions.size());
    libmesh_assert_equal_to(subfunctions.size(), index_maps.size());

    subfunctions.push_back(f.clone());
    index_maps.push_back(index_map);

    unsigned int max_index =
      *std::max_element(index_map.begin(), index_map.end());

    if (max_index >= reverse_index_map.size())
      reverse_index_map.resize
        (max_index+1, std::make_pair(libMesh::invalid_uint,
                                     libMesh::invalid_uint));

    for (auto j : index_range(index_map))
      {
        libmesh_assert_less(index_map[j], reverse_index_map.size());
        libmesh_assert_equal_to(reverse_index_map[index_map[j]].first,
                                libMesh::invalid_uint);
        libmesh_assert_equal_to(reverse_index_map[index_map[j]].second,
                                libMesh::invalid_uint);
        reverse_index_map[index_map[j]] =
          std::make_pair(subfunction_index, j);
      }

    // Now check for time dependence
    // We only check the function we just added instead of researching all subfunctions
    // If this is the first subfunction, then that determines the time-dependence.
    if (subfunctions.size() == 1)
      this->_is_time_dependent = subfunctions[0]->is_time_dependent();

    // Otherwise, we have more than 1 function already.
    // If _is_time_dependent is true, then one of the previous
    // subfunctions is time-dependent and thus this CompositeFunction
    // time-dependent. If _is_time_dependent is false, then the subfunction
    // just added determines the time-dependence.
    else if (!this->_is_time_dependent)
      this->_is_time_dependent = (subfunctions.back())->is_time_dependent();
  }

  virtual Output operator() (const Point & p,
                             const Real time = 0) override
  {
    return this->component(0,p,time);
  }

  virtual void operator() (const Point & p,
                           const Real time,
                           DenseVector<Output> & output) override
  {
    libmesh_assert_greater_equal (output.size(),
                                  reverse_index_map.size());

    // Necessary in case we have output components not covered by
    // any subfunctions
    output.zero();

    DenseVector<Output> temp;
    for (auto i : index_range(subfunctions))
      {
        temp.resize(cast_int<unsigned int>(index_maps[i].size()));
        (*subfunctions[i])(p, time, temp);
        for (auto j : index_range(temp))
          output(index_maps[i][j]) = temp(j);
      }
  }

  virtual Output component (unsigned int i,
                            const Point & p,
                            Real time) override
  {
    if (i >= reverse_index_map.size() ||
        reverse_index_map[i].first == libMesh::invalid_uint)
      return 0;

    libmesh_assert_less(reverse_index_map[i].first,
                        subfunctions.size());
    libmesh_assert_not_equal_to(reverse_index_map[i].second,
                                libMesh::invalid_uint);
    return subfunctions[reverse_index_map[i].first]->
      component(reverse_index_map[i].second,p,time);
  }

  virtual std::unique_ptr<FunctionBase<Output>> clone() const override
  {
    CompositeFunction * returnval = new CompositeFunction();
    for (auto i : index_range(subfunctions))
      returnval->attach_subfunction(*subfunctions[i], index_maps[i]);
    return std::unique_ptr<FunctionBase<Output>> (returnval);
  }

  unsigned int n_subfunctions () const
  {
    return subfunctions.size();
  }

  unsigned int n_components () const
  {
    return reverse_index_map.size();
  }

private:
  // list of functions which fill in our values
  std::vector<std::unique_ptr<FunctionBase<Output>>> subfunctions;

  // for each function, list of which global indices it fills in
  std::vector<std::vector<unsigned int>> index_maps;

  // for each global index, which local index of which function is it?
  std::vector<std::pair<unsigned int, unsigned int>> reverse_index_map;
};


} // namespace libMesh

#endif // LIBMESH_COMPOSITE_FUNCTION_H
