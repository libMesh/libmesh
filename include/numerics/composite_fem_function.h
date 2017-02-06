// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_COMPOSITE_FEM_FUNCTION_H
#define LIBMESH_COMPOSITE_FEM_FUNCTION_H

// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/libmesh.h"
#include "libmesh/point.h"

// C++ includes
#include <algorithm>
#include <utility>
#include <vector>

namespace libMesh
{

/**
 * FEMFunction which is a function of another function.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief FEMFunction which is a function of another function.
 */
template <typename Output=Number>
class CompositeFEMFunction : public FEMFunctionBase<Output>
{
public:
  explicit
  CompositeFEMFunction () {}

  ~CompositeFEMFunction ()
  {
    for (std::size_t i=0; i != subfunctions.size(); ++i)
      delete subfunctions[i];
  }

  // Attach a new subfunction, along with a map from the indices of
  // that subfunction to the indices of the global function.
  // (*this)(index_map[i]) will return f(i).
  void attach_subfunction (const FEMFunctionBase<Output> & f,
                           const std::vector<unsigned int> & index_map)
  {
    const unsigned int subfunction_index = subfunctions.size();
    libmesh_assert_equal_to(subfunctions.size(), index_maps.size());

    subfunctions.push_back(f.clone().release());
    index_maps.push_back(index_map);

    unsigned int max_index =
      *std::max_element(index_map.begin(), index_map.end());

    if (max_index >= reverse_index_map.size())
      reverse_index_map.resize
        (max_index+1, std::make_pair(libMesh::invalid_uint,
                                     libMesh::invalid_uint));

    for (std::size_t j=0; j != index_map.size(); ++j)
      {
        libmesh_assert_less(index_map[j], reverse_index_map.size());
        libmesh_assert_equal_to(reverse_index_map[index_map[j]].first,
                                libMesh::invalid_uint);
        libmesh_assert_equal_to(reverse_index_map[index_map[j]].second,
                                libMesh::invalid_uint);
        reverse_index_map[index_map[j]] =
          std::make_pair(subfunction_index, j);
      }
  }

  virtual Output operator() (const FEMContext & c,
                             const Point & p,
                             const Real time = 0) libmesh_override
  {
    return this->component(c,0,p,time);
  }

  virtual void operator() (const FEMContext & c,
                           const Point & p,
                           const Real time,
                           DenseVector<Output> & output) libmesh_override
  {
    libmesh_assert_greater_equal (output.size(),
                                  reverse_index_map.size());

    // Necessary in case we have output components not covered by
    // any subfunctions
    output.zero();

    DenseVector<Output> temp;
    for (std::size_t i=0; i != subfunctions.size(); ++i)
      {
        temp.resize(index_maps[i].size());
        (*subfunctions[i])(c, p, time, temp);
        for (std::size_t j=0; j != temp.size(); ++j)
          output(index_maps[i][j]) = temp(j);
      }
  }

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output component (const FEMContext & c,
                            unsigned int i,
                            const Point & p,
                            Real time) libmesh_override
  {
    if (i >= reverse_index_map.size() ||
        reverse_index_map[i].first == libMesh::invalid_uint)
      return 0;

    libmesh_assert_less(reverse_index_map[i].first,
                        subfunctions.size());
    libmesh_assert_not_equal_to(reverse_index_map[i].second,
                                libMesh::invalid_uint);
    return subfunctions[reverse_index_map[i].first]->
      component(c, reverse_index_map[i].second, p, time);
  }

  virtual UniquePtr<FEMFunctionBase<Output> > clone() const libmesh_override
  {
    CompositeFEMFunction * returnval = new CompositeFEMFunction();
    for (std::size_t i=0; i != subfunctions.size(); ++i)
      returnval->attach_subfunction(*subfunctions[i], index_maps[i]);
    return UniquePtr<FEMFunctionBase<Output> > (returnval);
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
  std::vector<FEMFunctionBase<Output> *> subfunctions;

  // for each function, list of which global indices it fills in
  std::vector<std::vector<unsigned int> > index_maps;

  // for each global index, which local index of which function is it?
  std::vector<std::pair<unsigned int, unsigned int> > reverse_index_map;
};


} // namespace libMesh

#endif // LIBMESH_COMPOSITE_FEM_FUNCTION_H
